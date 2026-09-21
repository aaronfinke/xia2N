from __future__ import annotations

import logging
import os
import pathlib
import shutil
import subprocess
from dataclasses import dataclass, field

from dials.array_family import flex
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dxtbx.model import Crystal, ExperimentList
from dxtbx.serialize import load

from xia2.Driver.timing import record_step
from xia2.Handlers.Files import FileHandler
from xia2.Handlers.Streams import banner
from xia2.Modules.Laue_TOF.laue_tof import (
    NXLAUETOF,
    NXSNSEVENT,
    BinningParams,
    ExportParams,
    FileInput,
    IndexingParams,
    IntegrationParams,
    LaueTOFSetup,
    SpotfindingParams,
)

xia2_logger = logging.getLogger(__name__)

# Not a dials.index method: the label used when an orientation is supplied,
# which routes dials.index into the known orientation indexer instead.
KNOWN_ORIENTATION = "known_orientation"

# Foreground/background mask methods, as named by dials.tof_integrate
MASKS = ("ellipse", "seed_skewness")


@dataclass
class IndexingResult:
    """The outcome of one attempt at indexing an exposure."""

    method: str
    n_indexed: int = 0
    n_strong: int = 0
    rmsd_px: float | None = None
    crystal: Crystal | None = None
    accepted: bool = False
    reason: str = ""

    @property
    def fraction_indexed(self) -> float:
        return self.n_indexed / self.n_strong if self.n_strong else 0.0

    def summary(self) -> str:
        if not self.n_strong:
            return f"{self.method}: {self.reason}"
        rmsd = f"{self.rmsd_px:.2f} px" if self.rmsd_px is not None else "n/a"
        return (
            f"{self.method}: {self.n_indexed}/{self.n_strong} indexed"
            f" ({self.fraction_indexed:.1%}), RMSD {rmsd}"
            + ("" if self.accepted else f" - rejected, {self.reason}")
        )


@dataclass
class ExposureResult:
    """The outcome of processing one crystal orientation."""

    name: str
    directory: pathlib.Path
    image: str
    n_strong: int = 0
    indexing: IndexingResult | None = None
    integrated: bool = False
    attempts: list[IndexingResult] = field(default_factory=list)

    @property
    def indexed(self) -> bool:
        return self.indexing is not None and self.indexing.accepted


@dataclass
class IntegrationChoices:
    """
    Integration decisions taken on the first orientation and reused for the rest.

    Every orientation of a run has to be integrated the same way for the
    batches to be comparable when they are scaled together, and retrying a
    method that has already failed once, or re-comparing the two masks on every
    orientation, only costs time.
    """

    method: str | None = None
    mask: str | None = None


def _executable(name: str, override: pathlib.Path | None = None) -> str:
    """Find a program, either where the user said it is or on $PATH."""
    if override:
        return os.fspath(override)
    found = shutil.which(name)
    if not found:
        raise ValueError(
            f"Unable to find {name} on $PATH. Either add it to $PATH or give its"
            f" location in the xia2.laue_tof phil scope."
        )
    return found


def _run_program(
    command: list[str], working_directory: pathlib.Path, name: str
) -> subprocess.CompletedProcess:
    """Run one of the external programs, raising on anything but success."""
    xia2_logger.debug("Running %s", " ".join(command))
    with record_step(name):
        result = subprocess.run(
            command, cwd=working_directory, capture_output=True, encoding="utf-8"
        )
    if result.returncode:
        raise ValueError(
            f"{name} returned error status:\n{result.stderr or result.stdout}"
        )
    return result


def _record_log(tag: str, logfile: pathlib.Path) -> None:
    if logfile.is_file():
        FileHandler.record_log_file(tag, str(logfile))


def bin_run(
    working_directory: pathlib.Path,
    image: str,
    file_class: str,
    nbins: int,
    params: BinningParams,
) -> pathlib.Path:
    """
    Reduce a raw data file to TOF bins with essnmx.

    The reduced file is written into the working directory, and reused if it is
    already there, as the reduction is the most expensive step of the pipeline.
    """
    output_file = working_directory / f"binned_{nbins}.h5"
    if output_file.is_file():
        xia2_logger.info(f"Using existing reduction {output_file}")
        return output_file

    if file_class == NXSNSEVENT:
        program = _executable("essmandi-reduce", params.essmandi_reduce)
    else:
        program = _executable("essnmx-reduce", params.essnmx_reduce)
    name = pathlib.Path(program).name

    command = [
        program,
        "--input-file",
        image,
        "--nbins",
        str(nbins),
        "--output-file",
        os.fspath(output_file),
        "--overwrite",
    ]
    if params.min_time_bin is not None:
        command += ["--min-time-bin", str(params.min_time_bin)]
    if params.max_time_bin is not None:
        command += ["--max-time-bin", str(params.max_time_bin)]
    if params.min_time_bin is not None or params.max_time_bin is not None:
        command += ["--time-bin-unit", params.time_bin_unit]
    if params.detector_ids:
        command += ["--detector-ids"] + [str(i) for i in params.detector_ids]
    command += params.extra_args

    xia2_logger.notice(banner(f"Reducing to {nbins} TOF bins"))  # type: ignore
    _run_program(command, working_directory, name)
    if not output_file.is_file():
        raise ValueError(f"{name} did not write the expected output {output_file}")
    return output_file


def run_import(
    working_directory: pathlib.Path,
    image: str,
    file_input: FileInput,
    path_type: str = "image",
) -> pathlib.Path:
    """Run dials.import on one exposure, and check that the data really is TOF."""
    experiments = working_directory / "imported.expt"
    command = [
        _executable("dials.import"),
        f"{path_type}={image}" if path_type != "image" else image,
        "output.experiments=imported.expt",
    ]
    if file_input.import_phil:
        command.insert(1, os.fspath(file_input.import_phil))
    if file_input.mask:
        command.append(f"mask={os.fspath(file_input.mask)}")

    xia2_logger.notice(banner("Importing"))  # type: ignore
    _run_program(command, working_directory, "dials.import")
    _record_log(
        f"{working_directory.name} import", working_directory / "dials.import.log"
    )

    expts = load.experiment_list(experiments, check_format=False)
    if not expts.all_tof():
        raise ValueError(
            f"{image} was not imported as time-of-flight data. Check that the"
            " file holds TOF bins, and that it is being read by the expected"
            " format class (see dials.show imported.expt)."
        )
    return experiments


def find_spots(working_directory: pathlib.Path, params: SpotfindingParams) -> int:
    """Run dials.find_spots, returning the number of strong spots found."""
    command = [
        _executable("dials.find_spots"),
        "imported.expt",
        "output.reflections=strong.refl",
        f"spotfinder.filter.min_spot_size={params.min_spot_size}",
        f"spotfinder.threshold.algorithm={params.threshold_algorithm}",
        f"spotfinder.mp.nproc={params.nproc}",
    ]
    if params.phil:
        command.insert(1, os.fspath(params.phil))
    if params.max_spot_size:
        command.append(f"spotfinder.filter.max_spot_size={params.max_spot_size}")
    if params.d_min:
        command.append(f"spotfinder.filter.d_min={params.d_min}")
    if params.d_max:
        command.append(f"spotfinder.filter.d_max={params.d_max}")

    xia2_logger.notice(banner("Spotfinding"))  # type: ignore
    _run_program(command, working_directory, "dials.find_spots")
    _record_log(
        f"{working_directory.name} find_spots",
        working_directory / "dials.find_spots.log",
    )

    strong = flex.reflection_table.from_file(working_directory / "strong.refl")
    n_strong = strong.size()
    xia2_logger.info(f"Found {n_strong} strong spots")
    if n_strong > params.max_strong:
        raise ValueError(
            f"{n_strong} strong spots found, more than spotfinding.max_strong"
            f" ({params.max_strong}). Large spot lists trip an assertion in the"
            " DIALS TOF reflection predictor during refinement; raise"
            " spotfinding.min_spot_size, or raise the limit if you are sure."
        )
    return n_strong


def _primitive_setting_a_matrix(crystal: Crystal) -> tuple[float, ...]:
    """
    The A matrix of a crystal, in the primitive setting. Needed currently
    due to a DIALS bug.

    dials.index expects known_symmetry.A_matrix in the primitive setting and
    contains the conversion, but IndexerKnownOrientation assigns it to the loop
    variable and discards it (DIALS 3.30), while the post-indexing conversion
    back to the centred setting still runs - so a centred cell comes back wrong.
    Converting here supplies the setting the indexer needs, leaving the skipped
    conversion nothing to do. Note this is tied to the defect: the crystal is
    built with known_symmetry.space_group, so a fixed DIALS would convert this
    matrix a second time. _assess_indexing's cell check is what would catch it.
    """
    cb_op = crystal.get_space_group().info().change_of_basis_op_to_primitive_setting()
    return crystal.change_basis(cb_op).get_A()


def _assess_indexing(
    working_directory: pathlib.Path, method: str, params: IndexingParams, n_strong: int
) -> IndexingResult:
    """
    Decide whether an indexing solution is good enough to keep.

    Indexing that 'worked' is not enough: a small, clean spot list indexes
    almost completely, and an alternative indexing can index plenty of
    reflections on the wrong cell. So check the cell, the number indexed and
    the positional residuals.
    """
    result = IndexingResult(method=method, n_strong=n_strong)
    experiments = working_directory / "indexed.expt"
    reflections = working_directory / "indexed.refl"
    if not (experiments.is_file() and reflections.is_file()):
        result.reason = "no solution found"
        return result

    expts = load.experiment_list(experiments, check_format=False)
    refl = flex.reflection_table.from_file(reflections)
    indexed = refl.select(refl.get_flags(refl.flags.indexed))
    result.n_indexed = indexed.size()

    if "xyzcal.px" in indexed and indexed.size():
        dx, dy, _ = (indexed["xyzobs.px.value"] - indexed["xyzcal.px"]).parts()
        result.rmsd_px = float(flex.mean(dx * dx + dy * dy) ** 0.5)

    if not expts.crystals():
        result.reason = "no crystal model"
        return result
    crystal = expts.crystals()[0]
    result.crystal = crystal

    if result.n_indexed < params.min_indexed:
        result.reason = (
            f"only {result.n_indexed} reflections indexed, fewer than"
            f" indexing.min_indexed ({params.min_indexed})"
        )
        return result

    if params.unit_cell and not crystal.get_unit_cell().is_similar_to(
        params.unit_cell, relative_length_tolerance=0.1, absolute_angle_tolerance=5.0
    ):
        result.reason = (
            f"refined cell {crystal.get_unit_cell()} does not match the"
            f" target cell {params.unit_cell}"
        )
        return result

    result.accepted = True
    return result


def index(
    working_directory: pathlib.Path,
    params: IndexingParams,
    n_strong: int,
    seed_crystal: Crystal | None = None,
) -> list[IndexingResult]:
    """
    Run dials.index, trying each method in turn until a solution is accepted.

    Near the limit of what a detector layout can support, indexing succeeds
    stochastically rather than reproducibly, so one recipe is not enough.
    Returns every attempt made, the last of which is the accepted solution if
    there is one.
    """
    if n_strong < params.min_spots:
        return [
            IndexingResult(
                method="none",
                n_strong=n_strong,
                reason=f"only {n_strong} strong spots, fewer than indexing.min_spots",
            )
        ]

    base_command = [
        _executable("dials.index"),
        "imported.expt",
        "strong.refl",
    ]
    if params.phil:
        base_command.insert(1, os.fspath(params.phil))
    if params.space_group:
        base_command.append(f"indexing.known_symmetry.space_group={params.space_group}")
    if params.unit_cell:
        cell = ",".join(f"{p}" for p in params.unit_cell.parameters())
        base_command.append(f"indexing.known_symmetry.unit_cell={cell}")
    if params.max_cell:
        base_command.append(f"indexing.max_cell={params.max_cell}")
    if params.d_min_start:
        base_command.append(
            f"indexing.refinement_protocol.d_min_start={params.d_min_start}"
        )
    base_command += [
        f"refinement.reflections.outlier.algorithm={params.outlier_algorithm}",
        f"refinement.reflections.outlier.tukey.iqr_multiplier={params.outlier_iqr_multiplier}",
        f"refinement.reflections.outlier.separate_panels={params.outlier_separate_panels}",
    ]
    if seed_crystal:
        a_matrix = ",".join(f"{v}" for v in _primitive_setting_a_matrix(seed_crystal))
        base_command.append(f"indexing.known_symmetry.A_matrix={a_matrix}")

    xia2_logger.notice(banner("Indexing"))  # type: ignore

    if seed_crystal:
        # Supplying an A matrix makes dials.index use the known orientation
        # indexer, which ignores indexing.method, so there is no ladder to walk.
        xia2_logger.info("Seeding with the crystal model of a previous exposure")
        methods = [KNOWN_ORIENTATION]
    else:
        methods = list(params.methods)

    attempts = []
    for method in methods:
        command = base_command + [f"output.log=dials.index.{method}.log"]
        if method != KNOWN_ORIENTATION:
            command.append(f"indexing.method={method}")
        try:
            _run_program(command, working_directory, "dials.index")
        except ValueError as e:
            xia2_logger.debug(f"Indexing with method={method} failed:\n{e}")
            result = IndexingResult(
                method=method, n_strong=n_strong, reason="dials.index failed"
            )
        else:
            result = _assess_indexing(working_directory, method, params, n_strong)
        _record_log(
            f"{working_directory.name} index ({method})",
            working_directory / f"dials.index.{method}.log",
        )
        xia2_logger.info(result.summary())
        attempts.append(result)
        if result.accepted:
            break
    return attempts


def refine(working_directory: pathlib.Path) -> None:
    """Run dials.refine on the indexing solution."""
    command = [
        _executable("dials.refine"),
        "indexed.expt",
        "indexed.refl",
    ]
    xia2_logger.notice(banner("Refining"))  # type: ignore
    _run_program(command, working_directory, "dials.refine")
    _record_log(
        f"{working_directory.name} refine", working_directory / "dials.refine.log"
    )


def _tof_integrate_command(
    params: IntegrationParams,
    method: str,
    mask: str,
    tag: str,
    experiments: str,
    reflections: str,
) -> list[str]:
    """The dials.tof_integrate command line for one method/mask combination."""
    command = [
        _executable("dials.tof_integrate"),
        experiments,
        reflections,
        f"method={method}",
        f"mask={mask}",
        f"integration_type={params.integration_type}",
        f"background_model={params.background_model}",
        f"ellipse_mask.scale={params.ellipse_mask_scale}",
        f"bbox_tof_padding={params.bbox_tof_padding}",
        f"bbox_xy_padding={params.bbox_xy_padding}",
        f"corrections.lorentz={params.lorentz}",
        f"mp.nproc={params.nproc}",
        f"output.experiments=integrated_{tag}.expt",
        f"output.reflections=integrated_{tag}.refl",
        f"output.log=tof_integrate_{tag}.log",
    ]
    if params.phil:
        command.insert(1, os.fspath(params.phil))
    if params.calculated_d_min:
        command.append(f"calculated.dmin={params.calculated_d_min}")
    if params.wavelength_range:
        low, high = params.wavelength_range
        command.append(f"wavelength_range={low},{high}")
    if params.incident_run:
        command.append(f"corrections.incident_run={os.fspath(params.incident_run)}")
    if params.empty_run:
        command.append(f"corrections.empty_run={os.fspath(params.empty_run)}")
    if params.absorption:
        for name, value in params.absorption.items():
            command.append(f"corrections.absorption.target_spectrum.{name}={value}")
    return command


def _integration_quality(reflections: pathlib.Path) -> tuple[int, float]:
    """
    How well an integration run went: how many reflections were integrated by
    summation, and their mean I/sigma.

    Summation is used for the comparison whatever the method, because it is the
    one set of intensities that every run produces.
    """
    table = flex.reflection_table.from_file(reflections)
    table = table.select(table.get_flags(table.flags.integrated_sum))
    if not table.size():
        return 0, 0.0
    variance = table["intensity.sum.variance"]
    table = table.select(variance > 0)
    if not table.size():
        return 0, 0.0
    i_over_sigma = table["intensity.sum.value"] / flex.sqrt(
        table["intensity.sum.variance"]
    )
    return table.size(), float(flex.mean(i_over_sigma))


def _profile_fraction(reflections: pathlib.Path) -> float:
    """
    The fraction of integrated reflections that a profile fit succeeded for.

    A profile run can finish cleanly having fitted almost nothing, leaving an
    intensity.prf column that covers a handful of reflections, so success of
    the program is not success of the method.
    """
    table = flex.reflection_table.from_file(reflections)
    if "intensity.prf.value" not in table:
        return 0.0
    n_summed = table.get_flags(table.flags.integrated_sum).count(True)
    if not n_summed:
        return 0.0
    return table.get_flags(table.flags.integrated_prf).count(True) / n_summed


def _run_tof_integrate(
    working_directory: pathlib.Path,
    params: IntegrationParams,
    method: str,
    mask: str,
    experiments: str,
    reflections: str,
) -> tuple[str, pathlib.Path]:
    """
    Integrate with one mask, falling back to summation if a profile method
    fails outright.

    Profile fitting is not a per-reflection risk only: a single shoebox that
    the fit cannot bound aborts dials.tof_integrate, so the whole run is lost.
    Returns the method that actually produced the output.
    """
    tag = f"{method}_{mask}"
    try:
        _run_program(
            _tof_integrate_command(params, method, mask, tag, experiments, reflections),
            working_directory,
            "dials.tof_integrate",
        )
    except ValueError as e:
        fallback = params.fallback_method
        if method == fallback or fallback == "none":
            raise
        xia2_logger.warning(
            f"Integration with method={method} failed, falling back to"
            f" {fallback}. The failure was:\n{str(e).strip().splitlines()[-1]}"
        )
        _record_log(
            f"{working_directory.name} integrate ({tag}, failed)",
            working_directory / f"tof_integrate_{tag}.log",
        )
        return _run_tof_integrate(
            working_directory, params, fallback, mask, experiments, reflections
        )

    _record_log(
        f"{working_directory.name} integrate ({tag})",
        working_directory / f"tof_integrate_{tag}.log",
    )
    integrated = working_directory / f"integrated_{tag}.refl"
    if not integrated.is_file():
        raise ValueError(f"dials.tof_integrate did not write {integrated.name}")

    if method != "summation":
        fraction = _profile_fraction(integrated)
        if fraction < params.min_profile_fraction:
            # The summation intensities in this file are still good, so there
            # is nothing to rerun - but say so, and stop paying for a profile
            # fit that is not delivering on the orientations still to come.
            xia2_logger.warning(
                f"{method} fitted only {fraction:.1%} of the reflections, fewer"
                f" than integration.min_profile_fraction"
                f" ({params.min_profile_fraction:.0%}), so the summation"
                " intensities will be used"
            )
            return "summation", integrated
        xia2_logger.info(f"{method} fitted {fraction:.1%} of the reflections")
    return method, integrated


def _record_lorentz(reflections: pathlib.Path, applied: bool) -> None:
    """
    Record whether the Lorentz correction was applied, per reflection.

    dials.tof_integrate applies the correction to the intensities themselves
    and leaves no trace of it in the output beyond a line in its log, so the
    file cannot be asked later whether it has been corrected. The Lorentz
    factor must be applied exactly once, and the scaling programs downstream
    have their own keyword for it, so carry the answer with the data.
    """
    table = flex.reflection_table.from_file(reflections)
    table["lorentz_applied"] = flex.bool(table.size(), applied)
    table.as_file(reflections)


def integrate(
    working_directory: pathlib.Path,
    params: IntegrationParams,
    choices: IntegrationChoices | None = None,
    experiments: str = "refined.expt",
    reflections: str = "refined.refl",
) -> pathlib.Path:
    """
    Run dials.tof_integrate, returning the integrated reflection file.

    With mask=both the data are integrated with each foreground mask and the
    better of the two is kept; that choice, and any fallback from a failed
    profile method, are recorded in choices and reused for the remaining
    orientations, so that every batch is integrated the same way.
    """
    choices = choices if choices is not None else IntegrationChoices()
    method = choices.method or params.method
    if choices.mask:
        masks = [choices.mask]
    elif params.mask == "both":
        masks = list(MASKS)
    else:
        masks = [params.mask]

    xia2_logger.notice(banner("Integrating"))  # type: ignore

    integrated_files = {}
    for mask in masks:
        method, integrated = _run_tof_integrate(
            working_directory, params, method, mask, experiments, reflections
        )
        integrated_files[mask] = integrated
    choices.method = method

    if len(integrated_files) > 1:
        quality = {
            mask: _integration_quality(f) for mask, f in integrated_files.items()
        }
        for mask, (n_integrated, i_over_sigma) in quality.items():
            xia2_logger.info(
                f"{mask}: {n_integrated} reflections integrated,"
                f" <I/sigma> {i_over_sigma:.1f}"
            )
        index = 1 if params.mask_metric == "i_over_sigma" else 0
        mask = max(quality, key=lambda m: quality[m][index])
        xia2_logger.info(f"Keeping the {mask} mask, by {params.mask_metric}")
    else:
        mask = masks[0]
    choices.mask = mask

    integrated = working_directory / "integrated.refl"
    shutil.copyfile(integrated_files[mask], integrated)
    shutil.copyfile(
        integrated_files[mask].with_suffix(".expt"),
        working_directory / "integrated.expt",
    )
    _record_lorentz(integrated, params.lorentz)
    FileHandler.record_data_file(str(working_directory / "integrated.expt"))
    FileHandler.record_data_file(str(integrated))
    return integrated


def combine(
    working_directory: pathlib.Path, results: list[ExposureResult]
) -> tuple[pathlib.Path, pathlib.Path]:
    """
    Concatenate the integrated output of every orientation into one pair of
    files, ready to export.

    Each orientation keeps its own experiment identifier and its own batch, so
    that the exported intensities can be scaled per setting downstream.
    """
    integrated = [result for result in results if result.integrated]
    if not integrated:
        raise ValueError(
            "No orientation was integrated, so there is nothing to combine"
        )
    working_directory.mkdir(parents=True, exist_ok=True)
    xia2_logger.notice(banner("Combining"))  # type: ignore

    experiments = ExperimentList()
    tables = []
    # Which orientation each of the combined experiments came from
    origins: list[str] = []
    for result in integrated:
        expts = load.experiment_list(
            result.directory / "integrated.expt", check_format=False
        )
        table = flex.reflection_table.from_file(result.directory / "integrated.refl")
        split = parse_multiple_datasets([table])
        experiments.extend(expts)
        tables.extend(split)
        origins.extend([result.name] * len(split))
    experiments, tables = assign_unique_identifiers(experiments, tables)

    # dials.export numbers the Laue batches by imageset_id, which is zero in
    # every per-orientation file. Renumber it to the position of the experiment
    # in the combined list, which is also where its imageset now sits, so that
    # each orientation exports as its own batch.
    for i, (table, origin) in enumerate(zip(tables, origins)):
        table["imageset_id"] = flex.int(table.size(), i)
        xia2_logger.info(f"Batch {i}: {origin}, {table.size()} reflections")

    if len(tables) > 1:
        combined = flex.reflection_table.concat(tables)
    else:
        combined = tables[0]
    _report_lorentz(combined)

    experiments_file = working_directory / "combined.expt"
    reflections_file = working_directory / "combined.refl"
    experiments.as_file(experiments_file)
    combined.as_file(reflections_file)
    xia2_logger.info(
        f"Combined {len(experiments)} experiment(s), {combined.size()} reflections,"
        f" into {experiments_file.name} and {reflections_file.name}"
    )
    FileHandler.record_data_file(str(experiments_file))
    FileHandler.record_data_file(str(reflections_file))
    return experiments_file, reflections_file


def _intensity_choice(reflections: pathlib.Path, params: ExportParams) -> str:
    """
    Which intensities to export.

    dials.export refuses its own auto choice for SHELX when both summation and
    profile-fitted intensities are present, which is exactly what a successful
    profile run produces, so the choice is always made here.
    """
    if params.intensity != "auto":
        return params.intensity
    fraction = _profile_fraction(reflections)
    if fraction >= params.min_profile_fraction:
        return "profile"
    if fraction:
        xia2_logger.info(
            f"Profile-fitted intensities cover only {fraction:.1%} of the"
            " reflections, so the summation intensities are exported"
        )
    return "sum"


def _report_lorentz(table: flex.reflection_table) -> None:
    """
    Say whether the combined data carry the Lorentz correction.

    The correction has to be applied exactly once, here or by the scaling
    program, so the answer decides what is passed downstream. Orientations that
    disagree cannot be scaled together.
    """
    if "lorentz_applied" not in table:
        xia2_logger.warning(
            "The combined data do not record whether the Lorentz correction was"
            " applied. It must be applied exactly once - check before scaling."
        )
        return
    applied = set(table["lorentz_applied"])
    if len(applied) > 1:
        raise ValueError(
            "Some orientations were integrated with the Lorentz correction and"
            " some without, so they cannot be scaled together. Reintegrate them"
            " the same way (integration.lorentz)."
        )
    if applied.pop():
        xia2_logger.info(
            "The Lorentz correction was applied during integration, so it must"
            " not be applied again when scaling."
        )
    else:
        xia2_logger.info(
            "The Lorentz correction was not applied during integration, so it"
            " must be applied when scaling."
        )


def export_shelx(
    working_directory: pathlib.Path,
    params: ExportParams,
    experiments: str = "combined.expt",
    reflections: str = "combined.refl",
) -> pathlib.Path:
    """Run dials.export format=shelx, returning the exported hkl file."""
    intensity = _intensity_choice(working_directory / reflections, params)
    command = [
        _executable("dials.export"),
        experiments,
        reflections,
        "format=shelx",
        f"intensity={intensity}",
        f"shelx.composition={params.composition}",
    ]
    if params.phil:
        command.insert(1, os.fspath(params.phil))

    xia2_logger.notice(banner("Exporting"))  # type: ignore
    _run_program(command, working_directory, "dials.export")
    _record_log("export", working_directory / "dials.export.log")

    hklout = working_directory / "dials.hkl"
    if not hklout.is_file():
        raise ValueError("dials.export did not write dials.hkl")
    xia2_logger.info(f"Exported the {intensity} intensities to {hklout.name}")
    FileHandler.record_data_file(str(hklout))
    ins = working_directory / "dials.ins"
    if ins.is_file():
        FileHandler.record_data_file(str(ins))
    return hklout


def _prepare_images(
    working_directory: pathlib.Path,
    image: str,
    file_class: str,
    setup: LaueTOFSetup,
) -> tuple[str, str]:
    """
    Get the files to index and to integrate from one input file.

    These are the same file unless the input needs reducing and two different
    binnings have been asked for, in which case the coarser reduction is
    indexed and the finer one integrated.
    """
    binning = setup.binning_params
    if "bin" not in setup.options.steps or image not in setup.files_to_bin:
        return image, image
    index_image = os.fspath(
        bin_run(working_directory, image, file_class, binning.index_bins, binning)
    )
    if binning.single_reduction:
        return index_image, index_image
    integrate_image = os.fspath(
        bin_run(working_directory, image, file_class, binning.integrate_bins, binning)
    )
    return index_image, integrate_image


def _refined_crystal(working_directory: pathlib.Path) -> Crystal:
    """The crystal model of a refined experiment."""
    expts = load.experiment_list(working_directory / "refined.expt", check_format=False)
    return expts.crystals()[0]


def _index_and_refine(
    working_directory: pathlib.Path,
    image: str,
    setup: LaueTOFSetup,
    result: ExposureResult,
    seed_crystal: Crystal | None,
    path_type: str = "image",
) -> bool:
    """
    Import, find spots, index and refine in one directory.

    Returns whether there is a refined solution to integrate.
    """
    run_import(working_directory, image, setup.file_input, path_type)

    if "find_spots" not in setup.options.steps:
        return False
    result.n_strong = find_spots(working_directory, setup.spotfinding_params)

    if "index" not in setup.options.steps:
        return False
    attempts = index(
        working_directory, setup.indexing_params, result.n_strong, seed_crystal
    )
    result.attempts = attempts
    result.indexing = attempts[-1] if attempts else None
    if not result.indexed:
        reason = result.indexing.reason if result.indexing else "no attempt made"
        xia2_logger.warning(f"{result.name} did not index: {reason}")
        return False

    if "refine" not in setup.options.steps:
        return False
    refine(working_directory)
    return True


def process_exposure(
    working_directory: pathlib.Path,
    image: str,
    setup: LaueTOFSetup,
    seed_crystal: Crystal | None = None,
    path_type: str = "image",
    choices: IntegrationChoices | None = None,
) -> ExposureResult:
    """
    Process one crystal orientation: bin, import, find spots, index, refine
    and integrate, in its own directory.

    If the indexing and integration binnings differ, the orientation is
    determined on the coarser reduction and then carried into the finer one,
    which is processed in a subdirectory.
    """
    working_directory.mkdir(parents=True, exist_ok=True)
    name = working_directory.name
    xia2_logger.notice(banner(f"Processing {name}"))  # type: ignore
    xia2_logger.info(f"Input: {image}")

    file_class = setup.input_classes.get(image, NXLAUETOF)
    result = ExposureResult(name=name, directory=working_directory, image=image)

    index_image, integrate_image = _prepare_images(
        working_directory, image, file_class, setup
    )

    if not _index_and_refine(
        working_directory, index_image, setup, result, seed_crystal, path_type
    ):
        return result

    if "integrate" not in setup.options.steps:
        return result

    if integrate_image == index_image:
        integrate(working_directory, setup.integration_params, choices)
        result.integrated = True
        return result

    # A finer reduction was asked for to integrate against. It is a different
    # file, with its own imageset, so it has to be processed in its own right -
    # but seeded with the orientation just refined, rather than indexed afresh.
    xia2_logger.info(
        f"Carrying the refined orientation into the"
        f" {setup.binning_params.integrate_bins} bin reduction"
    )
    integrate_directory = (
        working_directory / f"bins_{setup.binning_params.integrate_bins}"
    )
    integrate_directory.mkdir(parents=True, exist_ok=True)
    integrate_result = ExposureResult(
        name=f"{name}/{integrate_directory.name}",
        directory=integrate_directory,
        image=integrate_image,
    )
    if not _index_and_refine(
        integrate_directory,
        integrate_image,
        setup,
        integrate_result,
        _refined_crystal(working_directory),
    ):
        xia2_logger.warning(
            f"{name}: the refined orientation could not be transferred to the"
            " integration binning, so this orientation was not integrated"
        )
        return result

    integrate(integrate_directory, setup.integration_params, choices)
    result.integrated = True
    result.directory = integrate_directory
    return result


def _report(results: list[ExposureResult]) -> None:
    """Summarise what happened to each orientation."""
    xia2_logger.notice(banner("Summary"))  # type: ignore
    for result in results:
        if result.integrated:
            state = "integrated"
        elif result.indexed:
            state = "indexed, not integrated"
        elif result.n_strong:
            state = "not indexed"
        else:
            state = "no strong spots"
        line = f"{result.name}: {state}"
        if result.indexing and result.indexing.accepted:
            line += f" ({result.indexing.summary()})"
        xia2_logger.info(line)


def run_data_integration(
    root_working_directory: pathlib.Path, setup: LaueTOFSetup
) -> list[ExposureResult]:
    """
    Process every crystal orientation, each in its own directory.

    All the exposures are of the same crystal, so once one of them has indexed
    its crystal model is used to seed the rest. Any exposure that did not index
    on its own is then retried with that seed, which is what recovers the
    orientations that never index de novo.
    """
    inputs = [("image", image) for image in setup.file_input.images]
    inputs += [("template", t) for t in setup.file_input.templates]
    inputs += [("directory", d) for d in setup.file_input.directories]

    results: list[ExposureResult] = []
    seed_crystal: Crystal | None = None
    choices = IntegrationChoices()

    for i, (path_type, image) in enumerate(inputs, start=1):
        working_directory = root_working_directory / f"orientation_{i}"
        result = process_exposure(
            working_directory, image, setup, seed_crystal, path_type, choices
        )
        results.append(result)
        if (
            seed_crystal is None
            and setup.indexing_params.seed_from_first
            and result.indexed
        ):
            seed_crystal = result.indexing.crystal  # type: ignore[union-attr]

    # Retry anything that did not index, now that there is a model to seed with
    if seed_crystal and setup.indexing_params.seed_from_first:
        for i, result in enumerate(results):
            if result.indexed or not result.n_strong:
                continue
            xia2_logger.info(f"Retrying {result.name} with the seed crystal")
            retried = process_exposure(
                result.directory,
                result.image,
                setup,
                seed_crystal,
                inputs[i][0],
                choices,
            )
            results[i] = retried

    _report(results)

    if "combine" not in setup.options.steps:
        return results
    if not any(result.integrated for result in results):
        xia2_logger.warning(
            "No orientation was integrated, so there is nothing to combine or export"
        )
        return results

    scale_directory = root_working_directory / "scale"
    combine(scale_directory, results)
    if "export" in setup.options.steps:
        export_shelx(scale_directory, setup.export_params)
    return results
