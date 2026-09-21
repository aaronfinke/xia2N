from __future__ import annotations

import errno
import logging
import os
import pathlib
from dataclasses import dataclass, field

import h5py
import iotbx.phil
from cctbx import crystal, sgtbx, uctbx
from dials.util.system import CPU_COUNT
from libtbx import Auto, AutoType

from xia2.Modules.Laue_TOF.util import report_timing

xia2_logger = logging.getLogger(__name__)

# Input file classes, as identified from the NeXus entry. NXlauetof files are
# already TOF-binned and can go straight to dials.import; the others hold raw
# event data and must be reduced first (see BinningParams).
NXLAUETOF = "nxlauetof"
TOFRAW = "tofraw"
NXSNSEVENT = "nxsnsevent"
UNKNOWN = "unknown"

# Integration methods, as named by dials.tof_integrate
INTEGRATION_METHODS = ("summation", "profile_1d_ibix", "profile_3d_gutmann")
# Indexing methods that work with polychromatic (Laue-TOF) data
INDEXING_METHODS = (
    "fft3d",
    "fft1d",
    "real_space_grid_search",
    "low_res_spot_match",
    "pink_indexer",
)

phil_str = """
image = None
  .type = str
  .multiple = True
  .help = "Path to an image file. Give one per crystal orientation (setting)."
  .expert_level = 0
template = None
  .type = str
  .multiple = True
  .help = "The image sequence template"
  .expert_level = 0
directory = None
  .type = str
  .multiple = True
  .help = "A directory with images"
  .expert_level = 0
vanadium_run = None
  .type = path
  .help = "Incident run (e.g. vanadium), used to normalise the intensities by the"
          "incident spectrum during integration. Optional, but recommended."
  .expert_level = 0
empty_run = None
  .type = path
  .help = "Empty instrument run, used to correct the incident run for background."
          "Only used if vanadium_run is also given."
  .expert_level = 0
mask = None
  .type = path
  .help = "A mask to use for spotfinding and integration"
  .expert_level = 1
space_group = None
  .type = space_group
  .help = "Space group to be used for indexing and integration."
  .expert_level = 0
d_min = None
  .type = float
  .help = "High resolution cutoff for spotfinding."
  .expert_level = 1
d_max = None
  .type = float
  .help = "Low resolution cutoff for spotfinding."
  .expert_level = 1
multiprocessing {
  nproc = Auto
    .type = int
    .expert_level = 2
}
dials_import.phil = None
  .type = path
  .help = "Phil file to use for dials.import. Parameters defined in the"
          "xia2.laue_tof phil scope will take precedent over identical options"
          "defined in the phil file."
  .expert_level = 3
"""

binning_phil_str = """
binning {
  enabled = Auto
    .type = bool
    .help = "Reduce the raw data to TOF bins with essnmx before importing."
            "Auto means bin only those input files that are not already binned"
            "(i.e. that are not NXlauetof files). Set to False to use the input"
            "files exactly as given."
    .expert_level = 1
  index_bins = 50
    .type = int(value_min=1)
    .help = "Number of TOF bins in the reduction used for spotfinding and"
            "indexing. Coarser binning indexes more reliably."
    .expert_level = 1
  integrate_bins = 200
    .type = int(value_min=1)
    .help = "Number of TOF bins in the reduction used for integration. Finer"
            "binning resolves more reflections. If equal to index_bins, a single"
            "reduction is made and used for both."
    .expert_level = 1
  min_time_bin = None
    .type = float
    .help = "Lower edge of the time range to bin over, in time_bin_unit."
    .expert_level = 2
  max_time_bin = None
    .type = float
    .help = "Upper edge of the time range to bin over, in time_bin_unit."
    .expert_level = 2
  time_bin_unit = ms *us ns
    .type = choice
    .expert_level = 2
  detector_ids = None
    .type = ints
    .help = "Detector panels to reduce. Defaults to all panels."
    .expert_level = 2
  essnmx_reduce = None
    .type = path
    .help = "Path to the essnmx-reduce executable, used for NMX data."
            "Defaults to looking it up on $PATH."
    .expert_level = 3
  tof_padding = 100
    .type = float
    .help = "Padding (microseconds) added to each end of the time-of-flight"
            "range when histogramming event data in place, as dxtbx does for"
            "MANDI. This is the dxtbx default. The detector geometry that goes"
            "with it, panel size included, is read from the format class, not"
            "set here."
    .expert_level = 2
  extra_args = None
    .type = str
    .multiple = True
    .help = "Additional command line arguments, passed through to the reduction"
            "program unchanged e.g. extra_args=--chunk-size-events=100000"
    .expert_level = 3
}
"""

spotfinding_phil_str = """
spotfinding {
  min_spot_size = 15
    .type = int
    .help = "The minimum spot size to allow in spotfinding. Note that this is"
            "deliberately larger than the DIALS TOF default of 6, which admits"
            "small split peaks; raising it markedly improves the positional"
            "residuals of the indexing solution."
    .expert_level = 1
  max_spot_size = None
    .type = int
    .help = "The maximum spot size to allow in spotfinding."
    .expert_level = 2
  threshold_algorithm = dispersion *dispersion_extended radial_profile
    .type = choice
    .expert_level = 2
  max_strong = 35000
    .type = int
    .help = "If more than this many strong spots are found, stop rather than"
            "proceeding to indexing. Large spot lists trip an assertion in the"
            "DIALS TOF reflection predictor during refinement."
    .expert_level = 3
  phil = None
    .type = path
    .help = "Phil options file to use for spotfinding. Parameters defined in"
            "the xia2.laue_tof phil scope will take precedent over identical"
            "options defined in the phil file."
    .expert_level = 3
}
"""

indexing_phil_str = """
indexing {
  unit_cell = None
    .type = unit_cell
    .expert_level = 0
  method = None
    .type = str
    .help = "Indexing method to use. If None, each method given in"
            "indexing.ladder is tried in turn until one succeeds."
    .expert_level = 1
  ladder = fft3d real_space_grid_search fft1d
    .type = strings
    .help = "Indexing methods to try, in order, when indexing.method is not set."
            "fft1d is the method that succeeds on sparse, narrow-wedge data from"
            "distant detector layouts, where fft3d cannot populate its grid."
    .expert_level = 2
  max_cell = None
    .type = float
    .help = "Fixed maximum cell length (Angstrom) for the basis vector search."
            "If None, it is set to max_cell_multiplier times the longest edge of"
            "indexing.unit_cell, or left to the DIALS nearest-neighbour estimate"
            "if no unit cell is given. The estimate is unreliable on sparse"
            "Laue-TOF spot lists."
    .expert_level = 2
  max_cell_multiplier = 1.3
    .type = float(value_min=1.0)
    .expert_level = 3
  d_min_start = 4.0
    .type = float
    .help = "Resolution limit for the first macrocycle of indexing refinement."
    .expert_level = 2
  min_spots = 10
    .type = int
    .help = "Only attempt to index an exposure with at least this many strong spots."
    .expert_level = 2
  min_indexed = 50
    .type = int
    .help = "Reject an indexing solution with fewer than this many indexed"
            "reflections, and move on to the next method in the ladder. This is"
            "a floor: min_indexed_fraction of the strong spots is used when that"
            "is the larger number."
    .expert_level = 2
  min_indexed_fraction = 0.25
    .type = float(value_min=0, value_max=1)
    .help = "Reject an indexing solution that indexes a smaller fraction of the"
            "strong spots than this. A fraction travels between datasets in a"
            "way an absolute count does not: a good MANDI solution indexes 65%"
            "of 227 spots, which an absolute cut of 250 would have thrown away,"
            "while a bad solution on 2300 NMX spots can index 300 of them."
    .expert_level = 2
  max_rmsd_px = 5.0
    .type = float
    .help = "Reject an indexing solution whose positional RMSD is worse than"
            "this, in pixels. Indexing plenty of reflections on a cell of the"
            "right size is not enough: fft3d has been seen to take 39% of a"
            "MANDI spot list at 6.1 px where another method fitted the same data"
            "at 3.2 px. Note this RMSD covers every indexed reflection, so it"
            "runs higher than the RMSD_X/RMSD_Y dials.index reports after"
            "outlier rejection. None turns the test off."
    .expert_level = 2
  target_rmsd_px = 2.0
    .type = float
    .help = "An accepted solution this good stops the ladder. Otherwise every"
            "method in the ladder is tried and the best by RMSD is kept, since"
            "the first method to give an acceptable answer is not necessarily"
            "the one that fits best. None makes the first acceptable solution"
            "win, as before."
    .expert_level = 2
  seed_from_first = True
    .type = bool
    .help = "Seed the indexing of each exposure with the crystal model of the"
            "first exposure that indexed. All exposures are measurements of the"
            "same crystal, so this both enforces a consistent setting and"
            "recovers orientations that will not index on their own."
    .expert_level = 1
  outlier {
    algorithm = auto *mcd tukey sauter_poon null
      .type = choice
      .help = "Outlier rejection algorithm used in indexing refinement. The"
              "residuals of sparse Laue-TOF data carry a heavy tail that a"
              "tight cut is needed to remove."
      .expert_level = 2
    iqr_multiplier = 0.8
      .type = float
      .expert_level = 3
    separate_panels = True
      .type = bool
      .expert_level = 3
  }
  phil = None
    .type = path
    .help = "Phil options file to use for indexing. Parameters defined in the"
            "xia2.laue_tof phil scope will take precedent over identical options"
            "defined in the phil file."
    .expert_level = 3
}
"""

integration_phil_str = """
integration {
  method = summation *profile_1d_ibix profile_3d_gutmann
    .type = choice
    .help = "Integration method. Summation intensities are written whichever"
            "method is chosen; a profile method adds intensity.prf alongside"
            "them. Profile fitting is not robust on all data - it can abort the"
            "whole run on a bad shoebox - so see fallback_method."
    .expert_level = 1
  fallback_method = *summation none
    .type = choice
    .help = "What to run when a profile method fails outright. The failure and"
            "the fallback are logged, and the fallback is then used for the"
            "remaining orientations rather than failing again on each."
    .expert_level = 2
  min_profile_fraction = 0.5
    .type = float(value_min=0, value_max=1)
    .help = "The fraction of reflections a profile method has to fit for its"
            "intensities to be usable. A profile run can finish successfully"
            "having fitted almost nothing - profile_3d_gutmann fitted 13 of"
            "1359 reflections on NMX test data - which is worse than useless,"
            "because the profile column then exists but covers a handful of"
            "reflections. Below this fraction the summation intensities are"
            "used instead and the remaining orientations skip the profile fit."
    .expert_level = 2
  integration_type = *observed calculated
    .type = choice
    .help = "observed integrates only reflections observed during spotfinding,"
            "calculated integrates all reflections out to calculated_d_min."
    .expert_level = 2
  calculated_d_min = None
    .type = float
    .help = "Resolution limit used when integration_type=calculated."
    .expert_level = 2
  background_model = constant2d constant3d linear2d *linear3d
    .type = choice
    .expert_level = 2
  mask = ellipse seed_skewness *both
    .type = choice
    .help = "Foreground/background mask method. both integrates twice, once"
            "with each mask, and keeps the better of the two (see"
            "mask_metric). The masks partition the shoebox differently, so"
            "which one is better is data dependent: seed_skewness typically"
            "takes in more of the peak tails and gives higher I/sigma, at"
            "several times the run time. The choice made on the first"
            "orientation is reused for the rest, so that every batch is"
            "integrated the same way."
    .expert_level = 2
  mask_metric = *i_over_sigma n_integrated
    .type = choice
    .help = "How to choose between the two masks when mask=both. Neither is a"
            "substitute for merging statistics, which are only available after"
            "scaling: treat the choice as a starting point and check it there."
    .expert_level = 3
  ellipse_mask_scale = 1.0
    .type = float(value_min=0.5)
    .expert_level = 3
  wavelength_range = None
    .type = floats(size=2)
    .help = "Reflections outside this wavelength range (Angstrom) are not integrated."
    .expert_level = 2
  bbox_tof_padding = 2
    .type = int
    .expert_level = 3
  bbox_xy_padding = 1
    .type = int
    .expert_level = 3
  lorentz = True
    .type = bool
    .help = "Apply the TOF Lorentz correction during integration. Note that the"
            "DIALS default is False. The correction must be applied exactly"
            "once: if it is to be applied by a downstream scaling program"
            "instead, set this to False."
    .expert_level = 1
  absorption {
    enabled = False
      .type = bool
      .help = "Apply the spherical absorption correction. Requires the target"
              "spectrum parameters below to be set."
      .expert_level = 2
    sample_number_density = None
      .type = float
      .help = "Sample number density (atoms/A^3) of the crystal."
      .expert_level = 2
    sample_radius = None
      .type = float
      .help = "Sample radius (mm) of the crystal."
      .expert_level = 2
    scattering_x_section = None
      .type = float
      .help = "Sample scattering cross section (barns) of the crystal."
      .expert_level = 2
    absorption_x_section = None
      .type = float
      .help = "Sample absorption cross section (barns) of the crystal."
      .expert_level = 2
  }
  phil = None
    .type = path
    .help = "Phil options file to use for integration. Parameters defined in the"
            "xia2.laue_tof phil scope will take precedent over identical options"
            "defined in the phil file."
    .expert_level = 3
}
"""

output_phil_str = """
output {
  intensity = *auto sum profile
    .type = choice
    .help = "Which intensities to export. auto uses the profile-fitted"
            "intensities when they are present and cover at least"
            "integration.min_profile_fraction of the reflections, and the"
            "summation intensities otherwise. dials.export refuses its own auto"
            "when both columns exist, so one is always chosen here explicitly."
    .expert_level = 1
  composition = CH
    .type = str
    .help = "Chemical composition of the asymmetric unit, written into the SHELX"
            ".ins file that accompanies the exported intensities."
    .expert_level = 2
  phil = None
    .type = path
    .help = "Phil options file to use for the export. Parameters defined in the"
            "xia2.laue_tof phil scope will take precedent over identical options"
            "defined in the phil file."
    .expert_level = 3
}
"""

workflow_phil_str = """
workflow {
  steps = *bin *find_spots *index *refine *integrate *combine *export
    .type = choice(multi=True)
    .help = "Option to turn off particular steps. Multiple choices should be of"
            "the format steps=find_spots+index"
    .expert_level = 3
}
"""

full_phil_str = (
    phil_str
    + binning_phil_str
    + spotfinding_phil_str
    + indexing_phil_str
    + integration_phil_str
    + output_phil_str
    + workflow_phil_str
)


def _resolved_file(value: str | None) -> pathlib.Path | None:
    """Resolve a phil path option, checking that the file exists."""
    if not value:
        return None
    path = pathlib.Path(value).resolve()
    if not path.is_file():
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), os.fspath(path)
        )
    return path


def _has_event_data(group: h5py.Group, max_depth: int = 3) -> bool:
    """
    Whether a NeXus group holds unbinned event data.

    Raw data from both NMX and MANDI stores its events in NXevent_data groups,
    which a TOF binned file does not have.
    """
    if max_depth < 0:
        return False
    nx_class = group.attrs.get("NX_class", "")
    if isinstance(nx_class, bytes):
        nx_class = nx_class.decode()
    if nx_class == "NXevent_data" or "event_time_offset" in group:
        return True
    return any(
        _has_event_data(item, max_depth - 1)
        for item in group.values()
        if isinstance(item, h5py.Group)
    )


def _instrument_name(entry: h5py.Group) -> str:
    """The instrument name an entry declares, if any."""
    try:
        name = entry["instrument"]["name"][()]
    except (KeyError, TypeError, ValueError):
        return ""
    if hasattr(name, "__len__") and not isinstance(name, bytes):
        if not len(name):
            return ""
        name = name[0]
    return name.decode() if isinstance(name, bytes) else str(name)


def _entry_definition(entry: h5py.Group) -> str:
    """The NeXus application definition of an entry, if it declares one."""
    for key in ("definition", "definitions"):
        if key in entry:
            value = entry[key][()]
            if isinstance(value, bytes):
                return value.decode()
            if hasattr(value, "__len__") and len(value):
                first = value[0]
                return first.decode() if isinstance(first, bytes) else str(first)
            return str(value)
    return ""


def classify_input(filename: pathlib.Path) -> str:
    """
    Identify what kind of data a NeXus file holds.

    Returns NXLAUETOF for data that has already been reduced to TOF bins and can
    be passed straight to dials.import, TOFRAW or NXSNSEVENT for raw data that
    must be reduced first, or UNKNOWN for anything that is not a NeXus file we
    recognise (e.g. ISIS SXD files, which dials.import reads directly).
    """
    try:
        with h5py.File(filename, "r") as handle:
            for name in handle:
                entry = handle[name]
                if not isinstance(entry, h5py.Group):
                    continue
                nx_class = entry.attrs.get("NX_class", "")
                if isinstance(nx_class, bytes):
                    nx_class = nx_class.decode()
                definition = _entry_definition(entry)
                if nx_class == "NXlauetof" or definition == "NXlauetof":
                    return NXLAUETOF
                if nx_class == "NXsnsevent" or definition == "NXsnsevent":
                    return NXSNSEVENT
                if definition == "TOFRAW":
                    return TOFRAW
                if _has_event_data(entry):
                    # Raw data does not always declare an application
                    # definition, so fall back on whether it holds events.
                    if _instrument_name(entry) == "MANDI":
                        return NXSNSEVENT
                    return TOFRAW
    except OSError:
        # Not an HDF5 file at all
        return UNKNOWN
    except KeyError:
        return UNKNOWN
    return UNKNOWN


@dataclass
class FileInput:
    images: list[str] = field(default_factory=list)
    templates: list[str] = field(default_factory=list)
    directories: list[str] = field(default_factory=list)
    vanadium_run: pathlib.Path | None = None
    empty_run: pathlib.Path | None = None
    mask: pathlib.Path | None = None
    import_phil: pathlib.Path | None = None

    def resolve_paths(self) -> None:
        for filetype in (self.images, self.templates):
            for i, obj in enumerate(filetype):
                filetype[i] = str(pathlib.Path(obj).resolve())

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> FileInput:
        file_input = cls()
        if params.image:
            file_input.images = list(params.image)
        elif params.template:
            file_input.templates = list(params.template)
        elif params.directory:
            file_input.directories = [
                str(pathlib.Path(i).resolve()) for i in params.directory
            ]
        else:
            raise ValueError(
                "No input data identified (use image=, template= or directory=)"
            )
        file_input.resolve_paths()
        file_input.vanadium_run = _resolved_file(params.vanadium_run)
        file_input.empty_run = _resolved_file(params.empty_run)
        file_input.mask = _resolved_file(params.mask)
        file_input.import_phil = _resolved_file(params.dials_import.phil)
        if file_input.empty_run and not file_input.vanadium_run:
            raise ValueError(
                "An empty_run was given without a vanadium_run. The empty run is"
                " only used to correct the incident spectrum, so both are needed."
            )
        return file_input


@dataclass
class BinningParams:
    enabled: bool | AutoType = Auto
    index_bins: int = 50
    integrate_bins: int = 200
    min_time_bin: float | None = None
    max_time_bin: float | None = None
    time_bin_unit: str = "us"
    detector_ids: list[int] | None = None
    tof_padding: float = 100.0
    essnmx_reduce: pathlib.Path | None = None
    extra_args: list[str] = field(default_factory=list)

    @property
    def single_reduction(self) -> bool:
        """Whether one reduction serves both indexing and integration."""
        return self.index_bins == self.integrate_bins

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> BinningParams:
        binning = params.binning
        return cls(
            binning.enabled,
            binning.index_bins,
            binning.integrate_bins,
            binning.min_time_bin,
            binning.max_time_bin,
            binning.time_bin_unit,
            list(binning.detector_ids) if binning.detector_ids else None,
            binning.tof_padding,
            _resolved_file(binning.essnmx_reduce),
            [arg for arg in binning.extra_args if arg],
        )


@dataclass
class SpotfindingParams:
    min_spot_size: int = 15
    max_spot_size: int | None = None
    threshold_algorithm: str = "dispersion_extended"
    max_strong: int = 35000
    d_min: float | None = None
    d_max: float | None = None
    nproc: int = 1
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> SpotfindingParams:
        return cls(
            params.spotfinding.min_spot_size,
            params.spotfinding.max_spot_size,
            params.spotfinding.threshold_algorithm,
            params.spotfinding.max_strong,
            params.d_min,
            params.d_max,
            params.multiprocessing.nproc,
            _resolved_file(params.spotfinding.phil),
        )


@dataclass
class IndexingParams:
    space_group: sgtbx.space_group | None = None
    unit_cell: uctbx.unit_cell | None = None
    methods: list[str] = field(default_factory=lambda: ["fft3d"])
    max_cell: float | None = None
    d_min_start: float | None = 4.0
    min_spots: int = 10
    min_indexed: int = 50
    min_indexed_fraction: float = 0.25
    max_rmsd_px: float | None = 5.0
    target_rmsd_px: float | None = 2.0
    seed_from_first: bool = True
    outlier_algorithm: str = "mcd"
    outlier_iqr_multiplier: float = 0.8
    outlier_separate_panels: bool = True
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> IndexingParams:
        indexing = params.indexing
        if indexing.unit_cell and params.space_group:
            try:
                _ = crystal.symmetry(
                    unit_cell=indexing.unit_cell,
                    space_group_info=params.space_group,
                    assert_is_compatible_unit_cell=True,
                )
            except AssertionError as e:
                raise ValueError(e)

        if indexing.method:
            methods = [indexing.method]
        else:
            methods = list(indexing.ladder)
        if not methods:
            raise ValueError(
                "No indexing method given (set indexing.method or indexing.ladder)"
            )
        for method in methods:
            if method not in INDEXING_METHODS:
                raise ValueError(
                    f"Unknown indexing method {method}."
                    f" Choose from {', '.join(INDEXING_METHODS)}"
                )

        max_cell = indexing.max_cell
        if max_cell is None and indexing.unit_cell:
            # The nearest-neighbour estimate is unreliable on the sparse spot
            # lists that Laue-TOF data gives, so derive a fixed value instead.
            max_cell = indexing.max_cell_multiplier * max(
                indexing.unit_cell.parameters()[:3]
            )

        return cls(
            params.space_group,
            indexing.unit_cell,
            methods,
            max_cell,
            indexing.d_min_start,
            indexing.min_spots,
            indexing.min_indexed,
            indexing.min_indexed_fraction,
            indexing.max_rmsd_px,
            indexing.target_rmsd_px,
            indexing.seed_from_first,
            indexing.outlier.algorithm,
            indexing.outlier.iqr_multiplier,
            indexing.outlier.separate_panels,
            _resolved_file(indexing.phil),
        )


@dataclass
class IntegrationParams:
    method: str = "profile_1d_ibix"
    fallback_method: str = "summation"
    min_profile_fraction: float = 0.5
    integration_type: str = "observed"
    calculated_d_min: float | None = None
    background_model: str = "linear3d"
    mask: str = "both"
    mask_metric: str = "i_over_sigma"
    ellipse_mask_scale: float = 1.0
    wavelength_range: tuple[float, float] | None = None
    bbox_tof_padding: int = 2
    bbox_xy_padding: int = 1
    lorentz: bool = True
    incident_run: pathlib.Path | None = None
    empty_run: pathlib.Path | None = None
    absorption: dict[str, float] | None = None
    nproc: int = 1
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(
        cls, params: iotbx.phil.scope_extract, file_input: FileInput
    ) -> IntegrationParams:
        integration = params.integration
        if integration.method not in INTEGRATION_METHODS:
            raise ValueError(
                f"Unknown integration method {integration.method}."
                f" Choose from {', '.join(INTEGRATION_METHODS)}"
            )
        if integration.integration_type == "calculated" and not (
            integration.calculated_d_min
        ):
            raise ValueError(
                "integration.calculated_d_min must be set when"
                " integration.integration_type=calculated"
            )

        absorption = None
        if integration.absorption.enabled:
            absorption = {
                "sample_number_density": integration.absorption.sample_number_density,
                "sample_radius": integration.absorption.sample_radius,
                "scattering_x_section": integration.absorption.scattering_x_section,
                "absorption_x_section": integration.absorption.absorption_x_section,
            }
            missing = [name for name, value in absorption.items() if value is None]
            if missing:
                raise ValueError(
                    "integration.absorption.enabled=True requires the target"
                    f" spectrum parameters, but {', '.join(missing)} not set"
                )

        wavelength_range = None
        if integration.wavelength_range:
            wavelength_range = tuple(integration.wavelength_range)

        return cls(
            integration.method,
            integration.fallback_method,
            integration.min_profile_fraction,
            integration.integration_type,
            integration.calculated_d_min,
            integration.background_model,
            integration.mask,
            integration.mask_metric,
            integration.ellipse_mask_scale,
            wavelength_range,
            integration.bbox_tof_padding,
            integration.bbox_xy_padding,
            integration.lorentz,
            file_input.vanadium_run,
            file_input.empty_run,
            absorption,
            params.multiprocessing.nproc,
            _resolved_file(integration.phil),
        )


@dataclass
class ExportParams:
    intensity: str = "auto"
    min_profile_fraction: float = 0.5
    composition: str = "CH"
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> ExportParams:
        return cls(
            params.output.intensity,
            params.integration.min_profile_fraction,
            params.output.composition,
            _resolved_file(params.output.phil),
        )


@dataclass
class AlgorithmParams:
    steps: list[str] = field(default_factory=list)
    nproc: int = 1

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> AlgorithmParams:
        return cls(list(params.workflow.steps), params.multiprocessing.nproc)


@dataclass
class LaueTOFSetup:
    """The fully resolved input of a run, as processed from the phil scope."""

    file_input: FileInput
    binning_params: BinningParams
    spotfinding_params: SpotfindingParams
    indexing_params: IndexingParams
    integration_params: IntegrationParams
    export_params: ExportParams
    options: AlgorithmParams
    # Input file class, keyed on the image path, for files given with image=
    input_classes: dict[str, str] = field(default_factory=dict)

    @property
    def files_to_bin(self) -> list[str]:
        """The input files that must be reduced before dials.import can read them."""
        if self.binning_params.enabled is False:
            return []
        if self.binning_params.enabled is True:
            return list(self.input_classes)
        return [
            filename
            for filename, file_class in self.input_classes.items()
            if file_class in (TOFRAW, NXSNSEVENT)
        ]


def _classify_file_input(file_input: FileInput) -> dict[str, str]:
    """Identify the data type of each input file given with image=."""
    input_classes = {}
    for image in file_input.images:
        path = pathlib.Path(image)
        if not path.is_file():
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), os.fspath(path)
            )
        input_classes[image] = classify_input(path)
    return input_classes


def _log_setup(setup: LaueTOFSetup) -> None:
    """Report what was resolved from the phil scope, before any work is done."""
    for filename, file_class in setup.input_classes.items():
        if file_class == NXLAUETOF:
            detail = "already TOF binned"
        elif file_class in (TOFRAW, NXSNSEVENT):
            detail = f"raw data ({file_class}), needs reducing"
        else:
            detail = "not a binnable NeXus file, will be imported as given"
        xia2_logger.info(f"Input {filename}: {detail}")

    to_bin = setup.files_to_bin
    binning = setup.binning_params
    if "bin" not in setup.options.steps:
        if to_bin:
            xia2_logger.warning(
                f"{len(to_bin)} input file(s) hold raw data, but the bin step is"
                " switched off. dials.import is unlikely to be able to read them."
            )
    elif to_bin:
        if binning.single_reduction:
            xia2_logger.info(
                f"Reducing {len(to_bin)} file(s) to {binning.index_bins} TOF bins"
            )
        else:
            xia2_logger.info(
                f"Reducing {len(to_bin)} file(s) to {binning.index_bins} TOF bins"
                f" for indexing and {binning.integrate_bins} for integration"
            )

    indexing = setup.indexing_params
    if len(indexing.methods) > 1:
        xia2_logger.info(
            "Indexing methods to try, in order: " + ", ".join(indexing.methods)
        )
    else:
        xia2_logger.info(f"Indexing method: {indexing.methods[0]}")
    if indexing.max_cell:
        xia2_logger.info(f"Using fixed max_cell of {indexing.max_cell:.1f} Angstrom")

    integration = setup.integration_params
    if integration.incident_run:
        xia2_logger.info(
            f"Normalising by the incident spectrum from {integration.incident_run}"
        )
        if not integration.empty_run:
            xia2_logger.warning(
                "No empty_run given, so the incident spectrum will not be"
                " background corrected."
            )
    else:
        xia2_logger.warning(
            "No vanadium_run given: the intensities will not be normalised by the"
            " incident spectrum."
        )
    xia2_logger.info(
        "TOF Lorentz correction will "
        + ("be applied during integration" if integration.lorentz else "not be applied")
    )
    if integration.mask == "both":
        xia2_logger.info(
            "Integrating the first orientation with each foreground mask and"
            f" keeping the better by {integration.mask_metric}"
        )
    xia2_logger.info(
        f"Integration method: {integration.method}"
        + (
            ""
            if integration.method == "summation"
            or integration.fallback_method == "none"
            else f", falling back to {integration.fallback_method} if it fails"
        )
    )


def setup_from_phil(params: iotbx.phil.scope_extract) -> LaueTOFSetup:
    """Turn the phil scope into the set of parameter objects the steps take."""
    if params.multiprocessing.nproc is Auto:
        params.multiprocessing.nproc = CPU_COUNT

    file_input = FileInput.from_phil(params)
    setup = LaueTOFSetup(
        file_input=file_input,
        binning_params=BinningParams.from_phil(params),
        spotfinding_params=SpotfindingParams.from_phil(params),
        indexing_params=IndexingParams.from_phil(params),
        integration_params=IntegrationParams.from_phil(params, file_input),
        export_params=ExportParams.from_phil(params),
        options=AlgorithmParams.from_phil(params),
        input_classes=_classify_file_input(file_input),
    )
    _log_setup(setup)
    return setup


@report_timing
def run_xia2_laue_tof(
    root_working_directory: pathlib.Path, params: iotbx.phil.scope_extract
) -> LaueTOFSetup:
    """
    Run the Laue-TOF processing pipeline.

    Each crystal orientation is processed independently, in its own
    subdirectory: reduce to TOF bins (if the input is not already binned),
    import, find spots, index, refine and integrate.
    """
    setup = setup_from_phil(params)

    try:
        from xia2.Modules.Laue_TOF.data_integration import run_data_integration
    except ImportError as e:
        raise NotImplementedError(
            "The processing steps are not implemented yet; this module currently"
            " only resolves and validates the input parameters."
        ) from e

    run_data_integration(root_working_directory, setup)
    return setup
