"""
Write an unmerged Laue MTZ from combined, integrated Laue-TOF data.

dials.export refuses format=mtz for time-of-flight data, and its MTZ writer
would in any case drop the per-reflection wavelength, so the file that pointless
and lawless need is written here instead.

Everything in the table below has been seen to fail silently, so the file is
checked after it is written rather than trusted:

    LAMBDA          per-observation wavelength; without it lawless stops on LAUE
    M/ISYM = 1      unreduced indices; anything else and pointless reads the
                    file as merged, collapses the observations and drops columns
    batch records   one per orientation, LDTYPE = 3, a cell, non-zero PHIRANGE,
                    BSCALE = 1 and the refined UMAT; a missing UMAT makes
                    SCALES SECONDARY and ABSORPTION refine flat, silently
    datasets        HKL_base as dataset 0, data in dataset 1, batch LBSETID
                    pointing at the data set

Both sets of intensities are written when the data carry both: summation as
I/SIGI and profile-fitted as IPR/SIGIPR, as MOSFLM labels them. Reflections
without a profile fit get the MTZ missing-number flag in IPR/SIGIPR rather than
a zero, so nothing downstream mistakes a failed fit for a measurement.
"""

from __future__ import annotations

import logging
import math
import pathlib

import gemmi
import numpy
from dials.array_family import flex
from dxtbx.model import Crystal, Experiment

xia2_logger = logging.getLogger(__name__)

# MTZ batch header word offsets, as used by CCP4 and gemmi. Taken from a working
# implementation rather than re-derived (see lauenorm2mtz.py, addbatches.py).
I_NWORDS, I_NINTGR, I_NREALS = 0, 1, 2
I_LBCELL = 4  # 4..9   cell refinement flags
I_MISFLG = 10
I_JUMPAX = 11
I_NCRYST = 12
I_LCRFLG = 13
I_LDTYPE = 14  # 1 = 2D oscillation, 2 = 3D area detector, 3 = Laue
I_JSCAXS = 15
I_NBSCAL = 16
I_NGONAX = 17
I_LBMFLG = 18
I_NDET = 19
I_LBSETID = 20

F_UMAT = 6  # 6..14
F_PHISTT = 36
F_PHIEND = 37
F_SCANAX = 38  # 38..40
F_BSCALE = 43
F_PHIRANGE = 47
F_E1 = 59  # 59..61  goniostat axis 1
F_SOURCE = 80  # 80..82
F_S0 = 83  # 83..85
F_ALAMBD = 86
F_DELAMB = 87
F_DX = 93  # 93..94
F_THETA = 95  # 95..96
F_DETLM = 97  # 97..104

# Batch numbers are 1-based in MTZ, while imageset_id/id are 0-based
FIRST_BATCH = 1


def _batch(
    number: int,
    title: str,
    cell: gemmi.UnitCell,
    crystal: Crystal,
    beam_direction: tuple[float, float, float],
    wavelength: float,
    half_band: float,
    dataset_id: int = 1,
) -> gemmi.Mtz.Batch:
    """
    One batch record for one stationary exposure.

    The crystal is stationary, so the exposure's orientation lives entirely in
    UMAT. phi is given a one degree pseudo-range per batch because lawless
    treats a batch with PHIRANGE = 0 as having no valid phi, and every
    observation then gets ROT = 0.
    """
    batch = gemmi.Mtz.Batch()
    batch.number = number
    batch.title = title[:70]
    batch.ints[I_NWORDS] = 185
    batch.ints[I_NINTGR] = 29
    batch.ints[I_NREALS] = 156
    for i in range(6):
        batch.ints[I_LBCELL + i] = 0
    batch.ints[I_MISFLG] = 0  # no separate missetting angles: U is complete
    batch.ints[I_JUMPAX] = 0
    batch.ints[I_NCRYST] = 1
    batch.ints[I_LCRFLG] = 0
    batch.ints[I_LDTYPE] = 3  # Laue
    batch.ints[I_JSCAXS] = 1  # scan axis = E1
    batch.ints[I_NBSCAL] = 0
    batch.ints[I_NGONAX] = 1
    batch.ints[I_LBMFLG] = 0
    batch.ints[I_NDET] = 1
    batch.ints[I_LBSETID] = dataset_id
    batch.axes = ["PHI"]
    batch.cell = cell

    # MTZBAT.umat is column-major, so the flat array holds U transposed
    umat = numpy.asarray(crystal.get_U(), dtype=float).reshape(3, 3)
    for i, value in enumerate(umat.T.reshape(9)):
        batch.floats[F_UMAT + i] = float(value)

    batch.floats[F_PHISTT] = float(number)
    batch.floats[F_PHIEND] = float(number) + 1.0
    batch.floats[F_PHIRANGE] = 1.0
    batch.floats[F_SCANAX + 0] = 1.0
    batch.floats[F_E1 + 2] = 1.0
    batch.floats[F_BSCALE] = 1.0
    # SOURCE and S0 are anti-parallel to the beam
    for i in range(3):
        batch.floats[F_SOURCE + i] = -float(beam_direction[i])
        batch.floats[F_S0 + i] = -float(beam_direction[i])
    batch.floats[F_ALAMBD] = wavelength
    batch.floats[F_DELAMB] = half_band
    batch.floats[F_DX + 0] = 100.0
    batch.floats[F_THETA + 0] = 0.0
    batch.floats[F_DETLM + 0] = 0.0
    batch.floats[F_DETLM + 1] = 100.0
    batch.floats[F_DETLM + 2] = 0.0
    batch.floats[F_DETLM + 3] = 100.0
    batch.wavelength = wavelength
    batch.dataset_id = dataset_id
    return batch


def _beam_direction(experiment: Experiment) -> tuple[float, float, float]:
    """The unit vector from the sample towards the detector."""
    direction = experiment.beam.get_unit_s0()
    length = math.sqrt(sum(value * value for value in direction))
    if not length:
        return (0.0, 0.0, 1.0)
    return tuple(value / length for value in direction)  # type: ignore[return-value]


def _has_profile_intensities(table: flex.reflection_table) -> bool:
    return (
        "intensity.prf.value" in table and "intensity.prf.variance" in table
    ) and table.get_flags(table.flags.integrated_prf).count(True) > 0


def write_unmerged_mtz(
    experiments_file: pathlib.Path,
    reflections_file: pathlib.Path,
    mtz_file: pathlib.Path,
    project: str = "xia2.laue_tof",
    crystal_name: str = "crystal",
    dataset_name: str = "laue_tof",
) -> pathlib.Path:
    """
    Write the combined integrated data as an unmerged Laue MTZ.

    Summation intensities go to I/SIGI and profile-fitted intensities, when the
    data carry them, to IPR/SIGIPR, so that one file serves either choice
    downstream and the two can be compared after scaling.
    """
    from dxtbx.serialize import load

    experiments = load.experiment_list(experiments_file, check_format=False)
    table = flex.reflection_table.from_file(reflections_file)

    integrated = table.select(
        table.get_flags(table.flags.integrated_sum)
        | table.get_flags(table.flags.integrated_prf)
    )
    if not integrated.size():
        raise ValueError(f"{reflections_file} holds no integrated reflections")
    if "wavelength_cal" not in integrated:
        raise ValueError(
            f"{reflections_file} has no wavelength_cal column, so it cannot be"
            " written as a Laue MTZ. Was it integrated by dials.tof_integrate?"
        )

    profile = _has_profile_intensities(integrated)
    cell = gemmi.UnitCell(*experiments[0].crystal.get_unit_cell().parameters())
    space_group = gemmi.SpaceGroup(
        experiments[0].crystal.get_space_group().type().lookup_symbol()
    )

    mtz = gemmi.Mtz(with_base=True)
    mtz.title = "Unmerged Laue time-of-flight data from xia2.laue_tof"
    mtz.cell = cell
    mtz.spacegroup = space_group
    mtz.add_dataset(dataset_name)
    mtz.datasets[1].project_name = project
    mtz.datasets[1].crystal_name = crystal_name
    mtz.datasets[1].cell = cell

    # M/ISYM and BATCH belong with H, K, L in the base dataset; pointless reads
    # the file as merged if M/ISYM is missing or not 1 throughout.
    for label, column_type in (("M/ISYM", "Y"), ("BATCH", "B")):
        mtz.add_column(label, column_type, dataset_id=0)
    columns: list[tuple[str, str]] = [("I", "J"), ("SIGI", "Q")]
    if profile:
        columns += [("IPR", "J"), ("SIGIPR", "Q")]
    columns += [("LAMBDA", "R"), ("ROT", "R")]
    detector_positions = "xyzobs.px.value" in integrated
    if detector_positions:
        columns += [("XDET", "R"), ("YDET", "R")]
    for label, column_type in columns:
        mtz.add_column(label, column_type, dataset_id=1)

    miller = integrated["miller_index"]
    wavelength = integrated["wavelength_cal"].as_numpy_array()
    batch_of_reflection = (
        integrated["imageset_id"] if "imageset_id" in integrated else integrated["id"]
    )
    batch_numbers = batch_of_reflection.as_numpy_array() + FIRST_BATCH

    rows = numpy.zeros((integrated.size(), len(columns) + 5), dtype=numpy.float32)
    rows[:, 0:3] = numpy.array([list(index) for index in miller], dtype=numpy.float32)
    rows[:, 3] = 1.0  # M/ISYM
    rows[:, 4] = batch_numbers
    column = 5

    intensity = integrated["intensity.sum.value"].as_numpy_array()
    sigma = numpy.sqrt(integrated["intensity.sum.variance"].as_numpy_array())
    rows[:, column] = intensity
    rows[:, column + 1] = sigma
    column += 2
    if profile:
        prf = integrated["intensity.prf.value"].as_numpy_array().copy()
        prf_sigma = numpy.sqrt(integrated["intensity.prf.variance"].as_numpy_array())
        fitted = integrated.get_flags(integrated.flags.integrated_prf).as_numpy_array()
        # A reflection the profile fit failed for has no profile intensity, so
        # flag it as missing instead of writing whatever is in the column.
        prf[~fitted] = numpy.nan
        prf_sigma[~fitted] = numpy.nan
        rows[:, column] = prf
        rows[:, column + 1] = prf_sigma
        column += 2
    rows[:, column] = wavelength
    column += 1
    rows[:, column] = 0.0  # ROT: the crystal is stationary within an exposure
    column += 1
    if detector_positions:
        x, y, _ = integrated["xyzobs.px.value"].parts()
        rows[:, column] = x.as_numpy_array()
        rows[:, column + 1] = y.as_numpy_array()

    mtz.set_data(rows)

    for index, experiment in enumerate(experiments):
        number = index + FIRST_BATCH
        selection = batch_numbers == number
        if not selection.any():
            xia2_logger.warning(f"Batch {number} has no reflections")
            band_centre, half_band = 1.0, 0.0
        else:
            batch_wavelengths = wavelength[selection]
            band_centre = float(batch_wavelengths.mean())
            half_band = 0.5 * float(batch_wavelengths.max() - batch_wavelengths.min())
        mtz.batches.append(
            _batch(
                number,
                f"orientation {number}",
                gemmi.UnitCell(*experiment.crystal.get_unit_cell().parameters()),
                experiment.crystal,
                _beam_direction(experiment),
                band_centre,
                half_band,
            )
        )
    mtz.datasets[1].wavelength = sum(b.wavelength for b in mtz.batches) / len(
        mtz.batches
    )

    mtz.history = [
        "Written by xia2.laue_tof from combined integrated Laue-TOF data",
        "I/SIGI are summation intensities"
        + (", IPR/SIGIPR profile fitted" if profile else ""),
    ]
    mtz.write_to_file(str(mtz_file))

    labels = "I/SIGI" + (" and IPR/SIGIPR" if profile else "")
    xia2_logger.info(
        f"Wrote {integrated.size()} observations in {len(mtz.batches)} batch(es)"
        f" to {mtz_file.name}, with {labels}"
    )
    validate_unmerged_mtz(mtz_file, expected_batches=len(experiments))
    return mtz_file


def check_wavelength_column(mtz_file: pathlib.Path) -> str:
    """
    The name of the per-observation wavelength column in an MTZ, or "".

    Used to check that a program in the scaling tail has not dropped the
    wavelength: pointless builds before 1.13.6 do, which makes LAUE fatal in
    lawless. lawless takes the first of these labels that it finds.
    """
    mtz = gemmi.read_mtz_file(str(mtz_file))
    labels = {column.label.upper(): column.label for column in mtz.columns}
    for candidate in ("LAMBDA", "LAM", "WAVELENGTH"):
        if candidate in labels:
            return labels[candidate]
    return ""


def validate_unmerged_mtz(
    mtz_file: pathlib.Path, expected_batches: int | None = None
) -> None:
    """
    Check the things that fail silently downstream.

    Raises ValueError naming the specific cause, because the symptoms - a
    collapsed merged file, one batch reported for everything, a flat secondary
    correction - do not point back here.
    """
    mtz = gemmi.read_mtz_file(str(mtz_file))
    problems = []

    if not check_wavelength_column(mtz_file):
        problems.append("no LAMBDA column, so lawless cannot run LAUE")

    labels = [column.label for column in mtz.columns]
    for required in ("I", "SIGI", "M/ISYM", "BATCH"):
        if required not in labels:
            problems.append(f"no {required} column")

    if "M/ISYM" in labels:
        misym = numpy.asarray(mtz.column_with_label("M/ISYM").array)
        if misym.size and not bool((misym == 1).all()):
            problems.append(
                "M/ISYM is not 1 throughout, so pointless will read the file as"
                " merged and collapse the observations"
            )

    if not mtz.batches:
        problems.append(
            "no batch records, so pointless will report one batch whatever the"
            " BATCH column says"
        )
    else:
        if expected_batches is not None and len(mtz.batches) != expected_batches:
            problems.append(
                f"{len(mtz.batches)} batch records for {expected_batches} orientations"
            )
        for batch in mtz.batches:
            if batch.ints[I_LDTYPE] != 3:
                problems.append(f"batch {batch.number} is not LDTYPE 3 (Laue)")
            if not batch.floats[F_PHIRANGE]:
                problems.append(f"batch {batch.number} has PHIRANGE 0")
            if batch.floats[F_BSCALE] != 1.0:
                problems.append(f"batch {batch.number} has BSCALE != 1")
            umat = numpy.array(
                [batch.floats[F_UMAT + i] for i in range(9)], dtype=float
            ).reshape(3, 3)
            if numpy.allclose(umat, numpy.eye(3)):
                problems.append(
                    f"batch {batch.number} has an identity UMAT, which makes"
                    " SCALES SECONDARY and ABSORPTION refine flat"
                )
            if numpy.linalg.det(umat) <= 0.9:
                problems.append(
                    f"batch {batch.number} has an invalid UMAT (det"
                    f" {numpy.linalg.det(umat):.3f})"
                )

    if problems:
        raise ValueError(
            f"{mtz_file} is not a usable unmerged Laue MTZ:\n  " + "\n  ".join(problems)
        )
