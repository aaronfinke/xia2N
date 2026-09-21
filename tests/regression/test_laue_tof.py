from __future__ import annotations

import os
import pathlib
import shutil
import subprocess

import h5py
import pytest
from dxtbx.serialize import load

from xia2.Modules.Laue_TOF.laue_tof import NXSNSEVENT, classify_input

# The MANDI test data is a 3 GB raw event file, too big for dials_data, so the
# MANDI tests run only where a copy exists. Default location is work/mandi in
# the repository (gitignored); XIA2_LAUE_TOF_MANDI overrides it.
MANDI_DEFAULT = (
    pathlib.Path(__file__).parents[2] / "work" / "mandi" / "MANDI_13378.nxs.h5"
)


@pytest.fixture
def mandi_raw() -> pathlib.Path:
    """The raw MANDI event file, or skip if there is no copy of it here."""
    path = pathlib.Path(os.environ.get("XIA2_LAUE_TOF_MANDI", MANDI_DEFAULT))
    if not path.is_file():
        pytest.skip(
            f"MANDI test data not found at {path}. Copy MANDI_13378.nxs.h5 there,"
            " or point XIA2_LAUE_TOF_MANDI at it, to run the MANDI tests."
        )
    return path


def mandi_histogram(image: pathlib.Path) -> tuple[int, list[str]]:
    """The TOF bin count and the panels holding histogram data in an event file."""
    with h5py.File(image, "r") as fh:
        entry = fh["entry"]
        if "time_of_flight" not in entry:
            return 0, []
        panels = [
            name
            for name in entry
            if name.endswith("_events") and "spectra" in entry[name]
        ]
        return len(entry["time_of_flight"]) - 1, panels


def test_classify_mandi_raw(mandi_raw):
    """Raw MANDI event data is recognised as needing binning."""
    assert classify_input(mandi_raw) == NXSNSEVENT


def test_mandi_histogram_and_import(regression_test, mandi_raw, tmp_path):
    """
    Bin a raw MANDI run and import the result.

    MANDI data is binned in place: dxtbx histograms the events into the event
    file itself (FormatMANDI.add_histogram_data_to_nxs_file), which is what
    FormatMANDI reads, and the event data is kept. The first run over a file
    writes the histogram, which takes a while; later runs reuse it.

    Note this test therefore modifies the data file it is given.
    """
    result = subprocess.run(
        [
            shutil.which("xia2.laue_tof"),
            f"image={os.fspath(mandi_raw)}",
            "binning.index_bins=50",
            "binning.integrate_bins=50",
            "workflow.steps=bin",
        ],
        cwd=tmp_path,
        capture_output=True,
        encoding="utf-8",
    )

    n_bins, panels = mandi_histogram(mandi_raw)
    assert n_bins, result.stdout
    assert len(panels) > 1

    with h5py.File(mandi_raw, "r") as fh:
        # remove_event_data=False: the events must still be there afterwards
        assert "event_time_offset" in fh["entry"][panels[0]]

    imported = tmp_path / "orientation_1" / "imported.expt"
    assert imported.is_file(), result.stdout
    expts = load.experiment_list(imported, check_format=False)
    assert expts.all_tof()
    assert len(expts[0].detector) > 1
    assert expts[0].scan.get_num_images() == n_bins
