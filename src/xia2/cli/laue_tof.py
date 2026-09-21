"""
xia2.laue_tof: A processing pipeline for neutron Laue time-of-flight data,
using tools from the DIALS package.

The data are expected to be a single crystal measured as one or more stationary
exposures at different orientations, each exposure a stack of TOF bins. Each
orientation is processed independently: (bin) -> import -> find_spots -> index
-> refine -> tof_integrate.

With an already TOF-binned file (NXlauetof), run e.g.:
    xia2.laue_tof image=run_binned.h5 unit_cell=x space_group=y
Raw event data (NMX TOFRAW, MANDI NXsnsevent) is reduced to TOF bins first, with
essnmx, before importing:
    xia2.laue_tof image=run_events.h5 binning.nbins=50
Give one image= per crystal orientation, e.g.:
    xia2.laue_tof image=setting1.h5 image=setting2.h5 unit_cell=x space_group=y
An incident (vanadium) run and an empty-instrument run are optional but
recommended; they are passed on to dials.tof_integrate:
    xia2.laue_tof image=run.h5 vanadium_run=vanadium.nxs empty_run=empty.nxs

Refer to the individual DIALS program documentation for more details on
time-of-flight processing in DIALS.
"""

from __future__ import annotations

import logging
import pathlib
import sys
import traceback

import iotbx.phil
from dials.util.options import ArgumentParser

import xia2.Driver.timing
import xia2.Handlers.Streams
from xia2.Applications.xia2_main import write_citations
from xia2.Handlers.Citations import Citations
from xia2.Handlers.Files import cleanup
from xia2.Modules.Laue_TOF.laue_tof import full_phil_str, run_xia2_laue_tof

phil_scope = iotbx.phil.parse(full_phil_str)

xia2_logger = logging.getLogger(__name__)


def run(args=sys.argv[1:]):
    """
    Parse the command line input, setup logging and run the Laue-TOF pipeline.
    """
    Citations.cite("dials-integration")
    Citations.cite("xia2")
    parser = ArgumentParser(
        usage="xia2.laue_tof image=run_binned.h5 unit_cell=x space_group=y",
        read_experiments=False,
        read_reflections=False,
        phil=phil_scope,
        check_format=False,
        epilog=__doc__,
    )
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    xia2.Handlers.Streams.setup_logging(
        logfile="xia2.laue_tof.log", debugfile="xia2.laue_tof.debug.log"
    )
    # remove the xia2 handler from the dials logger.
    dials_logger = logging.getLogger("dials")
    dials_logger.handlers.clear()

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        xia2_logger.info("The following parameters have been modified:\n%s", diff_phil)

    cwd = pathlib.Path.cwd()
    try:
        with cleanup(cwd):
            run_xia2_laue_tof(cwd, params)
    except ValueError as e:
        xia2_logger.error(f"Error: {e}")
        sys.exit(0)
    except FileNotFoundError as e:
        xia2_logger.error(e)
        sys.exit(0)
    except Exception as e:
        with (cwd / "xia2-error.txt").open(mode="w") as fh:
            traceback.print_exc(file=fh)
        xia2_logger.error("Error: %s", str(e))
        xia2_logger.info(traceback.format_exc())
        xia2_logger.warning(
            "Please send the contents of xia2.laue_tof.log and xia2-error.txt to xia2.support@gmail.com"
        )
        sys.exit(1)

    write_citations(program="xia2.laue_tof")
