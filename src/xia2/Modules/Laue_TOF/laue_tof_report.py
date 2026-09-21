"""
The run report, in the format the rest of xia2 uses.

`xia2.Modules.Report` already produces the page a crystallographer expects -
CC1/2 and I/sigma against resolution, completeness, multiplicity, Rmerge against
batch, second moments, the L test, the Wilson plot, the multiplicity images -
from an unmerged MTZ with BATCH and I/SIGI. That is exactly what lawless writes,
so the report here is that same page, rendered from the same template as
xia2.report and multiplex.

What the standard page cannot know about is the Laue part, so one panel is added
to it: the fitted wavelength normalisation curve, the wavelength band each
orientation contributed, and the parameters lawless fitted. Those come from the
LAMBDANORM file and the XML that lawless writes beside its log.
"""

from __future__ import annotations

import json
import logging
import os
import pathlib
from collections import OrderedDict

import gemmi
import numpy
from jinja2 import ChoiceLoader, Environment, PackageLoader

from xia2.Handlers.Files import FileHandler
from xia2.Handlers.Streams import banner
from xia2.Modules.Analysis import phil_scope
from xia2.Modules.Laue_TOF.laue_tof_mtz import F_ALAMBD, F_DELAMB
from xia2.Modules.Laue_TOF.laue_tof_scale import lawless_xml
from xia2.Modules.Report import Report
from xia2.XIA2Version import Version

xia2_logger = logging.getLogger(__name__)

PREFIX = "xia2.laue_tof"

# The fitted parameters worth putting in the report, per normalisation method
NORMALISATION_PARAMETERS = (
    ("ReferenceWavelength", "Reference wavelength (A)"),
    ("LambdaMin", "Lowest wavelength fitted (A)"),
    ("LambdaMax", "Highest wavelength fitted (A)"),
    ("Kernel", "Kernel"),
    ("LengthScale", "Length scale (A)"),
    ("SigmaF", "Sigma f"),
    ("NoiseInflation", "Noise inflation"),
    ("TrainingBins", "Training bins"),
    ("Degree", "Chebyshev degree"),
)


def parse_lambdanorm(path: pathlib.Path) -> dict[str, list[list[float]]]:
    """
    The curve and the binned observations from a lawless LAMBDANORM file.

    The file is a gnuplot script with the numbers in named inline data blocks,
    $LAMBDANORM (wavelength, w, relative uncertainty, w low, w high) and
    $LAMBDABINS (wavelength, w, sigma), so the blocks are read and the gnuplot
    around them ignored.
    """
    blocks: dict[str, list[list[float]]] = {}
    current: str | None = None
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("$") and "<<" in stripped:
            current = stripped[1:].split()[0]
            blocks[current] = []
            continue
        if current and stripped in ("EOD", "EOF"):
            current = None
            continue
        if not current or not stripped or stripped.startswith("#"):
            continue
        try:
            blocks[current].append([float(value) for value in stripped.split()])
        except ValueError:
            continue
    return blocks


def normalisation_plot(lambdanorm: pathlib.Path) -> dict | None:
    """The wavelength normalisation curve, with its uncertainty band."""
    blocks = parse_lambdanorm(lambdanorm)
    curve = blocks.get("LAMBDANORM")
    if not curve:
        return None
    rows = numpy.asarray(curve, dtype=float)
    wavelength, scale = rows[:, 0].tolist(), rows[:, 1].tolist()
    data = []
    if rows.shape[1] >= 5:
        low, high = rows[:, 3].tolist(), rows[:, 4].tolist()
        # One filled trace for the band: up the upper edge and back along the lower
        data.append(
            {
                "x": wavelength + wavelength[::-1],
                "y": high + low[::-1],
                "type": "scatter",
                "fill": "toself",
                "fillcolor": "rgba(31,119,180,0.2)",
                "line": {"color": "rgba(0,0,0,0)"},
                "name": "1 sigma",
                "hoverinfo": "skip",
            }
        )
    data.append(
        {
            "x": wavelength,
            "y": scale,
            "type": "scatter",
            "mode": "lines",
            "name": "w(lambda)",
        }
    )
    bins = blocks.get("LAMBDABINS")
    if bins:
        binned = numpy.asarray(bins, dtype=float)
        trace = {
            "x": binned[:, 0].tolist(),
            "y": binned[:, 1].tolist(),
            "type": "scatter",
            "mode": "markers",
            "name": "binned observations",
        }
        if binned.shape[1] >= 3:
            trace["error_y"] = {"type": "data", "array": binned[:, 2].tolist()}
        data.append(trace)
    return {
        "data": data,
        "layout": {
            "title": "Wavelength normalisation",
            "xaxis": {"title": "Wavelength (A)"},
            "yaxis": {"title": "w(lambda)", "rangemode": "tozero"},
        },
    }


def wavelength_band_plot(mtz_file: pathlib.Path) -> dict | None:
    """The wavelength band each orientation contributed, from the batch headers."""
    mtz = gemmi.read_mtz_file(os.fspath(mtz_file))
    if not mtz.batches:
        return None
    numbers = [batch.number for batch in mtz.batches]
    centres = [batch.floats[F_ALAMBD] for batch in mtz.batches]
    widths = [batch.floats[F_DELAMB] for batch in mtz.batches]
    return {
        "data": [
            {
                "x": numbers,
                "y": centres,
                "error_y": {"type": "data", "array": widths},
                "type": "scatter",
                "mode": "markers",
                "name": "band centre",
            }
        ],
        "layout": {
            "title": "Wavelength band per orientation",
            "xaxis": {"title": "Batch"},
            "yaxis": {"title": "Wavelength (A)", "rangemode": "tozero"},
        },
    }


def intensity_comparison_plot(unmerged_mtz: pathlib.Path) -> dict | None:
    """
    Summation against profile-fitted intensities, when the data has both.

    Taken from the unscaled unmerged file, which is the only one that carries
    both: what lawless scales is one of them.
    """
    mtz = gemmi.read_mtz_file(os.fspath(unmerged_mtz))
    labels = [column.label for column in mtz.columns]
    if "I" not in labels or "IPR" not in labels:
        return None
    summation = numpy.asarray(mtz.column_with_label("I").array)
    profile = numpy.asarray(mtz.column_with_label("IPR").array)
    both = numpy.isfinite(summation) & numpy.isfinite(profile)
    if not both.any():
        return None
    return {
        "data": [
            {
                "x": summation[both].tolist(),
                "y": profile[both].tolist(),
                "type": "scatter",
                "mode": "markers",
                "marker": {"size": 3, "opacity": 0.4},
                "name": f"{int(both.sum())} reflections",
            }
        ],
        "layout": {
            "title": "Profile-fitted against summation intensities (unscaled)",
            "xaxis": {"title": "I (summation)"},
            "yaxis": {"title": "I (profile)"},
        },
    }


def normalisation_table(xmlout: pathlib.Path) -> list[list[str]] | None:
    """What lawless fitted, as a table for the report."""
    if not xmlout.is_file():
        return None
    root = lawless_xml(xmlout)
    if root is None:
        return None
    for tag in ("WavelengthNormalisationGPR", "WavelengthNormalisationChebyshev"):
        block = root.find(tag)
        if block is None:
            continue
        rows = [
            ["Wavelength normalisation", tag.replace("WavelengthNormalisation", "")]
        ]
        for name, label in NORMALISATION_PARAMETERS:
            element = block.find(name)
            if element is not None and element.text and element.text.strip():
                rows.append([label, element.text.strip()])
        return rows
    return None


def _laue_graphs(
    scaled_unmerged_mtz: pathlib.Path,
    unmerged_mtz: pathlib.Path | None,
    lawless_directory: pathlib.Path | None,
) -> OrderedDict:
    """The Laue-specific plots, in the order they should appear."""
    graphs: OrderedDict = OrderedDict()
    if lawless_directory:
        lambdanorm = lawless_directory / "LAMBDANORM"
        if lambdanorm.is_file():
            plot = normalisation_plot(lambdanorm)
            if plot:
                graphs["wavelength_normalisation"] = plot
    band = wavelength_band_plot(scaled_unmerged_mtz)
    if band:
        graphs["wavelength_band"] = band
    if unmerged_mtz and unmerged_mtz.is_file():
        comparison = intensity_comparison_plot(unmerged_mtz)
        if comparison:
            graphs["intensity_comparison"] = comparison
    return graphs


def generate_report(
    working_directory: pathlib.Path,
    scaled_unmerged_mtz: pathlib.Path,
    unmerged_mtz: pathlib.Path | None = None,
    lawless_directory: pathlib.Path | None = None,
    prefix: str = PREFIX,
) -> pathlib.Path:
    """
    Write the html and json report for a finished run.

    The standard xia2 analysis comes from Report, and one extra panel carries the
    Laue-specific plots and the normalisation parameters.
    """
    xia2_logger.notice(banner("Writing the report"))  # type: ignore
    params = phil_scope.extract()
    params.dose.batch = []
    params.batch = []

    report = Report.from_unmerged_mtz(
        scaled_unmerged_mtz, params, report_dir=os.fspath(working_directory)
    )
    overall_stats_table, merging_stats_table, stats_plots = (
        report.resolution_plots_and_stats()
    )

    json_data: dict = {}
    json_data.update(stats_plots)
    json_data.update(report.batch_dependent_plots())
    json_data.update(report.intensity_stats_plots(run_xtriage=False))

    resolution_graphs = OrderedDict(
        (key, json_data[key])
        for key in (
            "cc_one_half",
            "i_over_sig_i",
            "second_moments",
            "wilson_intensity_plot",
            "completeness",
            "multiplicity_vs_resolution",
        )
        if key in json_data
    )
    batch_graphs = OrderedDict(
        (key, json_data[key])
        for key in ("scale_rmerge_vs_batch", "i_over_sig_i_vs_batch")
        if key in json_data
    )
    misc_graphs = OrderedDict(
        (key, json_data[key])
        for key in ("cumulative_intensity_distribution", "l_test", "multiplicities")
        if key in json_data
    )
    laue_graphs = _laue_graphs(scaled_unmerged_mtz, unmerged_mtz, lawless_directory)
    json_data.update(laue_graphs)

    laue_tables = []
    if lawless_directory:
        table = normalisation_table(lawless_directory / "lawless.xml")
        if table:
            laue_tables.append(table)

    environment = Environment(
        loader=ChoiceLoader(
            [PackageLoader("xia2", "templates"), PackageLoader("dials", "templates")]
        )
    )
    template = environment.get_template("laue_tof_report.html")
    html = template.render(
        page_title="xia2.laue_tof report",
        filename=os.fspath(scaled_unmerged_mtz.resolve()),
        space_group=report.intensities.space_group_info().symbol_and_number(),
        unit_cell=str(report.intensities.unit_cell()),
        mtz_history=[line.strip() for line in report.mtz_object.history()],
        xtriage_success=None,
        xtriage_warnings=None,
        xtriage_danger=None,
        overall_stats_table=overall_stats_table,
        merging_stats_table=merging_stats_table,
        cc_half_significance_level=params.cc_half_significance_level,
        resolution_graphs=resolution_graphs,
        batch_graphs=batch_graphs,
        misc_graphs=misc_graphs,
        laue_graphs=laue_graphs,
        laue_tables=laue_tables,
        styles={},
        xia2_version=Version,
        log_text="",
    )

    html_file = working_directory / f"{prefix}-report.html"
    json_file = working_directory / f"{prefix}-report.json"
    html_file.write_bytes(html.encode("utf-8", "xmlcharrefreplace"))
    json_file.write_text(json.dumps(json_data, indent=None))
    FileHandler.record_html_file(f"{prefix} report", str(html_file))
    FileHandler.record_data_file(str(json_file))
    xia2_logger.info(
        f"Wrote {html_file.name}, with {len(laue_graphs)} Laue plot(s) beside the"
        " standard analysis"
    )
    return html_file
