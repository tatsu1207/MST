"""
MST-Pipeline — PDF report generation for SourceTracker results.
"""
import json
from datetime import datetime
from pathlib import Path

import pandas as pd
from fpdf import FPDF


def _hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    h = hex_color.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _section_heading(pdf: FPDF, title: str):
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_text_color(44, 62, 80)
    pdf.cell(0, 7, title, new_x="LMARGIN", new_y="NEXT")
    pdf.set_draw_color(44, 62, 80)
    pdf.set_line_width(0.4)
    pdf.line(pdf.l_margin, pdf.get_y(), pdf.l_margin + 40, pdf.get_y())
    pdf.ln(3)
    pdf.set_text_color(0, 0, 0)


def _draw_stacked_bar_chart(pdf: FPDF, mp: pd.DataFrame, group_colors: dict):
    """Draw a stacked bar chart using fpdf2 drawing primitives."""
    avail_w = pdf.w - pdf.l_margin - pdf.r_margin
    n_samples = len(mp)
    if n_samples == 0:
        return

    # Adaptive sizing — shorter chart for fewer samples
    y_axis_w = 22
    chart_w = avail_w - y_axis_w
    chart_h = min(70, max(40, n_samples * 18))

    # Bar width: cap so single bars aren't comically wide
    max_bar_w = min(30, chart_w * 0.35 / max(n_samples, 1))
    bar_w = min(max((chart_w - 10) / (n_samples * 1.6), 16), max_bar_w)
    bar_gap = bar_w * 0.5
    total_bars_w = n_samples * bar_w + (n_samples - 1) * bar_gap
    x0 = pdf.l_margin + y_axis_w
    y0 = pdf.get_y()
    x_start = x0 + (chart_w - total_bars_w) / 2

    # Y-axis ticks + grid
    pdf.set_line_width(0.15)
    for tick in (0.0, 0.25, 0.5, 0.75, 1.0):
        y_tick = y0 + chart_h * (1 - tick)
        pdf.set_font("Helvetica", "", 6)
        pdf.set_text_color(90, 90, 90)
        pdf.set_xy(pdf.l_margin, y_tick - 1.5)
        pdf.cell(y_axis_w - 2, 3, f"{tick:.0%}", align="R")
        pdf.set_draw_color(210, 210, 210)
        pdf.set_dash_pattern(dash=0.8, gap=0.8)
        pdf.line(x0, y_tick, x0 + chart_w, y_tick)
    pdf.set_dash_pattern()

    # Baseline
    pdf.set_draw_color(160, 160, 160)
    pdf.set_line_width(0.25)
    pdf.line(x0, y0 + chart_h, x0 + chart_w, y0 + chart_h)

    # Bars
    for i, sample in enumerate(mp.index):
        x_bar = x_start + i * (bar_w + bar_gap)
        y_bottom = y0 + chart_h

        for col in mp.columns:
            val = mp.loc[sample, col]
            if val < 0.002:
                continue
            seg_h = val * chart_h
            y_top = y_bottom - seg_h

            r, g, b = _hex_to_rgb(group_colors.get(col, "#888888"))
            pdf.set_fill_color(r, g, b)
            pdf.set_draw_color(255, 255, 255)
            pdf.set_line_width(0.2)
            pdf.rect(x_bar, y_top, bar_w, seg_h, style="DF")

            # Percentage label inside segment
            if seg_h > 6 and val >= 0.05:
                pdf.set_text_color(255, 255, 255)
                fs = 6 if bar_w >= 25 else 5
                pdf.set_font("Helvetica", "B", fs)
                pdf.set_xy(x_bar, y_top + (seg_h - 3) / 2)
                pdf.cell(bar_w, 3, f"{val:.0%}", align="C")

            y_bottom = y_top

        # Sample label
        pdf.set_text_color(40, 40, 40)
        label_fs = min(7, max(5, int(bar_w / 3.5)))
        pdf.set_font("Helvetica", "", label_fs)
        label = str(sample)
        pdf.set_xy(x_bar - 4, y0 + chart_h + 1.5)
        pdf.cell(bar_w + 8, 4, label, align="C")

    # Legend — draw below chart
    legend_y = y0 + chart_h + 9
    pdf.set_y(legend_y)

    swatch_w, swatch_h = 6, 3
    cols_per_row = min(4, len(mp.columns))
    item_w = avail_w / cols_per_row

    for i, col in enumerate(mp.columns):
        row = i // cols_per_row
        col_idx = i % cols_per_row
        x = pdf.l_margin + col_idx * item_w
        y = legend_y + row * 6

        r, g, b = _hex_to_rgb(group_colors.get(col, "#888888"))
        pdf.set_fill_color(r, g, b)
        pdf.set_draw_color(r, g, b)
        pdf.rect(x, y + 0.8, swatch_w, swatch_h, style="F")
        pdf.set_text_color(40, 40, 40)
        pdf.set_font("Helvetica", "", 7)
        pdf.set_xy(x + swatch_w + 1.5, y)
        pdf.cell(item_w - swatch_w - 3, 5, col)

    n_legend_rows = (len(mp.columns) + cols_per_row - 1) // cols_per_row
    pdf.set_y(legend_y + n_legend_rows * 6 + 4)
    pdf.set_text_color(0, 0, 0)
    pdf.set_draw_color(0, 0, 0)
    pdf.set_line_width(0.2)


def _draw_table(pdf: FPDF, df: pd.DataFrame, title: str, fmt: str = ".4f"):
    """Draw a formatted data table with auto-sized columns."""
    _section_heading(pdf, title)

    avail_w = pdf.w - pdf.l_margin - pdf.r_margin
    n_data_cols = len(df.columns)

    # Measure widths: sample column gets enough room for longest name
    pdf.set_font("Helvetica", "", 6.5)
    max_name_w = max((pdf.get_string_width(str(s)) for s in df.index), default=15)
    sample_col_w = min(max(max_name_w + 4, 28), 65)
    data_col_w = (avail_w - sample_col_w) / n_data_cols

    row_h = 5

    # Header
    pdf.set_font("Helvetica", "B", 6.5)
    pdf.set_fill_color(44, 62, 80)
    pdf.set_text_color(255, 255, 255)
    pdf.cell(sample_col_w, row_h, "Sample", border=1, fill=True, align="C")
    for c in df.columns:
        pdf.cell(data_col_w, row_h, str(c), border=1, fill=True, align="C")
    pdf.ln()

    # Data rows
    pdf.set_font("Helvetica", "", 6.5)
    for i, sample in enumerate(df.index):
        if i % 2 == 0:
            pdf.set_fill_color(243, 245, 249)
        else:
            pdf.set_fill_color(255, 255, 255)
        pdf.set_text_color(30, 30, 30)
        pdf.cell(sample_col_w, row_h, str(sample), border=1, fill=True, align="L")
        pdf.set_text_color(60, 60, 60)
        for c in df.columns:
            val = df.loc[sample, c]
            pdf.cell(data_col_w, row_h, f"{val:{fmt}}", border=1, fill=True, align="C")
        pdf.ln()

    pdf.set_text_color(0, 0, 0)
    pdf.ln(3)


def generate_st_report(
    run_id: int,
    run_dir: Path,
    run_name: str | None = None,
    dataset_name: str | None = None,
    variable_region: str | None = None,
    run_params: dict | None = None,
) -> Path:
    """Generate a PDF report for a SourceTracker run. Returns path to PDF."""
    from app.pipeline.sourcetracker import GROUP_COLORS

    mp = pd.read_csv(run_dir / "proportions.csv", index_col=0)
    stds = None
    stds_path = run_dir / "stds.csv"
    if stds_path.exists():
        stds = pd.read_csv(stds_path, index_col=0)

    status_info = {}
    status_file = run_dir / "status.json"
    if status_file.exists():
        try:
            status_info = json.loads(status_file.read_text())
        except Exception:
            pass

    # -- Build PDF -------------------------------------------------------------
    pdf = FPDF(orientation="P", unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=12)
    pdf.add_page()

    # Title banner
    pdf.set_fill_color(44, 62, 80)
    pdf.rect(0, 0, pdf.w, 24, style="F")
    pdf.set_text_color(255, 255, 255)
    pdf.set_font("Helvetica", "B", 18)
    pdf.set_y(5)
    pdf.cell(0, 8, "Source Tracking Report", align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 8)
    pdf.set_text_color(180, 200, 220)
    pdf.cell(0, 5, f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}", align="C", new_x="LMARGIN", new_y="NEXT")

    pdf.set_y(30)
    pdf.set_text_color(0, 0, 0)

    # -- Run Information -------------------------------------------------------
    _section_heading(pdf, "Run Information")

    info_rows = []
    if run_name:
        info_rows.append(("Run Name", run_name))
    if dataset_name:
        info_rows.append(("Dataset", dataset_name))
    if variable_region:
        info_rows.append(("Variable Region", variable_region))
    info_rows.append(("Samples", str(len(mp))))
    info_rows.append(("Source Groups", str(len(mp.columns))))

    if run_params:
        for key in ("feature_mode", "src_depth", "snk_depth", "restarts", "burnin", "draws"):
            val = run_params.get(key)
            if val is not None:
                label = key.replace("_", " ").title()
                info_rows.append((label, str(val)))

    alignment_matched = status_info.get("alignment_matched")
    alignment_total = status_info.get("alignment_total")
    if alignment_matched is not None and alignment_total is not None:
        pct = (alignment_matched / alignment_total * 100) if alignment_total else 0
        info_rows.append(("Alignment", f"{alignment_matched}/{alignment_total} features ({pct:.1f}%)"))

    col_w = (pdf.w - pdf.l_margin - pdf.r_margin) / 2
    for i, (label, value) in enumerate(info_rows):
        if i % 2 == 0:
            x_start = pdf.l_margin
            y_row = pdf.get_y()
        else:
            x_start = pdf.l_margin + col_w

        pdf.set_xy(x_start, y_row)
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_text_color(80, 80, 80)
        pdf.cell(28, 4.5, f"{label}:", align="R")
        pdf.set_font("Helvetica", "", 8)
        pdf.set_text_color(30, 30, 30)
        pdf.cell(col_w - 30, 4.5, f"  {value}")

        if i % 2 == 1 or i == len(info_rows) - 1:
            pdf.set_y(y_row + 5)

    pdf.ln(2)

    # Low alignment warnings
    low_samples = status_info.get("low_alignment_samples", [])
    sample_pcts = status_info.get("sample_alignment_pct", {})
    if low_samples:
        pdf.set_fill_color(255, 248, 235)
        pdf.set_draw_color(220, 160, 50)
        pdf.set_line_width(0.3)
        box_y = pdf.get_y()
        box_h = 6 + len(low_samples) * 4
        pdf.rect(pdf.l_margin, box_y, pdf.w - pdf.l_margin - pdf.r_margin, box_h, style="DF")
        pdf.set_xy(pdf.l_margin + 3, box_y + 1.5)
        pdf.set_font("Helvetica", "B", 7)
        pdf.set_text_color(160, 80, 0)
        pdf.cell(0, 3.5, f"Warning: {len(low_samples)} sample(s) excluded (low alignment <10%)")
        for s in low_samples:
            pdf.set_xy(pdf.l_margin + 6, pdf.get_y() + 4)
            pdf.set_font("Helvetica", "", 7)
            pct_val = sample_pcts.get(s, 0)
            pdf.cell(0, 3.5, f"- {s}: {pct_val}% reads retained")
        pdf.set_y(box_y + box_h + 3)
        pdf.set_text_color(0, 0, 0)

    # -- Chart -----------------------------------------------------------------
    _section_heading(pdf, "Source Contributions")
    _draw_stacked_bar_chart(pdf, mp, GROUP_COLORS)

    # -- Tables ----------------------------------------------------------------
    _draw_table(pdf, mp.round(4), "Mixing Proportions")

    if stds is not None and not stds.empty:
        _draw_table(pdf, stds.round(4), "Standard Deviations")

    # -- Footer ----------------------------------------------------------------
    pdf.ln(4)
    pdf.set_font("Helvetica", "I", 6)
    pdf.set_text_color(160, 160, 160)
    pdf.cell(0, 3, "MST-Pipeline  |  SourceTracker2 Gibbs Sampling", align="C")

    # Save
    pdf_path = run_dir / "source_tracking_report.pdf"
    pdf.output(str(pdf_path))

    return pdf_path
