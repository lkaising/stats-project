"""Markdown formatter for the analysis report."""

import math


def _fmt_num(v, digits=6):
    if isinstance(v, bool):
        return str(v)
    if v is None:
        return "—"
    if isinstance(v, float):
        if math.isnan(v):
            return "NaN"
        if abs(v) >= 1e4 or (abs(v) > 0 and abs(v) < 1e-3):
            return f"{v:.3e}"
        return f"{v:.{digits}f}"
    return str(v)


def _fmt_p(p):
    if p is None:
        return "—"
    if isinstance(p, float) and math.isnan(p):
        return "NaN"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _table(headers, rows):
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def _assumption_section(assumptions):
    lines = ["## Assumption checks", ""]
    for metric in ("weber", "michelson"):
        a = assumptions[metric]
        lines.append(f"### {metric.capitalize()}")
        lines.append("")
        rows = [
            [
                s["condition"],
                _fmt_num(s["W"]),
                _fmt_p(s["p"]),
                "yes" if s["violated"] else "no",
            ]
            for s in a["shapiro_per_condition"]
        ]
        lines.append("**Shapiro-Wilk (per condition, on site-aggregated means, n=4)**")
        lines.append("")
        lines.append(_table(["condition", "W", "p", "violated"], rows))
        lines.append("")
        lines.append(
            f"- Normality violated for any condition: **{a['normality_violated']}**"
        )
        lines.append(
            f"- Friedman fallback triggered: **{a['friedman_fallback_triggered']}**"
        )
        lines.append("")
        m = a["mauchly"]
        lines.append("**Mauchly's test of sphericity (on the 4 × 5 site-aggregated matrix)**")
        lines.append("")
        if m["status"] == "noncomputable":
            note = m.get("note") or m.get("error") or ""
            lines.append(
                f"- Status: **noncomputable** — {note}"
            )
            lines.append(
                "- Per the locked contingency, Greenhouse-Geisser correction is "
                "applied conservatively if the RM ANOVA path is taken."
            )
        else:
            lines.append(
                f"- W={_fmt_num(m['W'])}, χ²={_fmt_num(m['chi2'])}, "
                f"df={_fmt_num(m['dof'])}, p={_fmt_p(m['p'])}"
            )
            lines.append(f"- Violated: **{m['violated']}** (reason: {m['reason']})")
        lines.append(
            f"- Sphericity correction applied: **{a['sphericity_correction_applied']}** "
            f"(reason: {a['sphericity_reason']})"
        )
        lines.append("")
    return "\n".join(lines)


def _means_table(means):
    rows = [
        [
            m["condition"],
            _fmt_num(m["mean"]),
            _fmt_num(m["sd"]),
            str(m["n_sites"]),
            _fmt_num(m["ci_low"]),
            _fmt_num(m["ci_high"]),
        ]
        for m in means
    ]
    return _table(
        ["condition", "site-level mean", "sd", "n_sites", "95% CI low", "95% CI high"],
        rows,
    )


def _omnibus_block(res):
    path = res["path"]
    o = res["omnibus"]
    effect = res["effect_size"]
    lines = []
    if path == "rm_anova":
        lines.append(
            f"**Path**: one-way repeated-measures ANOVA (site = repeated-measures unit, condition = within-unit factor)."
        )
        lines.append("")
        lines.append(
            f"- F({_fmt_num(o['df1'])}, {_fmt_num(o['df2'])}) = {_fmt_num(o['F'])}"
        )
        lines.append(f"- Uncorrected p = {_fmt_p(o['p_unc'])}")
        lines.append(
            f"- Greenhouse-Geisser corrected p = {_fmt_p(o['p_gg'])} (ε_GG = {_fmt_num(o['eps_gg'])})"
        )
        lines.append(
            f"- Reported p (GG applied = {o['sphericity_correction_applied']}): "
            f"**{_fmt_p(o['p_reported'])}**"
        )
        lines.append(
            f"- Effect size (partial η²): {_fmt_num(effect['value'])}"
        )
    else:
        lines.append(
            f"**Path**: Friedman test (nonparametric fallback — normality violated)."
        )
        lines.append("")
        lines.append(
            f"- Q = {_fmt_num(o['Q'])}, df = {_fmt_num(o['dof'])}, p = {_fmt_p(o['p'])}"
        )
        lines.append(
            f"- Effect size (Kendall's W): {_fmt_num(effect['value'])}"
        )
    lines.append(f"- Significant at α = {o['alpha']}: **{o['significant']}**")
    return "\n".join(lines)


def _posthoc_block(res):
    posthoc = res["posthoc"]
    if not posthoc:
        return "- Post-hoc comparisons were not run (omnibus non-significant)."
    path = res["path"]
    if path == "rm_anova":
        headers = [
            "pair",
            "t",
            "df",
            "mean_diff",
            "95% CI low",
            "95% CI high",
            "p_raw",
            "p_holm",
            "reject (α=0.05)",
        ]
        rows = [
            [
                f"{r['a']} vs {r['b']}",
                _fmt_num(r["t"]),
                str(r["df"]),
                _fmt_num(r["mean_diff"]),
                _fmt_num(r["ci_low"]),
                _fmt_num(r["ci_high"]),
                _fmt_p(r["p_raw"]),
                _fmt_p(r["p_holm"]),
                str(r["reject_holm_at_0.05"]),
            ]
            for r in posthoc
        ]
    else:
        headers = [
            "pair",
            "Wilcoxon W",
            "median_diff",
            "p_raw",
            "p_holm",
            "reject (α=0.05)",
        ]
        rows = [
            [
                f"{r['a']} vs {r['b']}",
                _fmt_num(r["wilcoxon_W"]),
                _fmt_num(r["median_diff"]),
                _fmt_p(r["p_raw"]),
                _fmt_p(r["p_holm"]),
                str(r["reject_holm_at_0.05"]),
            ]
            for r in posthoc
        ]
    return _table(headers, rows)


def _a1_section(title, res):
    lines = [f"## {title}", ""]
    lines.append(
        f"Outcome: {res['metric']}. Unit of analysis: site-level mean (n={res['n_sites']} sites)."
    )
    lines.append(
        f"Fallback used (RM ANOVA → Friedman): **{res['assumption_path']['fallback_used']}**."
    )
    lines.append("")
    lines.append("### Condition means (site-aggregated, 95% CI)")
    lines.append("")
    lines.append(_means_table(res["condition_means"]))
    lines.append("")
    lines.append("### Omnibus")
    lines.append("")
    lines.append(_omnibus_block(res))
    lines.append("")
    lines.append("### Post-hoc pairwise comparisons")
    lines.append("")
    if res["posthoc"]:
        path_note = (
            "Wilcoxon signed-rank with Holm correction (nonparametric path)."
            if res["path"] == "friedman"
            else "Paired t-tests with Holm correction (parametric path)."
        )
        lines.append(path_note)
        lines.append("")
        if res["path"] == "friedman" and res["n_sites"] == 4:
            lines.append(
                "> Note: with n=4 site-level observations, the two-sided "
                "Wilcoxon signed-rank floor is p = 0.125; Holm-adjusted p-values "
                "therefore cannot reach α = 0.05 at the pairwise level in this "
                "design regardless of effect size. The Friedman omnibus remains "
                "the primary evidence for a condition effect."
            )
            lines.append("")
    else:
        lines.append("")
    lines.append(_posthoc_block(res))
    lines.append("")
    lines.append(f"**Condition ranking (desc. by mean):** " + " > ".join(res["ranking_desc_by_mean"]))
    lines.append("")
    return "\n".join(lines)


def _a2_section(res):
    lines = ["## A2 — Condition × site interaction (exploratory)", ""]
    lines.append(
        f"Outcome: {res['metric']}. Unit of analysis: {res['unit_of_analysis']} "
        f"(n = {res['n_obs']}). Model: {res['model']}."
    )
    lines.append("")
    lines.append(
        "Primary quantity: the `condition * site` interaction F-test. "
        "Treated as exploratory given the single-subject design and 4 sites."
    )
    lines.append("")
    rows = [
        [
            e["source"],
            _fmt_num(e.get("ss", float("nan"))),
            str(e.get("df", "—")),
            _fmt_num(e.get("ms", float("nan"))),
            _fmt_num(e.get("F", float("nan"))),
            _fmt_p(e.get("p", float("nan"))),
            _fmt_num(e.get("partial_eta_sq", float("nan"))),
        ]
        for e in res["anova_table"]
    ]
    lines.append(_table(
        ["source", "SS", "df", "MS", "F", "p", "partial η²"], rows
    ))
    lines.append("")
    inter = res["interaction_primary"]
    if inter is not None:
        lines.append(
            f"**Interaction test:** F({inter['df']}, "
            f"{res['anova_table'][-1]['df']}) = {_fmt_num(inter['F'])}, "
            f"p = {_fmt_p(inter['p'])}, partial η² = {_fmt_num(inter['partial_eta_sq'])}."
        )
        lines.append("")
    return "\n".join(lines)


def _a3_section(res):
    lines = ["## A3 — Condition-level repeatability (CV)", ""]
    lines.append(
        f"Outcome: per-cell Weber CV (4 × 5 matrix of cell CVs from 12 trials each). "
        f"Mode: **{res['a3_mode']}**."
    )
    lines.append("")
    if res["unstable_cells"]:
        lines.append(
            f"- Unstable cells (|cell mean| < CV_MEAN_EPS = 1e-6): {res['unstable_cells']}"
        )
        lines.append(
            "- Friedman test was not run; per the locked fallback, only "
            "per-condition bootstrap CIs on CV are reported."
        )
    else:
        f = res["friedman"]
        lines.append(
            f"- Friedman Q = {_fmt_num(f['Q'])}, df = {_fmt_num(f['dof'])}, "
            f"p = {_fmt_p(f['p'])}"
        )
        lines.append(f"- Kendall's W: {_fmt_num(f['kendalls_w'])}")
        lines.append(f"- Significant at α = {f['alpha']}: **{f['significant']}**")
    lines.append("")
    lines.append("### Per-condition CV with bootstrap 95% CIs")
    lines.append("")
    rows = [
        [
            b["condition"],
            _fmt_num(b["cv_point_estimate"]),
            _fmt_num(b["ci_low"]),
            _fmt_num(b["ci_high"]),
            str(b["iterations_kept"]),
        ]
        for b in res["bootstrap_cis"]
    ]
    lines.append(_table(
        ["condition", "CV point estimate", "95% CI low", "95% CI high", "iterations kept"],
        rows,
    ))
    lines.append("")
    return "\n".join(lines)


def _michelson_section(michelson_res, comparison):
    lines = [_a1_section("Michelson sensitivity check", michelson_res)]
    lines.append("### Ranking comparison (Weber vs Michelson)")
    lines.append("")
    lines.append(f"- Weber:     {' > '.join(comparison['weber_ranking'])}")
    lines.append(f"- Michelson: {' > '.join(comparison['michelson_ranking'])}")
    lines.append(f"- Rankings match: **{comparison['rankings_match']}**")
    lines.append("")
    return "\n".join(lines)


def _header(params, n_rows):
    lines = [
        "# Formal Analysis Report — NIR Vein-Contrast (Synthetic Data)",
        "",
        "> **Scope.** This report exercises the locked A1 / A2 / A3 / Michelson "
        "sensitivity pipeline on the **single synthetic dataset** generated by "
        "`generator/`. Results are pipeline outputs on simulated data and are "
        "**not** substantive scientific claims about NIR imaging performance. "
        "Real experimental data will be analyzed with the same pipeline on the "
        "completed platform.",
        "",
        f"- Input: `generator/output/synthetic_dataset.csv`",
        f"- Rows (trials): `{n_rows}`",
    ]
    if "subject_id" in params:
        lines.append(f"- Subject: `{params['subject_id']}`")
    if "seed" in params:
        lines.append(f"- Seed: `{params['seed']}`")
    if params.get("sites"):
        lines.append(f"- Sites: {', '.join(params['sites'])}")
    if params.get("conditions"):
        lines.append(f"- Conditions: {', '.join(params['conditions'])}")
    lines.append("")
    return "\n".join(lines)


def _closing(comparison):
    match_phrase = (
        "The Weber and Michelson condition rankings matched, "
        "consistent with the generator's calibrated ordering."
        if comparison["rankings_match"]
        else "The Weber and Michelson rankings **differed** — "
        "condition ordering is metric-dependent in the synthetic data."
    )
    return "\n".join([
        "## Synthetic-data framing",
        "",
        "- Every p-value, F-statistic, effect size, and CI above is a pipeline "
        "output computed on the single synthetic dataset. They do not support "
        "claims about NIR imaging performance on real subjects or platforms.",
        f"- {match_phrase}",
        "- Condition effects reflect the generator's calibrated structure: "
        "baseline intensities `μ_V(λ)` / `μ_B(λ)` at 850/940/1050 nm, small site "
        "offsets `δ_s`, modest condition-by-site interactions `η_{c,s}`, and an "
        "additional registration-noise term `σ_R` on the ratiometric conditions.",
        "- On real data, results will remain specific to the subject, session, "
        "and system, as stated in the original proposal.",
        "",
    ])


def render_report(params, n_rows, assumptions, a1, a2, a3, michelson, comparison):
    sections = [
        _header(params, n_rows),
        _assumption_section(assumptions),
        _a1_section("A1 — Condition main effect (Weber)", a1),
        _a2_section(a2),
        _a3_section(a3),
        _michelson_section(michelson, comparison),
        _closing(comparison),
    ]
    return "\n".join(sections).rstrip() + "\n"
