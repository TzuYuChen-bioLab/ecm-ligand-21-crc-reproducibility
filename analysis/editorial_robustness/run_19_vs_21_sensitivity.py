#!/usr/bin/env python3
"""Reconstruct GSE39582 ECM core-19 and ECM/secreted-factor-21 scores.

The script mirrors the archived analysis rules as closely as possible without R:
GPL570 probes are mapped to unambiguous gene symbols, the highest-IQR probe is
selected using the 443 Discovery tumours only, and every gene is standardized
against that same Discovery reference. The analysis is a membership sensitivity
check; it does not reselect genes, weights, cut-points, or outcomes.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import platform
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-19v21")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import optimize, stats


SCRIPT_DIR = Path(__file__).resolve().parent
RAW = Path(os.environ.get("CBC_19V21_RAW_DIR", SCRIPT_DIR / "raw")).resolve()
OUT = Path(os.environ.get("CBC_19V21_OUTPUT_DIR", SCRIPT_DIR / "results")).resolve()
OUT.mkdir(parents=True, exist_ok=True)

MATRIX_FILE = RAW / "GSE39582_series_matrix.txt.gz"
GPL_FILE = RAW / "GPL570.annot.gz"
HGU_DB_FILE = RAW / "hgu133plus2_db" / "hgu133plus2.db" / "inst" / "extdata" / "hgu133plus2.sqlite"
ORG_HS_DB_FILE = RAW / "org_hs_db" / "org.Hs.eg.db" / "inst" / "extdata" / "org.Hs.eg.sqlite"
CLINICAL_FILE = Path(
    os.environ.get(
        "CBC_19V21_CLINICAL_FILE",
        SCRIPT_DIR / "resources" / "GSE39582_tumour_only_RFS_analysis_ready.csv",
    )
).resolve()
_default_provenance = (
    SCRIPT_DIR.parent
    / "mechanistic_triangulation_v1"
    / "resources"
    / "CRC_ECM_LIGAND_21_provenance_annotation.csv"
)
if not _default_provenance.exists():
    _default_provenance = SCRIPT_DIR / "resources" / "CRC_ECM_LIGAND_21_provenance_annotation.csv"
PROVENANCE_FILE = Path(
    os.environ.get("CBC_19V21_PROVENANCE_FILE", _default_provenance)
).resolve()

CORE_19 = [
    "COL1A1", "COL1A2", "COL4A1", "COL4A2", "COL4A5", "COL6A1", "COL6A2", "COL6A3",
    "FN1", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMC1", "THBS1", "THBS2", "HSPG2",
    "TNC", "TNXB",
]
ADDED_2 = ["PTN", "MDK"]
FULL_21 = CORE_19 + ADDED_2


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_bioconductor_annotation() -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Mirror AnnotationDbi::select(hgu133plus2.db, columns='SYMBOL')."""
    probes_by_gene: dict[str, list[str]] = defaultdict(list)
    audit_rows: list[dict[str, object]] = []
    with sqlite3.connect(ORG_HS_DB_FILE) as org_conn:
        gene_to_symbol = dict(
            org_conn.execute(
                "SELECT genes.gene_id, gene_info.symbol FROM genes JOIN gene_info USING(_id)"
            ).fetchall()
        )
    probe_symbols: dict[str, set[str]] = defaultdict(set)
    probe_multiplicity: dict[str, int] = {}
    with sqlite3.connect(HGU_DB_FILE) as hgu_conn:
        for probe, gene_id, is_multiple in hgu_conn.execute(
            "SELECT probe_id, gene_id, is_multiple FROM probes"
        ):
            symbol = gene_to_symbol.get(gene_id)
            if symbol:
                probe_symbols[probe].add(symbol)
            probe_multiplicity[probe] = max(probe_multiplicity.get(probe, 0), int(is_multiple))
    for probe, symbols in probe_symbols.items():
        unique_symbols = sorted(symbols)
        if len(unique_symbols) != 1:
            continue
        symbol = unique_symbols[0]
        if symbol in FULL_21:
            probes_by_gene[symbol].append(probe)
            audit_rows.append(
                {
                    "probe": probe,
                    "gene": symbol,
                    "annotation_unambiguous": True,
                    "database_is_multiple_flag": bool(probe_multiplicity.get(probe, 0)),
                    "annotation_source": "Bioconductor hgu133plus2.db 3.13.0 joined to org.Hs.eg.db 3.13.0",
                }
            )
    missing = [g for g in FULL_21 if not probes_by_gene[g]]
    if missing:
        raise RuntimeError(f"No unambiguous hgu133plus2.db probe for: {missing}")
    return dict(probes_by_gene), pd.DataFrame(audit_rows)


def parse_series_matrix(candidate_probes: set[str]):
    sample_ids: list[str] | None = None
    dataset_labels: list[str] | None = None
    expression_header: list[str] | None = None
    probe_values: dict[str, np.ndarray] = {}
    in_table = False
    with gzip.open(MATRIX_FILE, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                row = next(csv.reader([line], delimiter="\t"))
                sample_ids = [x.strip().strip('"') for x in row[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "dataset:" in line.lower():
                row = next(csv.reader([line], delimiter="\t"))
                dataset_labels = [x.strip().strip('"').split(":", 1)[1].strip() for x in row[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                row = next(csv.reader([line], delimiter="\t"))
                if expression_header is None:
                    expression_header = [x.strip().strip('"') for x in row]
                    continue
                probe = row[0].strip().strip('"')
                if probe in candidate_probes:
                    probe_values[probe] = np.asarray(
                        [float(x.strip().strip('"')) if x.strip().strip('"') not in {"", "null", "NA"} else np.nan for x in row[1:]],
                        dtype=float,
                    )
    if sample_ids is None or dataset_labels is None or expression_header is None:
        raise RuntimeError("Required GSE39582 matrix metadata were not found")
    matrix_samples = expression_header[1:]
    if sample_ids != matrix_samples:
        raise RuntimeError("Sample order differs between GEO metadata and expression table")
    if len(dataset_labels) != len(sample_ids):
        raise RuntimeError("Dataset label count does not match sample count")
    missing_probes = sorted(candidate_probes - set(probe_values))
    if missing_probes:
        raise RuntimeError(f"Candidate probes missing from series matrix: {missing_probes}")
    return sample_ids, dataset_labels, probe_values


def select_reference_probes(
    probes_by_gene: dict[str, list[str]],
    probe_values: dict[str, np.ndarray],
    discovery_mask: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, str]]:
    rows: list[dict[str, object]] = []
    selected: dict[str, str] = {}
    for gene in FULL_21:
        candidates = []
        for probe in probes_by_gene[gene]:
            values = probe_values[probe][discovery_mask]
            q25, q75 = np.nanpercentile(values, [25, 75], method="linear")
            iqr = float(q75 - q25)
            candidates.append((probe, iqr))
        candidates.sort(key=lambda x: (-x[1], x[0]))
        selected_probe = candidates[0][0]
        selected[gene] = selected_probe
        for rank, (probe, iqr) in enumerate(candidates, start=1):
            rows.append(
                {
                    "gene": gene,
                    "probe": probe,
                    "discovery_reference_iqr": iqr,
                    "iqr_rank_within_gene": rank,
                    "selected": probe == selected_probe,
                    "selection_rule": "highest IQR in 443 Discovery tumours; probe ID ascending for ties",
                }
            )
    return pd.DataFrame(rows), selected


def build_scores(
    sample_ids: list[str],
    dataset_labels: list[str],
    selected: dict[str, str],
    probe_values: dict[str, np.ndarray],
):
    labels = np.asarray(dataset_labels, dtype=object)
    discovery_mask = np.char.lower(labels.astype(str)) == "discovery"
    tumour_mask = np.isin(np.char.lower(labels.astype(str)), ["discovery", "validation"])
    gene_expr = np.vstack([probe_values[selected[g]] for g in FULL_21])
    centres = np.nanmean(gene_expr[:, discovery_mask], axis=1)
    scales = np.nanstd(gene_expr[:, discovery_mask], axis=1, ddof=1)
    if np.any(~np.isfinite(scales)) or np.any(scales <= 0):
        raise RuntimeError("Invalid Discovery-reference gene standard deviation")
    z = (gene_expr - centres[:, None]) / scales[:, None]
    score19 = np.nanmean(z[: len(CORE_19), :], axis=0)
    score2 = np.nanmean(z[len(CORE_19) :, :], axis=0)
    score21 = np.nanmean(z, axis=0)
    ref19 = score19[discovery_mask]
    ref21 = score21[discovery_mask]
    score19_z = (score19 - np.nanmean(ref19)) / np.nanstd(ref19, ddof=1)
    score21_z = (score21 - np.nanmean(ref21)) / np.nanstd(ref21, ddof=1)
    sample_df = pd.DataFrame(
        {
            "sample": sample_ids,
            "cohort_split_matrix": labels,
            "score19_raw": score19,
            "score21_raw_reconstructed": score21,
            "score19_z_discovery_reference": score19_z,
            "score21_z_discovery_reference": score21_z,
            "PTN_MDK_2_raw": score2,
        }
    )
    sample_df = sample_df.loc[tumour_mask].reset_index(drop=True)
    expression_df = pd.DataFrame(
        gene_expr[:, tumour_mask].T,
        columns=FULL_21,
    )
    expression_df.insert(0, "cohort_split", labels[tumour_mask])
    expression_df.insert(0, "sample", np.asarray(sample_ids, dtype=object)[tumour_mask])
    params = pd.DataFrame(
        {
            "gene": FULL_21,
            "subset": ["core_19"] * 19 + ["added_secreted_factor"] * 2,
            "selected_probe": [selected[g] for g in FULL_21],
            "discovery_centre": centres,
            "discovery_scale": scales,
        }
    )
    return sample_df, params, expression_df


def cox_efron(time, event, x, names):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    x = np.asarray(x, dtype=float)
    valid = np.isfinite(time) & (time > 0) & np.isin(event, [0, 1]) & np.all(np.isfinite(x), axis=1)
    time, event, x = time[valid], event[valid], x[valid]
    n, p = x.shape
    event_times = np.unique(time[event == 1])

    def objective(beta):
        xb = np.clip(x @ beta, -40, 40)
        risk = np.exp(xb)
        loglik = 0.0
        grad = np.zeros(p)
        hess = np.zeros((p, p))
        for t in event_times:
            dmask = (time == t) & (event == 1)
            rmask = time >= t
            d = int(dmask.sum())
            xd = x[dmask]
            rd = risk[dmask]
            xr = x[rmask]
            rr = risk[rmask]
            s0 = rr.sum()
            s1 = (rr[:, None] * xr).sum(axis=0)
            s2 = np.einsum("i,ij,ik->jk", rr, xr, xr)
            e0 = rd.sum()
            e1 = (rd[:, None] * xd).sum(axis=0)
            e2 = np.einsum("i,ij,ik->jk", rd, xd, xd)
            loglik += xb[dmask].sum()
            grad += xd.sum(axis=0)
            for ell in range(d):
                frac = ell / d
                den = s0 - frac * e0
                num1 = s1 - frac * e1
                num2 = s2 - frac * e2
                loglik -= math.log(den)
                grad -= num1 / den
                hess -= num2 / den - np.outer(num1, num1) / (den * den)
        return -loglik, -grad, -hess

    fit = optimize.minimize(
        lambda b: objective(b)[0],
        np.zeros(p),
        jac=lambda b: objective(b)[1],
        method="BFGS",
        options={"gtol": 1e-9, "maxiter": 1000},
    )
    beta = fit.x
    information = objective(beta)[2]
    covariance = np.linalg.pinv(information)
    se = np.sqrt(np.diag(covariance))
    z = beta / se
    pvalue = 2 * stats.norm.sf(np.abs(z))
    result = pd.DataFrame(
        {
            "term": names,
            "beta": beta,
            "se": se,
            "HR": np.exp(beta),
            "CI_low": np.exp(beta - 1.96 * se),
            "CI_high": np.exp(beta + 1.96 * se),
            "Wald_z": z,
            "P": pvalue,
        }
    )
    result.attrs.update(
        n=n,
        events=int(event.sum()),
        converged=bool(fit.success or np.linalg.norm(objective(beta)[1]) < 1e-5),
        optimizer_message=str(fit.message),
        log_likelihood=float(-objective(beta)[0]),
    )
    return result, beta, valid


def harrell_c_index(time, event, risk_score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk_score = np.asarray(risk_score, dtype=float)
    concordant = tied = comparable = 0.0
    n = len(time)
    for i in range(n):
        for j in range(i + 1, n):
            if time[i] == time[j]:
                continue
            if time[i] < time[j] and event[i] == 1:
                shorter, longer = i, j
            elif time[j] < time[i] and event[j] == 1:
                shorter, longer = j, i
            else:
                continue
            comparable += 1
            if risk_score[shorter] > risk_score[longer]:
                concordant += 1
            elif risk_score[shorter] == risk_score[longer]:
                tied += 1
    return (concordant + 0.5 * tied) / comparable if comparable else np.nan


def run_survival(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort in ["Discovery", "Validation"]:
        dat = merged.loc[merged["cohort_split"] == cohort].copy()
        for signature, score_col in [
            ("ECM_CORE_19", "score19_z_discovery_reference"),
            ("ECM_LIGAND_21", "score21_z_discovery_reference"),
        ]:
            fit, beta, valid = cox_efron(
                dat["time"].to_numpy(),
                dat["event"].to_numpy(),
                dat[[score_col]].to_numpy(),
                ["score_per_1_SD"],
            )
            target = fit.iloc[0]
            cindex = harrell_c_index(
                dat.loc[valid, "time"],
                dat.loc[valid, "event"],
                dat.loc[valid, score_col].to_numpy() * beta[0],
            )
            rows.append(
                {
                    "cohort": cohort,
                    "model": "univariable Cox; continuous score per Discovery-reference SD",
                    "signature": signature,
                    "n": fit.attrs["n"],
                    "events": fit.attrs["events"],
                    "HR": target.HR,
                    "CI_low": target.CI_low,
                    "CI_high": target.CI_high,
                    "Wald_z": target.Wald_z,
                    "P": target.P,
                    "Harrell_C": cindex,
                    "converged": fit.attrs["converged"],
                    "interpretation_guardrail": "Sensitivity comparison only; no outcome-based gene selection or weighting.",
                }
            )
    return pd.DataFrame(rows)


def make_origin_audit() -> pd.DataFrame:
    provenance = pd.read_csv(PROVENANCE_FILE, encoding="utf-8-sig")
    provenance["membership_layer"] = np.where(
        provenance["gene"].isin(CORE_19), "ECM_CORE_19", "PTN_MDK_added_2"
    )
    provenance["operational_reason_for_21_gene_scope"] = np.where(
        provenance["gene"].isin(CORE_19),
        "Structural core-matrisome component",
        "Matrisome-associated secreted growth factor extending the panel from ECM structure to ECM-linked paracrine signalling",
    )
    provenance["historical_selection_status"] = (
        "Exact original inclusion/exclusion history is not recoverable; classification is a retrospective annotation of the unchanged frozen panel."
    )
    provenance["outcome_used_for_membership"] = False
    provenance["current_data_used_to_reweight"] = False
    return provenance


def make_figure(origin: pd.DataFrame, merged: pd.DataFrame, cox: pd.DataFrame, agreement: dict):
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9})
    fig = plt.figure(figsize=(10.5, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.20])
    fig.subplots_adjust(left=0.145, right=0.97, top=0.90, bottom=0.13, hspace=0.42, wspace=0.50)

    ax_a = fig.add_subplot(gs[0, :])
    category_order = ["Collagens", "ECM glycoproteins", "Proteoglycans", "Secreted factors"]
    colors = {
        "Collagens": "#3B6FB6",
        "ECM glycoproteins": "#67A9CF",
        "Proteoglycans": "#9ECAE1",
        "Secreted factors": "#D95F02",
    }
    y_positions = {cat: i for i, cat in enumerate(category_order[::-1])}
    for cat in category_order:
        genes = origin.loc[origin["matrisome_category"] == cat, "gene"].tolist()
        y = y_positions[cat]
        x = np.arange(len(genes)) * 1.06
        for xi, gene in zip(x, genes):
            ax_a.text(
                xi, y, gene, ha="center", va="center", fontsize=7.2,
                color="white" if cat in {"Collagens", "Secreted factors"} else "#17324D",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.35", facecolor=colors[cat], edgecolor="white", linewidth=0.8),
            )
    ax_a.set_yticks([y_positions[c] for c in category_order], category_order)
    ax_a.set_xlim(-1.25, 10.8)
    ax_a.set_ylim(-0.7, 3.7)
    ax_a.set_xticks([])
    ax_a.set_title("A  Retrospective functional decomposition of the unchanged 21-gene panel", loc="left", fontweight="bold")
    ax_a.text(
        10.7, y_positions["Secreted factors"],
        "PTN and MDK extend\n19 structural core-matrisome genes\nto ECM-linked secreted signalling",
        ha="right", va="center", fontsize=8.1, color="#8C2D04",
    )
    for spine in ax_a.spines.values():
        spine.set_visible(False)

    ax_b = fig.add_subplot(gs[1, 0])
    palette = {"Discovery": "#2166AC", "Validation": "#B2182B"}
    for cohort in ["Discovery", "Validation"]:
        d = merged.loc[merged["cohort_split"] == cohort]
        ax_b.scatter(
            d["score19_z_discovery_reference"], d["score21_z_discovery_reference"],
            s=14, alpha=0.58, label=cohort, color=palette[cohort], linewidth=0,
        )
    lims = [
        min(merged["score19_z_discovery_reference"].min(), merged["score21_z_discovery_reference"].min()),
        max(merged["score19_z_discovery_reference"].max(), merged["score21_z_discovery_reference"].max()),
    ]
    ax_b.plot(lims, lims, linestyle="--", color="#555555", linewidth=1)
    ax_b.set_xlim(lims)
    ax_b.set_ylim(lims)
    ax_b.set_xlabel("ECM_CORE_19 score (Discovery-reference SD)")
    ax_b.set_ylabel("ECM_LIGAND_21 score\n(Discovery-reference SD)")
    ax_b.set_title("B  Score agreement", loc="left", fontweight="bold")
    ax_b.text(
        0.03, 0.97,
        f"Pearson r = {agreement['pearson_r']:.3f}\nRisk-group agreement = {agreement['risk_group_agreement']*100:.1f}%",
        transform=ax_b.transAxes, ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#BBBBBB"),
    )
    ax_b.legend(frameon=False, loc="lower right", fontsize=8)
    ax_b.grid(color="#E6E6E6", linewidth=0.6)

    ax_c = fig.add_subplot(gs[1, 1])
    plot_rows = []
    for cohort_i, cohort in enumerate(["Discovery", "Validation"]):
        for sig_i, signature in enumerate(["ECM_CORE_19", "ECM_LIGAND_21"]):
            row = cox[(cox["cohort"] == cohort) & (cox["signature"] == signature)].iloc[0]
            y = 3 - cohort_i * 2 - sig_i * 0.62
            plot_rows.append((y, cohort, signature, row))
    for y, cohort, signature, row in plot_rows:
        color = "#4D4D4D" if signature == "ECM_CORE_19" else "#D95F02"
        marker = "o" if signature == "ECM_CORE_19" else "s"
        ax_c.errorbar(
            row.HR, y,
            xerr=np.array([[row.HR - row.CI_low], [row.CI_high - row.HR]]),
            fmt=marker, color=color, ecolor=color, capsize=3, markersize=6,
        )
    ax_c.axvline(1, color="#666666", linestyle="--", linewidth=1)
    ax_c.set_xscale("log")
    ax_c.set_xlim(0.62, 2.30)
    compact_labels = [
        f"{cohort} | {'19-gene' if signature == 'ECM_CORE_19' else '21-gene'}"
        for _, cohort, signature, _ in plot_rows
    ]
    ax_c.set_yticks([x[0] for x in plot_rows], compact_labels)
    ax_c.tick_params(axis="y", labelsize=8)
    ax_c.set_xlabel("Hazard ratio per 1 SD (95% CI)")
    ax_c.set_title("C  Recurrence-free survival sensitivity", loc="left", fontweight="bold")
    ax_c.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    for y, _, _, row in plot_rows:
        ax_c.text(
            0.995, y,
            f"{row.HR:.2f} ({row.CI_low:.2f}-{row.CI_high:.2f})",
            transform=ax_c.get_yaxis_transform(), va="center", ha="right", fontsize=7.3,
            bbox=dict(facecolor="white", edgecolor="none", alpha=1.0, pad=1.5),
        )

    fig.suptitle(
        "Supplementary Figure S15. Functional origin and 19-vs-21 membership sensitivity",
        fontsize=12.5, fontweight="bold",
    )
    fig.text(
        0.09, 0.035,
        "Annotation describes biological scope; it does not reconstruct an undocumented historical selection rule.\n"
        "Cox models are univariable continuous-score sensitivities; similarity does not establish biological necessity or causality.",
        ha="left", va="bottom", fontsize=7.2, color="#444444",
    )
    fig.savefig(OUT / "Supplementary_Figure_S15_19_vs_21_sensitivity.png", dpi=450)
    fig.savefig(OUT / "Supplementary_Figure_S15_19_vs_21_sensitivity.pdf")
    plt.close(fig)


def main():
    required = [MATRIX_FILE, GPL_FILE, HGU_DB_FILE, ORG_HS_DB_FILE, CLINICAL_FILE, PROVENANCE_FILE]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")

    probes_by_gene, annotation_candidates = parse_bioconductor_annotation()
    candidate_probes = {p for probes in probes_by_gene.values() for p in probes}
    sample_ids, dataset_labels, probe_values = parse_series_matrix(candidate_probes)
    labels = np.asarray(dataset_labels, dtype=object)
    discovery_mask = np.char.lower(labels.astype(str)) == "discovery"
    label_counts = Counter(dataset_labels)
    if int(discovery_mask.sum()) != 443 or label_counts.get("validation", 0) != 123:
        raise RuntimeError(f"Unexpected cohort labels: {label_counts}")
    probe_selection, selected = select_reference_probes(probes_by_gene, probe_values, discovery_mask)
    sample_scores, score_parameters, selected_expression = build_scores(
        sample_ids, dataset_labels, selected, probe_values
    )

    clinical = pd.read_csv(CLINICAL_FILE)
    clinical["sample"] = clinical["sample"].astype(str)
    merged = clinical.merge(sample_scores, on="sample", how="left", validate="one_to_one")
    if merged["score19_raw"].isna().any():
        raise RuntimeError("Some RFS samples were absent from the public expression matrix")
    split_match = merged["cohort_split"].str.lower() == merged["cohort_split_matrix"].str.lower()
    if not split_match.all():
        raise RuntimeError("Cohort labels differ between clinical and expression files")

    delta = merged["score21_raw_reconstructed"] - merged["score_raw"]
    reconstruction_r = stats.pearsonr(merged["score21_raw_reconstructed"], merged["score_raw"]).statistic
    reconstruction_slope, reconstruction_intercept = np.polyfit(
        merged["score21_raw_reconstructed"], merged["score_raw"], 1
    )
    agreement_r = stats.pearsonr(
        merged["score19_z_discovery_reference"], merged["score21_z_discovery_reference"]
    ).statistic
    agreement_spearman = stats.spearmanr(
        merged["score19_z_discovery_reference"], merged["score21_z_discovery_reference"]
    ).statistic

    discovery = merged[merged["cohort_split"] == "Discovery"]
    cut19 = float(discovery["score19_raw"].median())
    cut21 = float(discovery["score21_raw_reconstructed"].median())
    merged["risk_group_19_discovery_median"] = np.where(merged["score19_raw"] >= cut19, "High", "Low")
    merged["risk_group_21_discovery_median"] = np.where(
        merged["score21_raw_reconstructed"] >= cut21, "High", "Low"
    )
    risk_agreement = float(
        (merged["risk_group_19_discovery_median"] == merged["risk_group_21_discovery_median"]).mean()
    )
    kappa_table = pd.crosstab(
        merged["risk_group_19_discovery_median"], merged["risk_group_21_discovery_median"]
    ).reindex(index=["Low", "High"], columns=["Low", "High"], fill_value=0)
    po = np.trace(kappa_table.to_numpy()) / kappa_table.to_numpy().sum()
    row_marg = kappa_table.sum(axis=1).to_numpy() / kappa_table.to_numpy().sum()
    col_marg = kappa_table.sum(axis=0).to_numpy() / kappa_table.to_numpy().sum()
    pe = float((row_marg * col_marg).sum())
    kappa = float((po - pe) / (1 - pe)) if pe < 1 else np.nan

    cox = run_survival(merged)
    origin = make_origin_audit()
    core_definition = origin.loc[origin["gene"].isin(CORE_19)].copy()
    core_definition.insert(0, "signature", "ECM_CORE_19")
    core_definition["weight_within_19_gene_score"] = 1 / 19

    agreement = {
        "pearson_r": float(agreement_r),
        "spearman_rho": float(agreement_spearman),
        "risk_group_agreement": risk_agreement,
        "cohen_kappa": kappa,
        "discovery_median_cutpoint_19": cut19,
        "discovery_median_cutpoint_21_reconstructed": cut21,
        "n_analysis_samples": int(len(merged)),
        "n_discovery_analysis_samples": int((merged["cohort_split"] == "Discovery").sum()),
        "n_validation_analysis_samples": int((merged["cohort_split"] == "Validation").sum()),
    }

    summary_rows = [
        ("Pearson correlation, score19 vs score21", agreement_r),
        ("Spearman correlation, score19 vs score21", agreement_spearman),
        ("Median absolute standardized score difference", float(np.median(np.abs(merged["score21_z_discovery_reference"] - merged["score19_z_discovery_reference"])))),
        ("95th percentile absolute standardized score difference", float(np.quantile(np.abs(merged["score21_z_discovery_reference"] - merged["score19_z_discovery_reference"]), 0.95))),
        ("Discovery-median high/low group agreement", risk_agreement),
        ("Cohen kappa for high/low group agreement", kappa),
        ("Reconstructed score21 vs archived score Pearson correlation", float(reconstruction_r)),
        ("Reconstructed-to-archived score slope", float(reconstruction_slope)),
        ("Reconstructed-to-archived score intercept", float(reconstruction_intercept)),
        ("Maximum absolute reconstructed-vs-archived score difference", float(np.max(np.abs(delta)))),
    ]
    agreement_df = pd.DataFrame(summary_rows, columns=["metric", "value"])

    probe_selection.to_csv(OUT / "GSE39582_19_vs_21_probe_selection.csv", index=False)
    annotation_candidates.to_csv(OUT / "Bioconductor_21_gene_candidate_probes.csv", index=False)
    score_parameters.to_csv(OUT / "GSE39582_19_vs_21_standardisation_parameters.csv", index=False)
    selected_expression.to_csv(OUT / "GSE39582_selected_21_gene_expression.csv", index=False)
    merged.to_csv(OUT / "GSE39582_19_vs_21_sample_scores_and_RFS.csv", index=False)
    cox.to_csv(OUT / "GSE39582_19_vs_21_Cox_results.csv", index=False)
    agreement_df.to_csv(OUT / "GSE39582_19_vs_21_agreement_summary.csv", index=False)
    kappa_table.to_csv(OUT / "GSE39582_19_vs_21_risk_group_cross_tab.csv")
    origin.to_csv(OUT / "ECM_LIGAND_21_origin_audit.csv", index=False, encoding="utf-8-sig")
    core_definition.to_csv(OUT / "ECM_CORE_19_definition.csv", index=False, encoding="utf-8-sig")

    audit = {
        "analysis_role": "Prespecified membership sensitivity: 19 structural core-matrisome genes versus the unchanged 21-gene panel adding PTN and MDK",
        "claim_guardrail": "This analysis tests dependence on PTN/MDK; it does not prove that the two genes were historically selected by current data, improve prediction, are necessary, or are causal.",
        "matrix_source": "NCBI GEO GSE39582 Series Matrix",
        "platform_annotation_source": "Bioconductor hgu133plus2.db 3.13.0 with org.Hs.eg.db 3.13.0; NCBI GPL570 retained as platform cross-reference",
        "matrix_sha256": sha256(MATRIX_FILE),
        "gpl_annotation_sha256": sha256(GPL_FILE),
        "hgu133plus2_db_sha256": sha256(HGU_DB_FILE),
        "org_hs_eg_db_sha256": sha256(ORG_HS_DB_FILE),
        "clinical_sha256": sha256(CLINICAL_FILE),
        "sample_counts_matrix": dict(label_counts),
        "selected_probes": selected,
        "score21_reconstruction": {
            "pearson_r": float(reconstruction_r),
            "slope": float(reconstruction_slope),
            "intercept": float(reconstruction_intercept),
            "max_absolute_difference": float(np.max(np.abs(delta))),
        },
        "score19_vs_score21": agreement,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": plt.matplotlib.__version__,
        },
    }
    (OUT / "GSE39582_19_vs_21_analysis_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    make_figure(origin, merged, cox, agreement)

    print(json.dumps({"agreement": agreement, "reconstruction": audit["score21_reconstruction"], "cox": cox.to_dict(orient="records")}, indent=2))


if __name__ == "__main__":
    main()
