#!/usr/bin/env python3
"""Conservative robustness analyses requested during editorial review."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


def paired_exact_sensitivities(patient_scores: Path) -> pd.DataFrame:
    data = pd.read_csv(patient_scores)
    rows: list[dict[str, object]] = []
    for (cell_class, comparison), part in data.groupby(
        ["cell_class", "comparison"], sort=False
    ):
        case_condition = comparison
        wide = part.pivot(index="patient", columns="condition", values="score")
        if "Normal" not in wide or case_condition not in wide:
            continue
        difference = (wide[case_condition] - wide["Normal"]).dropna().to_numpy(float)
        nonzero = difference[difference != 0]
        positive = int(np.sum(nonzero > 0))
        n = int(nonzero.size)
        sign_p = float(stats.binomtest(positive, n, 0.5).pvalue) if n else math.nan
        if n:
            wilcoxon = stats.wilcoxon(
                nonzero,
                zero_method="wilcox",
                alternative="two-sided",
                method="exact",
            )
            wilcoxon_p = float(wilcoxon.pvalue)
            absolute = np.abs(nonzero)
            observed = abs(float(np.mean(nonzero)))
            permuted = []
            for signs in itertools.product((-1.0, 1.0), repeat=n):
                permuted.append(abs(float(np.mean(absolute * np.asarray(signs)))))
            sign_flip_p = float(np.mean(np.asarray(permuted) >= observed - 1e-12))
        else:
            wilcoxon_p = math.nan
            sign_flip_p = math.nan
        t_test = stats.ttest_1samp(difference, 0.0) if difference.size > 1 else None
        rows.append(
            {
                "cell_class": cell_class,
                "comparison": f"{comparison} minus Normal",
                "n_pairs": int(difference.size),
                "mean_difference": float(np.mean(difference)),
                "median_difference": float(np.median(difference)),
                "positive_pairs": positive,
                "sign_concordance": positive / n if n else math.nan,
                "paired_t_two_sided_p": float(t_test.pvalue) if t_test else math.nan,
                "exact_sign_two_sided_p": sign_p,
                "exact_wilcoxon_two_sided_p": wilcoxon_p,
                "exact_sign_flip_mean_two_sided_p": sign_flip_p,
                "interpretation": (
                    "Directionally concordant; exact small-sample inference is primary for robustness."
                    if positive == n
                    else "Mixed paired directions; exact small-sample inference is non-confirmatory."
                ),
            }
        )
    return pd.DataFrame(rows)


def modified_hartung_knapp(
    effects: pd.DataFrame, label: str, cohorts: list[str]
) -> dict[str, object]:
    sub = effects.set_index("cohort").loc[cohorts].reset_index()
    yi = sub["logHR"].to_numpy(float)
    sei = sub["SE_logHR"].to_numpy(float)
    k = int(yi.size)
    # The supplied REML analysis estimated tau^2 = 0 for both three-cohort sets.
    tau2 = 0.0
    weight = 1.0 / (sei**2 + tau2)
    mean = float(np.sum(weight * yi) / np.sum(weight))
    q_residual = float(np.sum(weight * (yi - mean) ** 2))
    hk_scale = q_residual / (k - 1)
    modified_scale = max(1.0, hk_scale)
    se = math.sqrt(modified_scale / float(np.sum(weight)))
    df = k - 1
    critical = float(stats.t.ppf(0.975, df))
    low = mean - critical * se
    high = mean + critical * se
    statistic = mean / se
    p_value = float(2.0 * stats.t.sf(abs(statistic), df))
    return {
        "analysis": label,
        "cohorts": ";".join(cohorts),
        "k": k,
        "tau2_REML": tau2,
        "pooled_HR": math.exp(mean),
        "modified_HK_CI_low": math.exp(low),
        "modified_HK_CI_high": math.exp(high),
        "modified_HK_p": p_value,
        "df": df,
        "Q_residual": q_residual,
        "unmodified_HK_scale": hk_scale,
        "modified_HK_scale": modified_scale,
        "prediction_low_same_tau0_model": math.exp(low),
        "prediction_high_same_tau0_model": math.exp(high),
        "note": (
            "Modified Hartung-Knapp bounds the residual scale at 1; with tau^2=0, "
            "the same t-based interval is the model-based prediction interval."
        ),
    }


def meta_sensitivities(effect_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    effects = pd.read_csv(effect_file)
    primary = ["GSE17536", "GSE17537", "GSE33113"]
    swap = ["GSE14333", "GSE17537", "GSE33113"]
    summaries = pd.DataFrame(
        [
            modified_hartung_knapp(
                effects,
                "Primary: GSE17536 retained; GSE14333 excluded for overlap",
                primary,
            ),
            modified_hartung_knapp(
                effects,
                "Sensitivity: GSE14333 replaces GSE17536",
                swap,
            ),
        ]
    )
    loo_rows: list[dict[str, object]] = []
    for omitted in primary:
        retained = [cohort for cohort in primary if cohort != omitted]
        row = modified_hartung_knapp(
            effects, f"Primary leave-one-out: omit {omitted}", retained
        )
        row["omitted_cohort"] = omitted
        loo_rows.append(row)
    return summaries, pd.DataFrame(loo_rows)


def gridded_shift_result(
    part: pd.DataFrame, section: str, analysis: str, grid_size: int
) -> dict[str, object]:
    x = part["x"].to_numpy(float)
    y = part["y"].to_numpy(float)
    x_scaled = (x - x.min()) / max(float(x.max() - x.min()), np.finfo(float).eps)
    y_scaled = (y - y.min()) / max(float(y.max() - y.min()), np.finfo(float).eps)
    ix = np.minimum((x_scaled * grid_size).astype(int), grid_size - 1)
    iy = np.minimum((y_scaled * grid_size).astype(int), grid_size - 1)
    grouped = (
        part.assign(ix=ix, iy=iy)
        .groupby(["iy", "ix"], as_index=False)
        .agg(
            epithelial_target_score=("epithelial_target_score", "mean"),
            nearest_fibroblast_ECM_score=("nearest_fibroblast_ECM_score", "mean"),
            n_epithelial_bins=("section", "size"),
        )
    )
    target = np.full((grid_size, grid_size), np.nan)
    ecm = np.full((grid_size, grid_size), np.nan)
    target[grouped["iy"].to_numpy(int), grouped["ix"].to_numpy(int)] = grouped[
        "epithelial_target_score"
    ]
    ecm[grouped["iy"].to_numpy(int), grouped["ix"].to_numpy(int)] = grouped[
        "nearest_fibroblast_ECM_score"
    ]
    observed_mask = np.isfinite(target) & np.isfinite(ecm)
    observed = float(stats.spearmanr(target[observed_mask], ecm[observed_mask]).statistic)
    null: list[float] = []
    valid_counts: list[int] = []
    for dy in range(grid_size):
        for dx in range(grid_size):
            if dx == 0 and dy == 0:
                continue
            shifted = np.roll(ecm, shift=(dy, dx), axis=(0, 1))
            mask = np.isfinite(target) & np.isfinite(shifted)
            if int(mask.sum()) < 10:
                continue
            rho = stats.spearmanr(target[mask], shifted[mask]).statistic
            if np.isfinite(rho):
                null.append(float(rho))
                valid_counts.append(int(mask.sum()))
    null_array = np.asarray(null)
    one_sided_p = (1.0 + float(np.sum(null_array >= observed - 1e-12))) / (
        1.0 + float(null_array.size)
    )
    return {
        "analysis": analysis,
        "section": section,
        "grid_size": grid_size,
        "occupied_blocks": int(observed_mask.sum()),
        "median_epithelial_bins_per_block": float(grouped["n_epithelial_bins"].median()),
        "observed_block_spearman_rho": observed,
        "valid_nonzero_toroidal_shifts": int(null_array.size),
        "minimum_overlap_blocks_across_shifts": int(min(valid_counts)),
        "one_sided_add_one_shift_p": one_sided_p,
        "null_rho_median": float(np.median(null_array)),
        "null_rho_95pct_low": float(np.quantile(null_array, 0.025)),
        "null_rho_95pct_high": float(np.quantile(null_array, 0.975)),
        "limitation": (
            "Coarse-grid toroidal shifts preserve gridded autocorrelation but wrap tissue boundaries; "
            "this is a sensitivity analysis, not patient-level inference."
        ),
    }


def spatial_shift_sensitivities(figure_root: Path) -> pd.DataFrame:
    spatial_dir = figure_root / "Supplementary_Figure_S9"
    inputs = [
        (
            "Primary broad definition",
            spatial_dir / "GSE280315_primary_spatial_neighbour_source_data.csv",
        ),
        (
            "Strict exclusion of pericytes",
            spatial_dir / "GSE280315_strict_no_pericytes_spatial_neighbour_source_data.csv",
        ),
    ]
    rows: list[dict[str, object]] = []
    for analysis, path in inputs:
        data = pd.read_csv(path)
        for section, part in data.groupby("section", sort=True):
            for grid_size in (8, 12, 16):
                rows.append(
                    gridded_shift_result(part, str(section), analysis, grid_size)
                )
    return pd.DataFrame(rows)


def spatial_specimen_summary(figure_root: Path) -> pd.DataFrame:
    source = (
        figure_root
        / "Supplementary_Figure_S9"
        / "primary_vs_strict_no_pericytes_comparison.csv"
    )
    data = pd.read_csv(source)
    rows: list[dict[str, object]] = []
    for analysis, part in data.groupby("analysis", sort=False):
        positive = int(np.sum(part["spearman_rho"].to_numpy(float) > 0))
        n = int(part.shape[0])
        rows.append(
            {
                "analysis": analysis,
                "n_specimens": n,
                "positive_specimens": positive,
                "median_bin_level_rho": float(part["spearman_rho"].median()),
                "hypothetical_one_sided_sign_p_if_independent": float(
                    stats.binomtest(positive, n, 0.5, alternative="greater").pvalue
                ),
                "inference_status": (
                    "Descriptive only: specimen independence was not established; "
                    "the sign-test P value is an illustrative upper-level check, not a valid patient-level test."
                ),
            }
        )
    return pd.DataFrame(rows)


def tgfb_sign_sensitivity(figure_root: Path) -> pd.DataFrame:
    source = (
        figure_root
        / "Supplementary_Figure_S10"
        / "ECM_LIGAND_21_patient_paired_differences.csv"
    )
    data = pd.read_csv(source)
    delta = data["delta"].to_numpy(float)
    decreases = int(np.sum(delta < 0))
    n = int(delta.size)
    return pd.DataFrame(
        [
            {
                "analysis": "Primary TGF-beta blockade minus control",
                "n_pairs": n,
                "patients_with_decrease": decreases,
                "mean_difference": float(np.mean(delta)),
                "exact_one_sided_sign_p_for_decrease": float(
                    stats.binomtest(decreases, n, 0.5, alternative="greater").pvalue
                ),
                "exact_two_sided_sign_p": float(
                    stats.binomtest(decreases, n, 0.5).pvalue
                ),
                "interpretation": "Directional only; n=3 and exact inference is non-confirmatory.",
            }
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-source-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    patient_scores = (
        args.figure_source_root
        / "Supplementary_Figure_S7"
        / "ECM_LIGAND_21_patient_scores.csv"
    )
    effects = Path(__file__).parent / "resources" / "external_recurrence_effects_primary.csv"
    paired = paired_exact_sensitivities(patient_scores)
    meta, loo = meta_sensitivities(effects)
    spatial = spatial_shift_sensitivities(args.figure_source_root)
    spatial_summary = spatial_specimen_summary(args.figure_source_root)
    tgfb = tgfb_sign_sensitivity(args.figure_source_root)
    paired.to_csv(args.output_dir / "paired_exact_sensitivity.csv", index=False)
    meta.to_csv(args.output_dir / "meta_modified_hartung_knapp.csv", index=False)
    loo.to_csv(args.output_dir / "meta_leave_one_out_modified_hartung_knapp.csv", index=False)
    spatial.to_csv(args.output_dir / "spatial_toroidal_shift_sensitivity.csv", index=False)
    spatial_summary.to_csv(
        args.output_dir / "spatial_specimen_level_summary.csv", index=False
    )
    tgfb.to_csv(args.output_dir / "tgfb_exact_sign_sensitivity.csv", index=False)
    environment = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }
    (args.output_dir / "analysis_environment.json").write_text(
        json.dumps(environment, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
