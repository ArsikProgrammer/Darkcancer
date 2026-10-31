"""Generate figures for the preprint from real project data."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def fig1_workflow():
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.axis("off")
    boxes = [
        (0.02, 0.62, 0.18, 0.28, "Cohorts\nSCLC n=143\nLUAD n=230"),
        (0.24, 0.62, 0.16, 0.28, "Residue validation\n(transcript + UniProt)\n11-21% discarded"),
        (0.44, 0.62, 0.16, 0.28, "Disorder label\nAlphaFold pLDDT<50\n~95% coverage"),
        (0.64, 0.62, 0.17, 0.28, "Expected counts E\n96-context spectrum\nx coding-space bg"),
        (0.85, 0.62, 0.13, 0.28, "Model\nO ~ log(E)\n+ SBS4 + TMB"),
        (0.30, 0.12, 0.42, 0.30,
         "Controls: SBS1 predictor, silent outcome\nNull test: background resampling x10k\nEquivalence bounds (0.833, 1.2) fixed in advance"),
    ]
    for x, y, w, h, text in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, fill=True, fc="#eef2f7", ec="#334155", lw=1))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.5)
    for x1, x2, y in [(0.20, 0.24, 0.76), (0.40, 0.44, 0.76), (0.60, 0.64, 0.76), (0.81, 0.85, 0.76)]:
        ax.annotate("", xy=(x2, y), xytext=(x1, y),
                    arrowprops=dict(arrowstyle="->", color="#334155", lw=1))
    ax.annotate("", xy=(0.51, 0.42), xytext=(0.51, 0.62),
                arrowprops=dict(arrowstyle="->", color="#334155", lw=1))
    fig.tight_layout()
    fig.savefig(OUT / "fig1_workflow.png", dpi=300)
    plt.close(fig)


def fig2_smoking_scale():
    sig = pd.read_csv(PROC / "luad_signature_activities.csv", low_memory=False)
    sig["Tumor_Sample_Barcode"] = sig["Tumor_Sample_Barcode"].astype(str)
    sig["patient"] = sig["Tumor_Sample_Barcode"].str[:12]
    sm = pd.read_csv(RAW / "luad_smoking_data_clinical_patient.txt", sep="\t",
                     comment="#", dtype=str)
    sm = sm[["PATIENT_ID", "TOBACCO_SMOKING_HISTORY_INDICATOR"]].copy()
    sm["status"] = pd.to_numeric(sm["TOBACCO_SMOKING_HISTORY_INDICATOR"].replace(
        "[Not Available]", None), errors="coerce")
    df = sig.merge(sm, left_on="patient", right_on="PATIENT_ID", how="inner")
    order = [1, 3, 4, 2]
    labels = ["Never", "Former\n>15 y", "Former\n<=15 y", "Current"]
    data = [df.loc[df["status"] == s, "SBS4_fraction"].dropna().values for s in order]
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.boxplot(data, tick_labels=labels, widths=0.55, showfliers=False)
    ax.set_ylabel("SBS4 fraction (NNLS)")
    ax.set_xlabel("Smoking status (TCGA LUAD)")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_smoking_scale.png", dpi=300)
    plt.close(fig)


def fig3_forest():
    rows = [
        ("SCLC primary (per 1 SD SBS4)", 1.0084, 0.9849, 1.0326),
        ("LUAD primary (per 1 SD SBS4)", 0.9998, 0.9752, 1.0250),
        ("LUAD pack-years (per 1 SD)", 0.9991, 0.9775, 1.0213),
        ("LUAD never vs current (per 1 SD SBS4)", 1.0167, 0.9657, 1.0705),
        ("LUAD discard covariate", 0.9991, 0.9744, 1.0243),
        ("Control: SBS1 predictor (SCLC)", 0.9938, 0.9725, 1.0155),
        ("Control: silent outcome (SCLC)", 1.0317, 0.9825, 1.0835),
        ("Control: SBS1 predictor (LUAD)", 0.9789, 0.9423, 1.0170),
        ("Control: silent outcome (LUAD)", 1.0123, 0.9626, 1.0646),
    ]
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    ax.axvspan(0.833, 1.2, color="#e8f0e8", zorder=0, label="equivalence bounds")
    ax.axvline(1.0, color="#555", lw=0.8, ls="--")
    y = np.arange(len(rows))[::-1]
    for yi, (name, est, lo, hi) in zip(y, rows):
        ax.plot([lo, hi], [yi, yi], color="#1f4e79", lw=1.6)
        ax.plot(est, yi, "o", color="#1f4e79", ms=4.5)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.set_xlabel("O/E ratio per 1 SD (95% CI)")
    ax.set_xscale("log")
    ax.set_xlim(0.9, 1.2)
    ax.set_xticks([0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2])
    ax.set_xticklabels(["0.90", "0.95", "1.00", "1.05", "1.10", "1.15", "1.20"])
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_forest.png", dpi=300)
    plt.close(fig)


def fig4_raw_scatter():
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9))
    sets = [
        ("SCLC", ["annotated_george.csv", "annotated_jhu.csv", "annotated_clcgp.csv"],
         "tumor_signature_activities_refit_dedup.csv", "r = -0.12 (Pearson)"),
        ("LUAD", ["annotated_luad.csv"], "luad_signature_activities.csv", "r = +0.08 (Pearson)"),
    ]
    for ax, (name, mut_files, sig_file, title) in zip(axes, sets):
        frames = []
        for f in mut_files:
            d = pd.read_csv(PROC / f, low_memory=False)
            d["Tumor_Sample_Barcode"] = d["Tumor_Sample_Barcode"].astype(str)
            frames.append(d)
        mut = pd.concat(frames, ignore_index=True)
        mut = mut[mut["is_dark"].notna()]
        per = mut.groupby("Tumor_Sample_Barcode")["is_dark"].mean().rename("dark").reset_index()
        sig = pd.read_csv(PROC / sig_file, low_memory=False)
        sig["Tumor_Sample_Barcode"] = sig["Tumor_Sample_Barcode"].astype(str)
        df = per.merge(sig, on="Tumor_Sample_Barcode", how="inner")
        ax.scatter(df["SBS4_fraction"], df["dark"], s=10, alpha=0.55, color="#1f4e79")
        coef = np.polyfit(df["SBS4_fraction"], df["dark"], 1)
        xs = np.linspace(df["SBS4_fraction"].min(), df["SBS4_fraction"].max(), 10)
        ax.plot(xs, coef[0] * xs + coef[1], color="#b03a2e", lw=1)
        ax.set_title(f"{name}  ({title})", fontsize=8)
        ax.set_xlabel("SBS4 fraction")
        ax.set_ylabel("dark mutation share")
    fig.tight_layout()
    fig.savefig(OUT / "fig4_raw_scatter.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    fig1_workflow()
    fig2_smoking_scale()
    fig3_forest()
    fig4_raw_scatter()
    print("figures saved to", OUT)
