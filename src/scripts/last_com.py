import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import math
import seaborn as sns
from scipy.stats import pearsonr

sample_path = "/mnt/cbis/home/e1124735/Capstone/src/snakemake_Master/subworkflows/snakemake_dataprocess/config/samples.tsv"
# Your custom sample names
samples = pd.read_csv(sample_path,sep="\t", dtype=str, comment="#",usecols=['sample']).T.values.tolist()[0]

directory = "/mnt/gtklab01/xiaoqing/result_for_thesis/origin+decap"
result_file = os.path.join(directory,"salmon_comparison.csv")

summary_rows = []
plots_data = []

for sample in samples:
    full_dir = os.path.join("/mnt/gtklab01/xiaoqing/salmon", sample, "origin")
    # filt_dir = os.path.join("/mnt/gtklab01/xiaoqing/salmon", sample, "mashmap")
    decap_dir = os.path.join("/mnt/gtklab01/xiaoqing/salmon", sample, "full")

    try:
        with open(os.path.join(full_dir, "aux_info", "meta_info.json")) as f:
            meta_full = json.load(f)

        with open(os.path.join(decap_dir, "aux_info", "meta_info.json")) as f:
            meta_decap = json.load(f)

        # TPM tables
        df_full = pd.read_csv(os.path.join(full_dir, "quant.sf"), sep="\t")
        df_decap = pd.read_csv(os.path.join(decap_dir, "quant.sf"), sep="\t")

        df_merged = df_full[["Name", "TPM"]].merge(
            df_decap[["Name", "TPM"]],
            on="Name",
            suffixes=("_full", "_decap")
        )

        # Correlation
        corr, _ = pearsonr(df_merged["TPM_full"], df_merged["TPM_decap"])

        # Save plot
        plt.figure(figsize=(6,6))
        sns.scatterplot(data=df_merged, x="TPM_full", y="TPM_decap", alpha=0.5)
        plt.xscale("log")
        plt.yscale("log")
        plt.plot([1e-3, 1e5], [1e-3, 1e5], 'r--')
        plt.xlabel("TPM (Original Decoy)")
        plt.ylabel("TPM (Final Decoy)")
        plt.title(f"TPM Correlation: {sample} (r = {corr:.2f})")
        plt.tight_layout()
        plt.savefig(os.path.join(directory,f"TPM_correlation_{sample}.svg"))
        plt.close()
        plots_data.append((sample, df_merged.copy(), corr))

        # Summary table
        summary_rows.append({
            "Sample": sample,
            "num_processed_full": meta_full.get("num_processed", 0),
            "num_mapped_full": meta_full.get("num_mapped", 0),
            "num_decoy_fragments_full": meta_full.get("num_decoy_fragments", 0),
            "percent_mapped_full": meta_full.get("percent_mapped", 0),

            "num_processed_decap": meta_decap.get("num_processed", 0),
            "num_mapped_decap": meta_decap.get("num_mapped", 0),
            "num_decoy_fragments_decap": meta_decap.get("num_decoy_fragments", 0),
            "percent_mapped_decap": meta_decap.get("percent_mapped", 0),

            "TPM_correlation": round(corr, 7)
        })

        print(f"✅ Processed: {sample}")
    except Exception as e:
        print(f"❌ Skipped {sample} due to error: {e}")

# Save summary
df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(result_file, index=False)
print("\n📁 Summary table saved as salmon_comparison_summary_all_samples.csv")
print("📊 Correlation plots saved as TPM_correlation_<sample>.svg")

# Plot all subplots in one large figure
num_plots = len(plots_data)
cols = 4  # number of columns in the grid
rows = math.ceil(num_plots / cols)
fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))

# Flatten axes array in case it's 2D
axes = axes.flatten()

for i, (sample, df_merged, corr) in enumerate(plots_data):
    ax = axes[i]
    sns.scatterplot(data=df_merged, x="TPM_full", y="TPM_decap", alpha=0.5, ax=ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot([1e-3, 1e5], [1e-3, 1e5], 'r--')
    ax.set_title(f"{sample}\nr = {corr:.2f}")
    ax.set_xlabel("TPM (Original Decoy)")
    ax.set_ylabel("TPM (Final Decoy)")

# Turn off any unused axes
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

fig.suptitle("TPM Correlation Across Samples", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.97])

# Save combined figure
plt.savefig(os.path.join(directory,"TPM_correlation_combined.pdf"))
plt.savefig(os.path.join(directory,"TPM_correlation_combined.png"))
plt.close()
print(f"🖼️ Combined plot saved at {directory}")