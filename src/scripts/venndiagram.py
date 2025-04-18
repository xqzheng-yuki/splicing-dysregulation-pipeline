import os
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2

# Sample names
sample_path = "/mnt/cbis/home/e1124735/Capstone/src/snakemake_Master/subworkflows/snakemake_dataprocess/config/samples.tsv"
samples = pd.read_csv(sample_path,sep="\t", dtype=str, comment="#",usecols=['sample']).T.values.tolist()[0]

# File path templates
full_template = "/mnt/gtklab01/xiaoqing/salmon/{}/full/aux_info/unmapped_d.lst"
mash_template = "/mnt/gtklab01/xiaoqing/salmon/{}/mashmap/aux_info/unmapped_d.lst"

# Initialize plot
fig, axes = plt.subplots(2, 4, figsize=(20, 10), gridspec_kw={'hspace': 0.3})
axes = axes.flatten()

for i, sample in enumerate(samples):
    file_full = full_template.format(sample)
    file_mash = mash_template.format(sample)

    try:
        with open(file_full) as f:
            full_set = set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        print(f"File not found: {file_full}")
        full_set = set()

    try:
        with open(file_mash) as f:
            mash_set = set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        print(f"File not found: {file_mash}")
        mash_set = set()

    ax = axes[i]
    if not full_set and not mash_set:
        ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
        ax.axis('off')
    else:
        venn = venn2([full_set, mash_set], set_labels=('full', 'mashmap'), ax=ax)

        # Percentages
        shared = len(full_set & mash_set)
        mash_total = len(mash_set)
        full_total = len(full_set)

        if mash_total > 0:
            only_full_pct = 100 * (full_total - shared) / full_total
            shared_pct = 100 * shared / mash_total
        else:
            only_full_pct = 0
            shared_pct = 0

        title_text = f"{sample}\n{only_full_pct:.1f}% only shown in Decap\n{shared_pct:.1f}% shared based on mashmap number"
        ax.set_title(title_text, fontsize=11)

plt.tight_layout()
plt.savefig("/mnt/gtklab01/xiaoqing/result_for_thesis/venn/venn_unmapped_d.png", dpi=300)
plt.savefig("/mnt/gtklab01/xiaoqing/result_for_thesis/venn/venn_unmapped_d.pdf")