import os
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

# List of sample names
sample_path = "/mnt/cbis/home/e1124735/Capstone/src/snakemake_Master/subworkflows/snakemake_dataprocess/config/samples.tsv"
samples = pd.read_csv(sample_path,sep="\t", dtype=str, comment="#",usecols=['sample']).T.values.tolist()[0]

# Base path template
base_path = "/mnt/gtklab01/xiaoqing/salmon/{}/full/aux_info/unmapped_names.txt"

pie_data_full = {}
pie_data_no_u = {}

# Loop through each sample
for sample in samples:
    file_path = base_path.format(sample)
    try:
        with open(file_path) as f:
            prefixes = [line.strip().split()[1] for line in f if len(line.strip().split()) > 1]
        total = len(prefixes)
        counts = Counter(prefixes)

        # Full with u
        full = [counts.get(k, 0) / total for k in ['d', 'm1', 'm2', 'u']] if total else [0, 0, 0, 0]
        pie_data_full[sample] = full

        # Normalized without u
        no_u_total = sum(counts.get(k, 0) for k in ['d', 'm1', 'm2'])
        if no_u_total:
            no_u = [counts.get(k, 0) / no_u_total for k in ['d', 'm1', 'm2']]
        else:
            no_u = [0, 0, 0]
        pie_data_no_u[sample] = no_u

    except FileNotFoundError:
        print(f"Warning: File not found for {sample}")
        pie_data_full[sample] = [0, 0, 0, 0]
        pie_data_no_u[sample] = [0, 0, 0]
        
print("\n📊 Start drawing graph")   
     
# Plotting all 8 pies in one figure
fig1, axes1 = plt.subplots(2, 4, figsize=(16, 8))
axes1 = axes1.flatten()
labels_full = ['d', 'm1', 'm2', 'u']
colors_full = ['#ff7f50', '#404969', '#bde4f4', 'grey']

for ax, sample in zip(axes1, samples):
    data = pie_data_full[sample]
    if sum(data) == 0:
        ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
        ax.axis('off')
    else:
        ax.pie(data, labels=labels_full, colors=colors_full, autopct='%1.1f%%', startangle=90)
        ax.set_title(sample)

fig1.suptitle("Unmapped Read Composition (Including 'u') [Decap decoys]", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("/mnt/gtklab01/xiaoqing/result_for_thesis/unmapped_decap_piecharts_with_u.png", dpi=300)
plt.savefig("/mnt/gtklab01/xiaoqing/result_for_thesis/unmapped_decap_piecharts_with_u.pdf")
print("📊 Finished first one.")   
# --- Figure 2: Without u ---
fig2, axes2 = plt.subplots(2, 4, figsize=(16, 8))
axes2 = axes2.flatten()
labels_no_u = ['d', 'm1', 'm2']
colors_no_u = ['#ff7f50', '#404969', '#bde4f4']

for ax, sample in zip(axes2, samples):
    data = pie_data_no_u[sample]
    if sum(data) == 0:
        ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
        ax.axis('off')
    else:
        ax.pie(data, labels=labels_no_u, colors=colors_no_u, autopct='%1.1f%%', startangle=90)
        ax.set_title(sample)

fig2.suptitle("Unmapped Read Prefix Composition (Excluding 'u') [Decap decoys]", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("/mnt/gtklab01/xiaoqing/result_for_thesis/unmapped_decap_piecharts_without_u.png", dpi=300)
plt.savefig("/mnt/gtklab01/xiaoqing/result_for_thesis/unmapped_decap_piecharts_without_u.pdf")

print("📊 DONE.")