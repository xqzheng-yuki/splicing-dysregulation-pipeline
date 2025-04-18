import os
import pandas as pd
from collections import Counter

sample_path = "/mnt/cbis/home/e1124735/Capstone/src/snakemake_Master/subworkflows/snakemake_dataprocess/config/samples.tsv"
samples = pd.read_csv(sample_path,sep="\t", dtype=str, comment="#",usecols=['sample']).T.values.tolist()[0]

# Templates
full_d_template = "/mnt/gtklab01/xiaoqing/salmon/{}/full/aux_info/unmapped_d.lst"
mash_d_template = "/mnt/gtklab01/xiaoqing/salmon/{}/mashmap/aux_info/unmapped_d.lst"
unmapped_names_template = "/mnt/gtklab01/xiaoqing/salmon/{}/full/aux_info/unmapped_names.txt"

# Store all results
all_counts = {}

for sample in samples:
    try:
        # Load both unmapped_d.lst sets
        with open(full_d_template.format(sample)) as f:
            full_ids = set(line.strip() for line in f if line.strip())
        with open(mash_d_template.format(sample)) as f:
            mash_ids = set(line.strip() for line in f if line.strip())
        
        # Find unique IDs
        mash_only_ids = mash_ids - full_ids

        # Count types from full/unmapped_names.txt
        matched_types = []
        with open(unmapped_names_template.format(sample)) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0] in mash_only_ids:
                    matched_types.append(parts[1])  # column 2 is the type

        count = Counter(matched_types)
        all_counts[sample] = count

        # Print result per sample
        print(f"== {sample} ==")
        for t in ['d', 'm1', 'm2', 'u']:
            print(f"{t}: {count.get(t, 0)}")
        print()

    except FileNotFoundError as e:
        print(f"[Missing file for {sample}]: {e}")
        all_counts[sample] = {}

# Optionally write to TSV file
with open("/mnt/gtklab01/xiaoqing/result_for_thesis/shared_type_counts.tsv", "w") as out:
    out.write("Sample\td\tm1\tm2\tu\n")
    for sample in samples:
        c = all_counts[sample]
        out.write(f"{sample}\t{c.get('d',0)}\t{c.get('m1',0)}\t{c.get('m2',0)}\t{c.get('u',0)}\n")

print("\n✔ Done! Summary written\n")