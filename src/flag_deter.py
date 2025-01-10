import csv
import argparse
import os

def extract_reads_to_tsv(sam_file, id_file, output_tsv):
    # Load the IDs to a set for fast lookup, trimming any extra whitespace
    with open(id_file, 'r') as f:
        read_ids = set(line.strip() for line in f)

    with open(sam_file, 'r') as sam, open(output_tsv, 'w') as out_tsv:
        # Write header to the TSV file
        out_tsv.write("Read ID\tFLAG\n")
        
        for line in sam:
            if line.startswith('@'):  # Skip SAM header lines
                continue
            columns = line.strip().split('\t')
            read_id = columns[0]  # Read ID is in the first column
            flag = columns[1]     # FLAG is in the second column
            
            if read_id in read_ids:
                # Write the matching entry to the TSV file
                out_tsv.write(f"{read_id}\t{flag}\n")

    print(f"Extraction complete. Results saved in {output_tsv}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract reads from a SAM file based on an ID list and save to a TSV.")
    parser.add_argument("sam_file", type=str, help="Path to the SAM file.")
    parser.add_argument("id_file", type=str, help="Path to the file containing read IDs (one per line).")

    # Parse arguments
    args = parser.parse_args()

    # Automatically determine the output TSV file name
    id_file_name = os.path.splitext(os.path.basename(args.id_file))[0]  # Get the base name without extension
    output_tsv = f"extracted_{id_file_name}_reads.tsv"

    # Call the function with arguments
    extract_reads_to_tsv(args.sam_file, args.id_file, output_tsv)
