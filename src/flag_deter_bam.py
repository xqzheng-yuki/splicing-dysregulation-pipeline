import pysam
import argparse
import os

def extract_reads_to_tsv(bam_file, id_file, output_tsv):
    # Load the IDs to a set for fast lookup, trimming any extra whitespace
    with open(id_file, 'r') as f:
        read_ids = set(line.strip() for line in f)

    # Open BAM file using pysam
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_tsv, 'w') as out_tsv:
        # Write header to the TSV file
        out_tsv.write("Read ID\tFLAG\n")
        
        for read in bam:
            if read.query_name in read_ids:  # Check if the read ID matches
                # Write the matching entry to the TSV file
                out_tsv.write(f"{read.query_name}\t{read.flag}\n")

    print(f"Extraction complete. Results saved in {output_tsv}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract reads from a BAM file based on an ID list and save to a TSV.")
    parser.add_argument("bam_file", type=str, help="Path to the BAM file.")
    parser.add_argument("id_file", type=str, help="Path to the file containing read IDs (one per line).")

    # Parse arguments
    args = parser.parse_args()

    # Automatically determine the output TSV file name
    id_file_name = os.path.splitext(os.path.basename(args.id_file))[0]  # Get the base name without extension
    output_tsv = f"extracted_{id_file_name}_reads_ori.tsv"

    # Call the function with arguments
    extract_reads_to_tsv(args.bam_file, args.id_file, output_tsv)
