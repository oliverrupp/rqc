import sys
import re

# Function to process the input file and generate the output
def process_file(input_file, output_file):
    c = 1
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith("#"):
                continue

            x = line.split("\t")

            # Extract the Name value from the 9th field using regex
            match = re.search(r"Name=([^;]*)", x[8])
            if match:
                t = match.group(1)

                id_ = f"rRNA_{c}_{t}"

                # Modify the values in the list
                x[2] = "exon"
                x[8] = f'transcript_id "{id_}"; gene_id "{id_}"'

                # Write the modified line to the output file
                outfile.write("\t".join(x) + "\n")

                c += 1

# Main execution with file arguments
if __name__ == "__main__":
    input_fasta = snakemake.input[0]
    output_fasta = snakemake.output[0]
    process_file(input_fasta, output_fasta)
