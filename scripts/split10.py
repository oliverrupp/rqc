def read_fasta(filepath):
    sequences = {}
    with open(filepath, 'r') as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                # Keep only the first word before any whitespace
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    return sequences

def split_sequence(seq, parts=10):
    length = len(seq)
    if length < parts:
        raise ValueError(f"Sequence too short to split into {parts} parts.")
    part_length = length // parts
    if part_length < 150:
        raise ValueError(f"Sequence too short to split into {parts} parts.")
    return [seq[i * part_length:(i + 1) * part_length] for i in range(parts)]

def write_fasta(sequences, output_path):
    with open(output_path, 'w') as f:
        for header, subseqs in sequences.items():
            for i, subseq in enumerate(subseqs):
                new_header = f">{header}_q{i+1}"
                f.write(new_header + '\n')
                for j in range(0, len(subseq), 60):
                    f.write(subseq[j:j+60] + '\n')

def main(input_file, output_file):
    fasta_seqs = read_fasta(input_file)
    split_seqs = {}
    for header, seq in fasta_seqs.items():
        try:
            parts = split_sequence(seq, 10)
            split_seqs[header] = parts
        except ValueError as e:
            print(f"Skipping {header}: {e}")
    write_fasta(split_seqs, output_file)

# Example usage
if __name__ == '__main__':
    input_fasta = snakemake.input[0]
    output_fasta = snakemake.output[0]
    main(input_fasta, output_fasta)
