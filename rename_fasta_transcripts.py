import sys

def load_id_mappings(mapping_file):
    """Load old to new transcript ID mappings from a file.
       This will be produced by script create_transc2gene.pl  
    """

    id_map = {}
    with open(mapping_file, 'r') as mf:
        for line in mf:
            old_id, new_id = line.strip().split('\t')
            id_map[old_id] = new_id
    return id_map

def rename_fasta_headers(fasta_file, id_map, output_file):
    """Rename FASTA headers using the provided ID map."""
    with open(fasta_file, 'r') as ff, open(output_file, 'w') as of:
        for line in ff:
            if line.startswith('>'):
                # Extract old ID (before first space)
                header = line.split()[0]
                old_id = header[1:]  # Remove '>' from the header

                # Replace old ID with the new one if it exists in the map
                if old_id in id_map:
                    new_header = f'>{id_map[old_id]}'
                    # Preserve the rest of the header (after the first space)
                    rest_of_header = line[len(header):]
                    of.write(f"{new_header}{rest_of_header}")
                else:
                    # Write the original header if no mapping is found
                    of.write(line)
            else:
                # Write the sequence lines unchanged
                of.write(line)

if __name__ == "__main__":
    # Input arguments: mapping file, fasta file, output file
    if len(sys.argv) != 4:
        print("Usage: python rename_fasta_transcripts.py <mapping_file> <fasta_file> <output_file>")
        sys.exit(1)

    mapping_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    # Load the transcript ID mappings
    id_mappings = load_id_mappings(mapping_file)

    # Rename FASTA headers
    rename_fasta_headers(fasta_file, id_mappings, output_file)

    print(f"FASTA file '{fasta_file}' has been processed. Renamed headers saved to '{output_file}'.")
