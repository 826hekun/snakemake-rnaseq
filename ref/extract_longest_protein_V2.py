import argparse
import time
from Bio import SeqIO

# 解析 GFF3 文件，获取基因与其对应的 mRNA ID 关系
def parse_gene_mrna(gff_file):
    """Parse the GFF3 file to map each gene to its corresponding mRNA IDs."""
    gene_to_mrna = {}
    with open(gff_file, "r") as f:
        for line in f:
            if "\tmRNA\t" in line:
                fields = line.strip().split("\t")
                attributes = fields[8]
                mrna_id = [i.split("=")[1] for i in attributes.split(";") if i.startswith("ID=")][0]
                parent_id = [i.split("=")[1] for i in attributes.split(";") if i.startswith("Parent=")][0]
                if parent_id not in gene_to_mrna:
                    gene_to_mrna[parent_id] = []
                gene_to_mrna[parent_id].append(mrna_id)
    return gene_to_mrna

# 从 fasta 文件中提取蛋白序列
def get_protein_sequences(fasta_file):
    """Load the protein sequences from the provided fasta file."""
    return {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="Extract the longest protein sequence for each gene.")
    parser.add_argument("protein_fasta", help="Path to the protein fasta file.")
    parser.add_argument("gff_file", help="Path to the GFF3 file.")
    parser.add_argument("output_longest_protein", help="Output fasta file for the longest proteins.")
    parser.add_argument("output_protein_lengths", help="Output file detailing the lengths of proteins for each gene.")
    args = parser.parse_args()

    start_time = time.time()
    print("[INFO] Parsing the GFF3 file...")
    gene_to_mrna = parse_gene_mrna(args.gff_file)
    
    print("[INFO] Loading protein sequences...")
    protein_sequences = get_protein_sequences(args.protein_fasta)

    longest_proteins = []
    protein_lengths_detail = {}

    print("[INFO] Extracting the longest protein for each gene...")
    for gene, mrnas in gene_to_mrna.items():
        protein_lengths_detail[gene] = []
        longest_protein = None
        max_len = 0
        for mrna in mrnas:
            if mrna in protein_sequences:
                seq_len = len(protein_sequences[mrna].seq)
                protein_lengths_detail[gene].append((mrna, seq_len))
                if seq_len > max_len:
                    max_len = seq_len
                    longest_protein = protein_sequences[mrna]
        if longest_protein:
            longest_protein.id = gene
            longest_proteins.append(longest_protein)
    
    print("[INFO] Writing the longest proteins to the output file...")
    SeqIO.write(longest_proteins, args.output_longest_protein, "fasta")

    print("[INFO] Writing protein lengths details to the output file...")
    with open(args.output_protein_lengths, 'w') as out:
        for gene, lengths in protein_lengths_detail.items():
            for mrna, length in lengths:
                out.write(f"{gene}\t{mrna}\t{length}\n")

    end_time = time.time()
    print(f"[INFO] Done! Total time: {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()

