from Bio import SeqIO

# 1. 从 GFF3 文件中解析每个基因与其对应的所有 mRNA ID。
def parse_gene_mrna(gff_file):
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

# 2. 从蛋白序列 fasta 文件中读取每个 mRNA ID 对应的蛋白序列。
def get_protein_sequences(fasta_file):
    return {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

# 主程序开始
gff_file = "raw/SWO.v3.0.gene.model.gff3"
fasta_file = "raw/SWO.v3.0.protein.fa"
output_file = "longest_protein.fa"

gene_to_mrna = parse_gene_mrna(gff_file)
protein_sequences = get_protein_sequences(fasta_file)

# 3. 对每个基因，从其对应的所有 mRNA 蛋白序列中选择最长的一个，并保存到输出文件中。
longest_proteins = []
for gene, mrnas in gene_to_mrna.items():
    longest_protein = None
    max_len = 0
    for mrna in mrnas:
        if mrna in protein_sequences and len(protein_sequences[mrna].seq) > max_len:
            max_len = len(protein_sequences[mrna].seq)
            longest_protein = protein_sequences[mrna]
    if longest_protein:
        longest_protein.id = gene
        longest_proteins.append(longest_protein)

SeqIO.write(longest_proteins, output_file, "fasta")

