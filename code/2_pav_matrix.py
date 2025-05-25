#!/usr/bin/env python3

"""
Generate a binary gene presence matrix from annotated gene lists.

This script processes all `*.sorted_annotated_genes.txt` files in the current directory.
Each file corresponds to a pairwise genome comparison and contains a list of gene identifiers
considered present in the comparison. The output is a binary matrix indicating gene presence in 
genome A and absent in genome B(1) or shared gene (0) in one comparison.

Output:
    - gene_presence_matrix.tsv: a tab-separated file with genes as rows and comparisons as columns.

Example:
    $ python generate_gene_matrix.py

"""

import pandas as pd
import glob
from collections import defaultdict
import os


def load_genes(filename):

    with open(filename) as f:
        return set(line.strip() for line in f if line.strip())


annotation_files = glob.glob("*.sorted_annotated_genes.txt")

# Data structure to store gene presence
# Format: {gene_id: {comparison_name: 1}}
gene_presence = defaultdict(dict)

for filepath in annotation_files:
    # Extract comparison name from filename
    comparison_name = os.path.basename(filepath).replace(
        ".sorted_annotated_genes.txt", ""
    )

    # Load genes present in the comparison
    genes = load_genes(filepath)

    # Record gene presence
    for gene in genes:
        gene_presence[gene][comparison_name] = 1

all_comparisons = sorted(
    os.path.basename(f).replace(".sorted_annotated_genes.txt", "")
    for f in annotation_files
)

df = pd.DataFrame.from_dict(gene_presence, orient="index").fillna(0).astype(int)

df = df[all_comparisons]

output_file = "gene_presence_matrix.tsv"
df.to_csv(output_file, sep="\t")

print(f"🧬 Matrix saved to {output_file}")
print(f"🔢 Total genes: {df.shape[0]} | Comparisons: {df.shape[1]}")
print(df.head())
