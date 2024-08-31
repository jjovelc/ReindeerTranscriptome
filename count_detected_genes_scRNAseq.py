from scipy.io import mmread
import numpy as np
import sys

matrix = sys.argv[1]
features = sys.argv[2]
prefix = sys.argv[3]

# Read the matrix
mtx = mmread(matrix).tocsc()

# Determine which genes were detected (non-zero counts in any cell)
detected_genes = np.array(mtx.sum(axis=1) > 0).flatten()

# Get indices of detected genes
detected_gene_indices = np.where(detected_genes)[0]

# Map indices to gene names
# Read the features (genes) file
import pandas as pd
genes = pd.read_csv(features, header=None, sep="\t")

# Subset to detected genes
detected_gene_names = genes.iloc[detected_gene_indices]

# Save results
outfile = f'{prefix}_detected_genes.txt'
detected_gene_names.to_csv(outfile, index=False, header=False)
