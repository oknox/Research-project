import numpy as np
import pandas as pd

import os

os.getcwd()
os.chdir("TPM_norm/")

counts = pd.read_csv("20180501_agg_gene_level_cnts.tsv", sep='\t')
# index_col=0 included in read.csv sets gene names to index
index = counts.index
columns = counts.columns
values = counts.values

# Selecting single column
counts["SRX193500"]

# Selecting multiple columns by passing it a list
counts[["SRX193500", "SRX146431", "SRX1329269"]]

# Use .loc to select for rows = df.loc[['Niko', 'Penelope']]
# using slice notation = df.loc['Niko':'Dean']
# From beginning to aaron: = df.loc[:'Aaron']
# Niko to Christina stepping by 2 = df.loc['Niko':'Christina':2]
# Dean to the end = df.loc['Dean':]

# Selecting 2 rows and 3 columns
counts.loc[["FBgn0000003", "FBgn0000008"], ["SRX193500", "SRX146431",
           "SRX1329269"]]
# Select all rows and columns - just use : - not necessary for selecting
# all columns

# Select single row with .iloc
counts.iloc[2]
# Selecting multiple rows
counts.iloc[[5, 3, 4]]
# Selecting rows and columns simulaneously
counts.iloc[[2, 3], [0, 4]]

# Set index from a column after reading in data - makes rownames
counts_idx = counts.set_index("FBgn")
index = counts_idx.index

genes = counts["FBgn"]

# Need the gene names to be matching = go back to R or do in python
# Can make genes index in gene_lengths by using set_index like above
# Google how to make pd.Series
gene_lengths = pd.read_csv("mean_gene_lengths.csv")
gene_lengths.columns = ["Gene", "mean_lengths"]

# Find the genes that are present in counts and gene_lengths
merged_inner = pd.merge(left=counts, right=gene_lengths, left_on='FBgn',
                        right_on='Gene')

# Trying a different way - above way gives the genes in the same order
gene_lengths_sim = gene_lengths[gene_lengths.Gene.isin(genes)]

# Make a pd.Series where index matches df.index and values are gene lengths
gene_series = pd.Series(merged_inner['mean_lengths'].values,
                        index=merged_inner['Gene'])
print(gene_series)

# Modify merged_inner so genes are rownames and the extra columns are removed
counts_mod = merged_inner.drop(['Gene', 'mean_lengths'], 1)
counts_mod = counts_mod.set_index("FBgn")

def tpm(df, gene_length, scale_library=1e6, scale_length=1e3, log=None):
    """Transcripts Per Killobase Million.
    Calcualtes TPM which normalizes by library size and gene length.
    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with genes as rows and samples as columns.
    gene_length : pd.Series
        Series where index matches df.index and values are gene lengths.
    scale_library : int or float
        Scaling factor to scale the library size.
    scale_length : int or float
        Scaling factor to scale gene model lengths.
    log : None or function or str
        If a function is giving will apply this to the data before returning.
        If a string is given then it uses numpy log functions that correspond
        to the provided string.
    """
    rpk = (df.T / (gene_length / scale_length)).T
    totals = rpk.sum()

    if log is None:
        log = lambda x: x - 1
    elif log == 'log2':
        log = np.log2
    elif log == 'log10':
        log = np.log10
    elif log == 'ln':
        log = np.log

    return log((rpk / (totals / scale_library)) + 1)


res = tpm(counts_mod, gene_series, log='log10')
res.to_csv("Log10_tpm_normalised_counts.csv")
res_nolog = tpm(counts_mod, gene_series)

# sum the columns to see if all add up to the same amount
res.sum(axis=0)
res_nolog.sum(axis=0)
res_nolog.to_csv("Nonlogged_tpm_normalised_counts.csv")

#Find maximum of the columns 
colmax = res_nolog.max()
colmax.max()

#Find minimum of the columns 
colmin = res_nolog.min()
colmin.min()