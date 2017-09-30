import argparse
import feather
import magic
import os

parser = argparse.ArgumentParser(description='wrapper for magic')

parser.add_argument('--matx', dest='matx',help='Matx path')
parser.add_argument('--out', dest='out',help='Output path')
parser.add_argument('--maxCellSize', dest='maxCS',type=int,default=1000000,help='Max num of reads allow in a cell')
parser.add_argument('--minCellSize', dest='minCS',type=int,default=1,help='Min num of reads allow in a cell')

args = parser.parse_args()

# Load single-cell RNA-seq data
df = feather.read_dataframe(args.matx)
scdata = magic.mg.SCData(df, 'sc-seq')
#scdata = magic.mg.SCData.from_csv(os.path.expanduser(args.matx),data_type='sc-seq', normalize=False)

# Minimum molecules/cell value
CELL_MIN = args.minCS

# Maximum molecules/cell values
CELL_MAX = args.maxCS

# Minimum number of nonzero cells/gene 
# (None if no filtering desired)
GENE_NONZERO = None

# Minimum number of molecules/gene
# (None if no filtering desired)
GENE_MOLECULES = None

scdata.filter_scseq_data(filter_cell_min=CELL_MIN, filter_cell_max=CELL_MAX, 
                         filter_gene_nonzero=GENE_NONZERO, filter_gene_mols=GENE_MOLECULES)

scdata = scdata.normalize_scseq_data()

# MAGIC
scdata.run_magic(n_pca_components=20, random_pca=True, t=6, k=30,ka=10, epsilon=1, rescale_percent=99)

feather.write_dataframe(scdata.magic.data, args.out)


