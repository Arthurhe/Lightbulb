import argparse
import feather
import magic
#import os

parser = argparse.ArgumentParser(description='wrapper for magic')

parser.add_argument('--matx', dest='matx',help='Matx path')
parser.add_argument('--out', dest='out',help='Output path')
#parser.add_argument('--maxCellSize', dest='maxCS',type=int,default=1000000,help='Max num of reads allow in a cell')
#parser.add_argument('--minCellSize', dest='minCS',type=int,default=1,help='Min num of reads allow in a cell')

args = parser.parse_args()

# Load single-cell RNA-seq data
df = feather.read_dataframe(args.matx)
scdata = magic.mg.SCData(df, 'sc-seq')

# MAGIC
scdata.run_magic(n_pca_components=20, random_pca=True, t=6, k=30,ka=10, epsilon=1, rescale_percent=99)

#output
feather.write_dataframe(scdata.magic.data, args.out)


