#!/usr/bin/env python
import os
import glob
import numpy as np
import gseapy as gp
import pandas as pd
import argparse

args = argparse.ArgumentParser()
args.add_argument('--in_dir', help='Input directory', default='Data')
args.add_argument('--out_dir', help='Output directory', default='Output')
args.add_argument('--sig_dir', help='Signature directory', default='Signatures')
args.add_argument('-g', '--gse', help='GSE file', default=None)
args.add_argument('-m', '--gmt', help='GMT file', default=None)
args = args.parse_args()

# Read the GSE ID
def run_ssgsea(GSE, gmt):
    # Read the counts
    if GSE.endswith('.tsv'):
        df = pd.read_csv(f"./{args.in_dir}/{GSE}",sep='\t', index_col=0)
    elif GSE.endswith('.csv'):
        df = pd.read_csv(f"./{args.in_dir}/{GSE}", index_col=0)
    else:
        print(f'Unknown file format for {GSE}')
        return
    # # Remove space in the genenames and make it caps 
    df.index = df.index.str.replace(' ','')
    df.index = df.index.str.upper()
    # Get the gene sets
    geneSets = f'./{args.sig_dir}/{gmt}'
    # Calculate ssGSEA scores
    try:
        ss = gp.ssgsea(data=df,
                gene_sets=geneSets,
                outdir=None,
                sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
                no_plot=True, threads=os.cpu_count()-2)
    except Exception as error:
        print(f'Error occured for {GSE} and {gmt}')
        print(error)
        return
    else:
        # Convert to dataframe with each signature as a column
        df = ss.res2d.pivot(index='Name', columns='Term', values='NES')
        df = df.convert_dtypes()
        # Save the output
        GSE = GSE.split('.')[0].split('_')[0]
        gmt = gmt.split('.')[0]
        df.to_csv(f"./{args.out_dir}/{GSE}_{gmt}-ssGSEA.csv", index=True)
        
# List files in Output
if args.gse is not None:
    GSEs = [args.gse]
else:
    # Get the GSEs
    GSEs = glob.glob('*.tsv', root_dir=args.in_dir) + glob.glob('*.csv', root_dir=args.in_dir)
# Get the GMTs
if args.gmt is not None:
    GMTs = [args.gmt]
else:
    GMTs = glob.glob('*.gmt', root_dir=args.sig_dir)
# Make directory for figures and output
os.makedirs(args.out_dir, exist_ok=True)

for GSE in GSEs:
    for gmt in GMTs:
        print(f'Running ssGSEA for {GSE} and {gmt}')
        run_ssgsea(GSE, gmt)
