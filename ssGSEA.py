import os
import glob
import numpy as np
import gseapy as gp
import pandas as pd

# Read the GSE ID
def run_ssgsea(GSE, gmt):
    # Read the counts
    if GSE.endswith('.tsv'):
        df = pd.read_csv("./Data/"+GSE,sep='\t', index_col=0)
    elif GSE.endswith('.csv'):
        df = pd.read_csv("./Data/"+GSE, index_col=0)
    else:
        print(f'Unknown file format for {GSE}')
        return
    # # Remove space in the genenames and make it caps 
    df.index = df.index.str.replace(' ','')
    df.index = df.index.str.upper()
    # Get the gene sets
    geneSets = './Signatures/'+gmt
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
        df.to_csv("./Output/"+GSE+"_"+gmt+"-ssGSEA.csv", index=True)
        
# List files in Output
# Get the GSEs
GSEs = glob.glob('*.tsv', root_dir='Data')
# Get the GMTs
GMTs = glob.glob('*.gmt', root_dir='Signatures')
# Make directory for figures and output
os.makedirs('../Output/', exist_ok=True)

for GSE in GSEs:
    for gmt in GMTs:
        print(f'Running ssGSEA for {GSE} and {gmt}')
        run_ssgsea(GSE, gmt)
