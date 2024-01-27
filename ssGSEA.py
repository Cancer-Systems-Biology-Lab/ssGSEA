import os
import glob
import numpy as np
import gseapy as gp
import pandas as pd
from multiprocessing import cpu_count

# Read the GSE ID
def run_ssgsea(GSE, gmt):
    # Read the counts
    df = pd.read_csv("./Data/"+GSE+"_TPM.tsv",sep='\t', index_col=0)
    # # Remove space in the genenames and make it caps 
    df.index = df.index.str.replace(' ','')
    df.index = df.index.str.upper()
    # Get the gene sets
    geneSets = './Signatures/'+gmt+'.gmt'
    # Calculate ssGSEA scores
    try:
        ss = gp.ssgsea(data=df,
                gene_sets=geneSets,
                outdir=None,
                sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
                no_plot=True, threads=cpu_count()-2)
    except Exception as error:
        print(f'Error occured for {GSE} and {gmt}')
        print(error)
        return
    else:
        # Convert to dataframe with each signature as a column
        df = ss.res2d.pivot(index='Name', columns='Term', values='NES')
        df = df.convert_dtypes()
        # Save the output
        df.to_csv("./Output/"+GSE+"-"+gmt+"-ssGSEA.csv", index=True)
        
# List files in Output
# Get the GSEs
GSEs = glob.glob('Data/*_TPM.tsv')
GSEs = [os.path.basename(x).replace('_TPM.tsv','') for x in GSEs]
# Get the GMTs
GMTs = glob.glob('Signatures/*.gmt')
GMTs = [os.path.basename(x).replace(f'.gmt','') for x in GMTs]
# Make directory for figures and output
os.makedirs('../Output/', exist_ok=True)

for GSE in GSEs:
    for gmt in GMTs:
        print(f'Running ssGSEA for {GSE} and {gmt}')
        run_ssgsea(GSE, gmt)
