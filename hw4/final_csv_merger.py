#!/home/randhal/proj/genomic_env/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Merge the columns from the original dataframe with the 
# columns from snpnexus and fathmm
def add_columns(df_original, df_snpnexus, df_fathmm, outname):
    # For this case, after careful analysis dropi_duplicates is not required. 

    df_snpnexus = df_snpnexus.rename(columns={
        'id': 'id_nexus',
        'chrom': 'chrom_nexus',
        'position': 'position_nexus', 
        'ref': 'ref_nexus',
        'alt': 'alt_nexus',
        'Genomic Coordinates': 'Gen.Coor Nexus',
        'Overlapped or Nearest Genes': 'OvrlppdNrst.Genes Nexus',
        'SIFT': 'SIFT_NEXUS',
        'CADD - Combined Annotation Dependent Depletion':'CADD_NEXUS', 
        'Ensembl': 'Ensemble_nexus',
        'RefSeq':'RefSeq Nexus'
        })

    # Create an UUID to merge the files
    df_original['UUID'] = df_original['left'].astype(str)
    df_snpnexus['UUID'] = df_snpnexus['position_nexus'].astype(str)
    df_fathmm['UUID'] = df_fathmm['Position'].astype(str)

    merged_df = pd.merge(df_original, df_snpnexus, on='UUID', how='inner')
    final_df = pd.merge(merged_df, df_fathmm, on='UUID', how='inner')

    final_df = final_df.drop(columns = ['UUID', 'Unnamed: 0'])

    final_df.to_csv(outname, index=False)


if __name__=='__main__':
    # Initial dataframes. 
    data_normal = pd.read_csv("ALL_Normal.csv")
    data_tumor = pd.read_csv("ALL_Tumor.csv")

    # Data from snpnexus
    nexus_normal = pd.read_csv("normal_study/normal_study_pervariant.tsv", sep='\t')
    nexus_tumor  = pd.read_csv("tumor_study/tumor_study_pervariant.tsv", sep='\t')

    # Data from fathmm
    fathmm_normal = pd.read_csv("normal_FATHMM.txt", sep='\t')
    fathmm_tumor = pd.read_csv("tumor_FATHMM.txt", sep='\t')

    # merge all 
    add_columns(data_normal, nexus_normal, fathmm_normal, "ramirez_hw4_normal.csv") 
    add_columns(data_tumor, nexus_tumor, fathmm_tumor, "ramirez_hw4_tumor.csv") 
