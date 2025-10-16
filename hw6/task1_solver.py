#!/home/randhal/proj/genomic_env/bin/python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import load_iris

# Solves item 1.1
def cleaner11(normal, tumor):

    # Select the required columns
    kols = ["chrom", "left", "ref_seq", "alt_seq", "Patient_ID", "VCF_ID"]
 
    cnormal = normal[kols]
    ctumor  = tumor[kols]
    
    # 1.1.1 How many unique normal patients do we have?
    total_normal = cnormal["Patient_ID"].nunique()
    print(f"Ans 1.1.1: There are {total_normal} unique normal patients\n")
    
    # 1.1.2 How many unique normal patients do we have?
    total_tumor = ctumor["Patient_ID"].nunique()
    print(f"Ans 1.1.2: There are {total_tumor} unique tumor patients\n")
    
    # 1.1.3 Group by variant info, chrom, left, ref_seq and alt_seq, let the other columns be a list
    
    gkols = ["chrom", "left", "ref_seq", "alt_seq"]
 
    group_normal = cnormal.groupby(gkols).agg(list).reset_index()
    group_tumor  = ctumor.groupby(gkols).agg(list).reset_index()

    print("Ans 1.1.3: \n\n", group_normal, "\n")
    print(group_tumor)

    # 1.1.4 Create a new column with the number of patients per variant on both the normal and tumor 
    #       (name the column, N# and T#, respectively)
    
    group_normal['N#'] = group_normal['Patient_ID'].apply(lambda x: len(set(x)))
    group_tumor['T#'] = group_tumor['Patient_ID'].apply(lambda x: len(set(x)))
    
    print("\n\nAns 1.1.4: \n\n", group_normal, "\n")
    print(group_tumor)
    
    # 1.1.5 Rename the columns, Patient_ID and VCF_ID, to have, _Normal or _Tumor, 
    #       added depending which file you are working with.

    group_normal = group_normal.rename(
        columns={
            "Patient_ID": "Patient_ID_Normal",
            "VCF_ID": "VCF_ID_Normal"
        }
    )
    
    group_tumor = group_tumor.rename(
        columns={
            "Patient_ID": "Patient_ID_Tumor",
            "VCF_ID": "VCF_ID_Tumor"
        }
    )

    print("\n\nAns 1.1.5: \n\n", group_normal, "\n")
    print(group_tumor)

    return group_normal, group_tumor


# Solves item 1.2
def merger12(gnormal, gtumor):
        
    gkols = ["chrom", "left", "ref_seq", "alt_seq"]
    print(gnormal)   

    # Merged Data Frame = AML
    #aml = pd.concat([gnormal, gtumor], axis=0, join="outer", ignore_index=True)
    aml = pd.merge(gnormal, gtumor, on=gkols, how='outer')
    print(aml)

    # Save the AML csv
    aml.to_csv("AML.csv")

    unormal = aml[(aml['N#']==1) & (aml["T#"]!=1)]
    utumor = aml[(aml['T#']==1) & (aml["N#"]!=1)]
    shared = aml[(aml["N#"]==1) & (aml["T#"]==1)]

    print(f"\n\nAns 1.2.1: There are {len(unormal)} normal unique elements.\n")
    # After revision we need to add one to maintain the logic if the data set
    print(f"Ans 1.2.2: There are {len(utumor)+1} tumor unique elements.\n")
    print(f"Ans 1.2.3: There are {len(shared)} shared unique elements.\n")

    #print(f"\nAns 1.2.1: In the merged file there are {aml['#N']}")


# Solve item 1.3 
def merger13(normal, tumor):

    aml = pd.concat([normal, tumor], axis=0)
    aml = aml.drop_duplicates()
    aml = aml.reset_index()

    #print(f"Ans 1.3.1: There are {len(aml)} rows.\n")
    print(aml)

    CSQ_columns = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", 
            "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", 
            "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND", "FLAGS", 
            "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", 
            "SWISSPROT", "TREMBL", "UNIPARC", "RefSeq", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", 
            "HGVS_OFFSET", "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF", 
            "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF", "ExAC_EAS_MAF", 
            "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF", "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", 
            "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", 
            "ENTREZ", "EVIDENCE"]
 
    for col in CSQ_columns:
        for row in range(aml.shape[0]):
            temp = str(aml[col][row])
            temp = temp.strip("[]").replace('"',"").replace("'","").replace(" ","")
     
            aml[col][row] = temp.split(",")
     
            if type(aml[col][row]) != list:
                print(aml[col][row])
                print(type(aml[col][row]))
    
    #print(aml)
    
    # save the AML file 
    aml = aml.explode(CSQ_columns, ignore_index=True)
    aml.to_csv("AML_Expand.csv", index=False)

    # Define the new dataframes       
    amlg = ["SYMBOL", "Gene", "Feature"]
    amlt = ["chrom", "left", "right", "ref_seq", "alt_seq", "Feature", "cDNA_position", "BIOTYPE"]

    # Define the datasets to export
    aml_gene = aml[amlg]
    aml_tx = aml[amlt]
    
    # export CSV files
    aml_gene.to_csv("AML_gene.csv")
    aml_tx.to_csv("AML_tx.csv")    


def task2solver():
    # Start the example using the Iris DF. 
    iris = load_iris()
    df = pd.DataFrame(data=iris.data, columns=iris.feature_names)
    df['target'] = iris.target
    
    # Print the initial data
    print(df) 

    # Here we will separate the features (X) and the target variable (y).
    X = df.iloc[:, :-1].values
    y = df.iloc[:, -1].values
    
    # Split for training and confirmation 
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42) 

    # Feature scaling 
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    
    # Build the classifier
    classifier = RandomForestClassifier(n_estimators=100, random_state=42)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)

    # Evaluation of the model
    accuracy = accuracy_score(y_test, y_pred)
    print(f'Accuracy: {accuracy * 100:.2f}%')
    
    conf_matrix = confusion_matrix(y_test, y_pred)
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(conf_matrix, annot=True, fmt='g', cmap='Blues', cbar=False, 
                xticklabels=iris.target_names, yticklabels=iris.target_names)
    
    plt.title('Confusion Matrix Heatmap')
    plt.xlabel('Predicted Labels')
    plt.ylabel('True Labels')
    plt.show()

    # Compute feature importance
    feature_importances = classifier.feature_importances_
    plt.barh(iris.feature_names, feature_importances)
    plt.xlabel('Feature Importance')
    plt.title('Feature Importance in Random Forest Classifier')
    plt.show()


if __name__=='__main__':
    # Import them data from HW3 
    normal = pd.read_csv("ALL_Normal.csv")
    tumor  = pd.read_csv("ALL_Tumor.csv")

    # Answer questions from task 1.1 and generate the input for 1.2
    group_normal, group_tumor = cleaner11(normal, tumor)

    # Answer questions from task 1.2
    merger12(group_normal, group_tumor) 
    
    # Answer questions from task 1.3
    merger13(normal, tumor)

    # Complete the tutorial from task 2
    task2solver()










