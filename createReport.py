__author__ = 'christine'

import argparse
import csv
import pandas as pd
import numpy as np
import math

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

TCGAToGTEx = {'STAD': 'Stomach', 'KIRP': 'Kidney', 'THCA': 'Thyroid', 'PAAD': 'Pancreas', 'KICH': 'Kidney',
              'ESCA': 'Esophagus', 'CESC': 'Cervix Uteri', 'OV': 'Ovary', 'PCPG': 'Adrenal Gland',
              'BRCA': 'Breast', 'KIRC': 'Kidney', 'LUAD': 'Lung', 'LGG': 'Brain', 'BLCA': 'Bladder',
              'GBM': 'Brain', 'SKCM': 'Skin', 'ACC': 'Adrenal Gland', 'LIHC': 'Liver', 'TGCT': 'Testis',
              'COAD': 'Colon', 'LUSC': 'Lung', 'UCEC': 'Uterus', 'PRAD': 'Prostate', 'UCS': 'Uterus'}

convertEnsembl = {'ENSG00000145649': 'GZMA', 'ENSG00000180644': 'PRF1'}
cancerGenes = {'STAD': ['Cytolytic granule mediated cell apoptosis'], 'KIRP': ['Cytolytic granule mediated cell apoptosis'], 'THCA': ['Cytolytic granule mediated cell apoptosis'], 'PAAD': ['Cytolytic granule mediated cell apoptosis'], 'KICH': ['Cytolytic granule mediated cell apoptosis']}
cancerPathways = {'Cytolytic granule mediated cell apoptosis': "Effective natural anti-tumor immunity require a cytolytic response in which cytotoxic T lymphocytes and natural killer cells release..."}
pathwayGenes = {'Cytolytic granule mediated cell apoptosis': ['ENSG00000145649', 'ENSG00000180644']}

gtexData = pd.read_table('gtexTable.tsv', index_col=0)
tcgaData = pd.read_table('tcgaTable.tsv', index_col=0)

def openPatientResults():
    f = open(args.filename)
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    tableReader = csv.DictReader(f, fieldnames=header, delimiter="\t")
    return tableReader

def getStatus(gtexValue, tcgaValue, patientValue, gene):
    if patientValue > gtexValue and patientValue > tcgaValue:
        return "Status: Based on the results the patient is recommended for immunotherapy"
    else:
        return "Status: Based on the results the patient is not recommended for immunotherapy"

def getTPMValue(gene):
    gene = gene[:15]
    for row in openPatientResults():
        if gene == row["gene_id"][:15]:
            tpmValue = math.log(float(row["TPM"]) + .001, 2)
            return tpmValue

def getTCGA(gene, cancer):
    for key in tcgaData.index.values:
        if key[:15] == gene:
            meanTCGA = tcgaData.loc[key, cancer + " normal"]
            return str(meanTCGA)
        break

def getGTEx(gene, cancer):
    for key in gtexData.index.values:
        if key[:15] == gene:
            meanGTEx = gtexData.loc[key, TCGAToGTEx[cancer]]
            print(str(meanGTEx))
            return str(meanGTEx)

"""
@param: cancer type
writes an output file with information on each gene signature
compares patient's tpm value of gene signature with GTEx and TCGA normal tissue gene expression values
transforms patient's tpm value log2(x+.01) since GTEx and TCGA values are transformed
"""

def makeOutputFile(cancer):
    outputFile = open('report', 'w')
    openPatientResults()
    pathways = cancerGenes[cancer]
    for pathway in pathways:
        description = cancerPathways[pathway]
        genes = pathwayGenes[pathway]
        outputFile.write("Pathway: " + pathway + "\n")
        outputFile.write("Description of pathway: " + description + "\n")
        outputFile.write("Pathway genes: \n")
        outputFile.write("\tGENE\tGTEX_THRESHOLD\tTCGA_NORMAL_THRESHOLD\tOBSERVED\n")
        for gene in genes:
            print(gene)
            if cancer in TCGAToGTEx.keys():
                gtex = getGTEx(gene, cancer)
                print(gtex)
            else:
                gtex = "NA"
            outputFile.write("\t" + gene + "\t" + gtex + "\t" + getTCGA(cancer, gene) + "\t"
                                         + str(getTPMValue(gene)) + "\n")
        outputFile.write("\n")
        outputFile.write("\n")
    outputFile.close()


makeOutputFile('BRCA')

    #     gene = gene[:15]
    #     for row in openPatientResults():
    #         if gene == row["gene_id"][:15]:
    #             tpmValue = math.log(float(row["TPM"]) + .01, 2)
    #             outputFile.write("Pathway: " + "\n")
    #             outputFile.write("Description of pathway: \n")
    #             for key in tcgaData.index.values:
    #                 if key[:15] == gene:
    #                     meanTCGA = tcgaData.loc[key, cancer + " normal"]
    #                     if cancer in TCGAToGTEx.keys(): # if TCGA has comparable GTEx cancer
    #                         for key in gtexData.index.values:
    #                             if key[:15] == gene:
    #                                 meanGTEx = gtexData.loc[key, TCGAToGTEx[cancer]]
    #                                 outputFile.write(getStatus(meanTCGA, meanGTEx, ))
    #                                 outputFile.write("Pathway genes:\n")
    #                                 outputFile.write("\t" + gene + "\t" + str(meanGTEx) + "\t"
    #                                                  + str(meanTCGA) + "\t" + str(tpmValue) + "\n")
    #                                 outputFile.write("\tgene " + str(n) + "\tthreshold for gene " + str(n)
    #                                                  + " in gtex\tthreshold for gene " + str(n) +
    #                                                  " in TCGA normals\tobserved from patient file")
    #
    #                     else:
    #                         outputFile.write("\t" + gene + "\t" + str(meanTCGA) + "\t" + str(tpmValue) + "\n")
    #                         outputFile.write("\tgene " + str(n) + "\tthreshold for gene " + str(n) +
    #                                          " in TCGA normals\tobserved from patient file")
    #             outputFile.write("\n")
    #             outputFile.write("\n")
    #     n += 1
    # outputFile.close()

