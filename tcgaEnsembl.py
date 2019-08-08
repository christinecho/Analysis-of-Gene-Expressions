__author__ = 'christine'

import argparse
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('tablename')
args = parser.parse_args()

# read in table with tumor and cancer information
f = open(args.tablename)
csv_f = csv.reader(f, delimiter="\t")
header = next(csv_f)
tableReader = csv.DictReader(f, fieldnames=header, delimiter="\t")

"""
create structure to house all needed information: cancer, tumor/normal, tumor type (if relevant), tissue name
basic structure: {cancer: {normal: [], tumor: {primary: [], recurring: []}}}
tissues ending in "01" or "03" are primary. tissues ending in "02" or "04" are recurrent
tissues ending in "11" are normal
"""
tissues = {}
for row in tableReader:
    cancerType = row["disease"]
    sample = row["barcode"][:15]  # ignore "." and numbers after because gene versions may change
    lastChars = int(sample[-2:])
    if lastChars <= 11 and lastChars != 10:
        if lastChars == 11:
            key = "normal"
        else:
            key = "tumor"
            if lastChars == 1 or lastChars == 3:
                innerKey = "primary"
            elif lastChars == 2 or lastChars == 4:
                innerKey = "recurrent"
            else:
                continue  # check this
        if cancerType not in tissues.keys():
            if key == "tumor":
                tumor = {}
                tumorType = {}
                tumorType[innerKey] = [sample]
                tumor[key] = tumorType
                tissues[cancerType] = tumor
            else:  # if key is normal
                tumor = {}
                tumor[key] = [sample]
                tissues[cancerType] = tumor
        else:
            if key == "tumor":
                if key not in tissues[cancerType].keys():
                    tumor = {}
                    tumorType = {}
                    tumorType[innerKey] = [sample]
                    tumor[key] = tumorType
                    tissues[cancerType].update(tumor)
                else:
                    if innerKey not in tissues[cancerType][key].keys():
                        tumorType = {}
                        tumorType[innerKey] = [sample]
                        tissues[cancerType][key].update(tumorType)
                    else:
                        tissues[cancerType][key][innerKey].append(sample)
            else:  # if key is normal
                if key not in tissues[cancerType].keys():
                    tumor = {}
                    tumor[key] = [sample]
                    tissues[cancerType].update(tumor)
                else:
                    tissues[cancerType][key].append(sample)

# write tsv file of mean gene expression values from dictionary of dictionaries
def writeTsv(dict):
    my_df = pd.DataFrame(dict)
    my_df = my_df.transpose()  # change indexes to columns and vice versa
    my_df.to_csv("tcgaTable.tsv", sep="\t")

# plots the tumor and normal tissues of a gene as a boxplot for each histological type
def practicePlot(myDict, paired, gene):
    boxColors = ['lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', 'lightblue', '#ee8080', ]
    if paired == False:  # GTEx
        sortedDict = sorted(myDict)
        plotDataList = []
        for key in sortedDict:
            plotDataList.append(myDict[key])
        i = 0
        # for key in sortedDict:
        #     if i % 2 == 0:
        #         print(key[0:4])
        #         sortedDict[key[0:4]] = sortedDict.pop(key)
        #     else:
        #         sortedDict[""] = sortedDict.pop(key)
        #     i += 1
    fig = plt.figure()
    fig.suptitle(gene)
    bp = plt.boxplot(plotDataList, labels=sortedDict, vert=False, patch_artist=True)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='blue', marker='o', markersize=3)
    for patch, color in zip(bp['boxes'], boxColors):
        patch.set_facecolor(color)

    #add a legend
    plt.figtext(0.80, 0.08, 'Normal Tissues', backgroundcolor=boxColors[0], color='black', weight='roman', size='x-small')
    plt.figtext(0.80, 0.045, 'Tumor Tissues', backgroundcolor=boxColors[1], color='black', weight='roman', size='x-small')
    plt.show()

def getMean(df):
    return str(df.stack().mean())

"""
@param: genes we want to put in tsv file
creates data frame with genes as column headers
only create DataFrame with desired genes
returns dictionary of dictionaries. outer dictionary has gene as key
inner dictionary has cancer + tumor as key and mean gene expression as value
use this function to create tsv files of mean gene expression values
"""
def getMeanDict(geneArray):
    with open(args.filename, 'r') as f:
        csvData = csv.reader(f, delimiter="\t")
        dataDict = {}
        i = 0
        for row in csvData:
            if i == 0:
                columnHeaders = row[1:]
            else:
                if any(row[0][0:15] in s for s in geneArray):
                    id = row[0]
                    dataDict[id] = row[1:]
            i += 1
    data = pd.DataFrame(dataDict, index=columnHeaders)
    meanPlotDict = {}
    cancers = tissues.keys()
    print(cancers)
    for gene in list(data):
        for cancer in cancers:
            for key in tissues[cancer]:
                if key == "normal":
                    tissueArray = tissues[cancer][key]
                else:
                    for innerKey in tissues[cancer][key]:
                        if innerKey == "primary":
                            tissueArray = tissues[cancer][key][innerKey]
                tissueVals = data.loc[tissueArray, [gene]]
                tissueVals = tissueVals.dropna()
                tissueVals = tissueVals.apply(pd.to_numeric)
                if gene not in meanPlotDict.keys():
                    meanGene = {}
                    meanGene[cancer + " " + key] = getMean(tissueVals)
                    meanPlotDict[gene] = meanGene
                else:
                    meanGene = {}
                    meanGene[cancer + " " + key] = getMean(tissueVals)
                    meanPlotDict[gene].update(meanGene)
    writeTsv(meanPlotDict)

"""
@param: genes we want to put in tsv file
creates data frame with genes as column headers
only create DataFrame with desired genes
returns dictionary with cancer type + tumor as key and mean gene expression value as value
use this function to make dictionary to pass in graphing function
"""

def getPlotDict(geneArray):
    with open(args.filename, 'r') as f:
        csvData = csv.reader(f, delimiter="\t")
        dataDict = {}
        i = 0
        for row in csvData:
            if i == 0:
                columnHeaders = row[1:]
            else:
                if any(row[0][0:15] in s for s in geneArray):
                    id = row[0]
                    dataDict[id] = row[1:]
            i += 1
    data = pd.DataFrame(dataDict, index=columnHeaders)
    cancers = tissues.keys()
    genes = list(data)
    tissueWork = tissues
    i = 0
    for gene in genes:
        if i == 0:
            for cancer in cancers:
                for key in tissues[cancer]:
                    if key == "normal":
                        tissueArray = tissues[cancer][key]
                        tissueVals = data.loc[tissueArray, gene]
                        tissueVals = tissueVals.dropna()
                        myArr = tissueVals.values.tolist()
                        myArr = [float(i) for i in myArr]
                        tissueWork[cancer][key] = myArr
                    else:
                        for innerKey in tissueWork[cancer][key]:
                            tissueArray = tissues[cancer][key][innerKey]
                            tissueVals = data.loc[tissueArray, gene]
                            tissueVals = tissueVals.dropna()
                            myArr = tissueVals.values.tolist()
                            myArr = [float(i) for i in myArr]
                            tissueWork[cancer][key][innerKey] = myArr
        i += 1
    for cancer in tissueWork.keys():
        if "normal" in tissueWork[cancer].keys():
            print(cancer)
            p = stats.ttest_ind(tissueWork[cancer]['tumor']['primary'], tissueWork[cancer]['normal'])
            print(p)
        # practicePlot(tissueWork, True, gene)

# get the mean expression values for PRF1 and GZMA
# can pass the dictionary to the practicePlot function to create graph
def getCytolyticActivityDict(geneArray):
    with open(args.filename, 'r') as f:
        csvData = csv.reader(f, delimiter="\t")
        dataDict = {}
        i = 0
        for row in csvData:
            if i == 0:
                columnHeaders = row[1:]
            else:
                if any(row[0][0:15] in s for s in geneArray):
                    id = row[0]
                    dataDict[id] = row[1:]
            i += 1
    data = pd.DataFrame(dataDict, index=columnHeaders)
    cancers = tissues.keys()
    genes = list(data)
    cytolDict = {}
    cancerWithNormals = []
    for cancer in tissues.keys():
        if "normal" in tissues[cancer].keys():
            cancerWithNormals.append(cancer)
    for gene in genes:
        for cancer in cancers:
            if cancer in cancerWithNormals:
                for key in tissues[cancer]:
                    if key == "normal":
                        newKey = cancer + " Normal"
                        tissueArray = tissues[cancer][key]
                        tissueVals = data.loc[tissueArray, gene]
                        tissueVals = tissueVals.dropna()
                        myArr = tissueVals.values.tolist()
                        myArr = [float(i) for i in myArr]
                        if newKey in cytolDict.keys():
                            cytolDict[newKey].extend(myArr)
                        else:
                            cytolDict[newKey] = myArr
                    else:
                        for innerKey in tissues[cancer][key]:
                            if innerKey == "primary":
                                newKey = cancer + " Tumor"
                                tissueArray = tissues[cancer][key][innerKey]
                                tissueVals = data.loc[tissueArray, gene]
                                tissueVals = tissueVals.dropna()
                                myArr = tissueVals.values.tolist()
                                myArr = [float(i) for i in myArr]
                                if newKey in cytolDict.keys():
                                    cytolDict[newKey].extend(myArr)
                                else:
                                    cytolDict[newKey] = myArr
    # for cancer in tissues.keys():
    #     if "normal" in tissues[cancer].keys():
    #         p = stats.ttest_ind(cytolDict[cancer + " Tumor"], cytolDict[cancer + ' Normal'])
    #         print(p)
    for i in cancerWithNormals:
        if np.mean(cytolDict[i + " Normal"]) > np.mean(cytolDict[i + " Tumor"]):
            cytolDict.pop(i + " Normal")
            cytolDict.pop(i + " Tumor")
    practicePlot(cytolDict, False, "PRF1 and GZMA")


#get genes in file
lines = [line.rstrip('\n') for line in open('mhc_genes')]
getCytolyticActivityDict(["ENSG00000145649", "ENSG00000180644"])