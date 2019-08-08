__author__ = 'christine'

import numpy as np
import pandas as pd
import argparse
import csv
import matplotlib.pyplot as plt

#read file using argparse
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

#read in table
f = open('SraRunTable.txt')
csv_f = csv.reader(f, delimiter="\t")
header = next(csv_f)
tableReader = csv.DictReader(f, fieldnames=header, delimiter="\t")

#create dictionary in format {cancer type:{ normal: {sample}}}
tissues = {}
for row in tableReader:
    cancerType = row["histological_type_s"]
    isTumor = row["is_tumor_s"]
    sample = row["Sample_Name_s"]
    if cancerType not in tissues.keys():
        tissues[cancerType] = [sample]
    else:
        tissues[cancerType].append(sample)

def makeDataFrame(geneArray):
    with open(args.filename, 'r') as f:
        csvData = csv.reader(f, delimiter="\t")
        dataDict = {}
        i = 0
        for row in csvData:
            if i == 0:
                columnHeaders = row[1:]
            else:
                if any(row[0][0:15] in s for s in geneArray):  # ignore "." and numbers after b/c versions may differ
                    id = row[0]
                    dataDict[id] = row[1:]
            i += 1
    data = pd.DataFrame(dataDict, index=columnHeaders)
    return data

def getMean(df):
    return str(df.stack().mean())

def getMedian(df):
    return str(df.stack().median())

def getStdev(df):
    return str(df.stack().std())

def practicePlot(myDict, paired, gene):
    if paired == False:
        colors = ['lightblue', '#ee8080', '#ffbd86', 'magenta', '#52a4ff', '#aabbcc', "#2cddbe", "#ae1534",
        "#691744", "#e4fdc4", "#ffad30", "#008080", '#ffc3a0', '#ff6666', '#c6e2ff', '#808080', "#3b5998",
        "#c39797", "#ffdab9", "#F63054", "#b0e0e6", "#999999", "#9ea6cf", "#c96dc2", "#03b18b", "#d84d92",
        "#98baa4", "#fe9153", "#d4b47a", "#54b9a7", "#fc3373"]
        sortedDict = sorted(myDict, key=lambda x:np.mean(myDict[x]))
        plotDataList = []
        for key in sortedDict:
            plotDataList.append(myDict[key])
    if paired == True:
        sortedDict = sorted(myDict.keys(), key=lambda x:np.mean(myDict.keys[x]))
    fig = plt.figure()
    fig.suptitle(gene)
    bp = plt.boxplot(plotDataList, labels=sortedDict, vert=False, patch_artist=True)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='blue', marker='o', markersize=3)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    plt.show()

def writeTsv(dict):
    my_df = pd.DataFrame(dict)
    my_df = my_df.transpose()
    my_df.to_csv("gtexTable.tsv", sep="\t")

def getValues(geneArray):
    with open(args.filename, 'r') as f:
        csvData = csv.reader(f, delimiter="\t")
        dataDict = {}
        i = 0
        for row in csvData:
            if i == 0:
                columnHeaders = row[1:]
                print(len(columnHeaders))
            else:
                if any(row[0][0:15] in s for s in geneArray):
                    id = row[0]
                    dataDict[id] = row[1:]
            i += 1
    data = pd.DataFrame(dataDict, index=columnHeaders)
    plotDict = {}
    meanPlotDict = {}
    cancers = tissues.keys()
    for gene in list(data):
        for cancer in cancers:
            if cancer != "<not provided>" and cancer != "Bone Marrow":
                tissueArray = tissues[cancer]
                tissueVals = data.loc[tissueArray, [gene]]
                tissueValsGene = tissueVals.dropna()
                tissueValsGene = tissueValsGene.apply(pd.to_numeric)
                if gene not in meanPlotDict.keys():
                    meanGene = {}
                    meanGene[cancer] = getMean(tissueValsGene)
                    meanPlotDict[gene] = meanGene
                else:
                    meanGene = {}
                    meanGene[cancer] = getMean(tissueValsGene)
                    meanPlotDict[gene].update(meanGene)
                myArr = tissueValsGene[gene].tolist()
                plotDict[cancer] = myArr
        practicePlot(plotDict, False, gene)
    # writeTsv(meanPlotDict)

lines = [line.rstrip('\n') for line in open('mhc_genes')]
getValues(lines)