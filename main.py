import sys
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt
from datetime import datetime

## Create list of lists of 0s (no SNP) and 1s (SNP at location) in relation to reference Sars-CoV-2 genome
# Matrix: rows are genetic sequence positions, columns are different sequences
# 1s represent SNPs and indels at their first position
# This script will make an insertion a 1 at its first position

refGenomeFile = sys.argv[1]
inputs = sys.argv[2:] or sys.stdin
refGenome = ""
with open(refGenomeFile) as refFile:
    for line in refFile:
        if line[0] != ">":
            refGenome += line

indSeq = []
allSeqs = []
fileDict = {}
fileName = ""
fileOrder = []
rCount = 0
rows = []  #SM
columns = []  #SM
data = []
lenInputs = 0

# REFERENCE IS NOT IN MATRIX, JUST AS ITS OWN STRING FROM REF_FILE
index = 1
col = 0
for inputLine in inputs:
    lenInputs += 1
    lineSplit = inputLine.split()
    if lineSplit[-1] != fileName: # if we see a new file
        if indSeq:
            col += 1 # SM
            for i in range(len(refGenome) - index):  # add zeros to the remainder of the previous file
                indSeq.append(0)
            allSeqs.append(indSeq)
            fileOrder.append(fileName)
            fileDict[fileName] = rCount
            rCount = 0
        fileName = lineSplit[-1]
        indSeq = []
        index = 0
    rows.append(int(lineSplit[0]) - 1)  # SM
    columns.append(col)  # SM
    data.append(1)  # SM
    # check for position of mutation and add 0s and 1s accordingly
    if lineSplit[0] == index:
        indSeq.append(1)
    else:
        j = 0
        for i in range(int(lineSplit[0]) - index):
            indSeq.append(0) # fill in the gaps with 0s
            j += 1
        #index += (int(lineSplit[0]) - index)
        index += j
        if int(lineSplit[0]) == index:
            indSeq.append(1)
            index += 1
    rCount += int(lineSplit[6])

for i in range(len(refGenome) - index):  # add zeros to the remainder of the last file
    indSeq.append(0)
allSeqs.append(indSeq)
fileOrder.append(fileName)
fileDict[fileName] = rCount

# use np to make matrix from list of lists
completeMatrix = np.array(allSeqs)
#np.savetxt("Sars-CoV-2_Matrix_Data_EG", completeMatrix, fmt="%d")
#np.savetxt("Sars-CoV-2_FileOrder", fileOrder, fmt="%s")

fileDictionary = open('Sars_fileDict', 'wt')
for key, value in fileDict.items():
    fileDictionary.write(str(key) + "\tTotal [R] Value: " + str(value) + "\n")
fileDictionary.close()

# Sort files by date
fileByDate = []
fileDatesSorted = []
for i in range(len(fileOrder)):
    splitStr = fileOrder[i].split("|")
    fileDate = splitStr[2]
    fileByDate.append(fileDate)
fileByDate.sort(key=lambda date: datetime.strptime(date, "%Y-%m-%d"))
for i in range(len(fileByDate)):
    for j in range(len(fileOrder)):
        if fileByDate[i] in fileOrder[j] and fileOrder[j] not in fileDatesSorted:
            fileDatesSorted.append(fileOrder[j])

## PANDAS
#seqsSeries = pd.Series(allSeqs)
#fileSeries = pd.Series(fileOrder)
#concatMatrix = pd.concat([fileSeries, seqsSeries], axis=1)
#concatMatrix.to_csv("PandasMatrix.csv", sep="\t")

## SPARSE MATRIX
rows = np.array(rows)
columns = np.array(columns)
data = np.array(data)
#sparseMatrix = csr_matrix((data, (columns, rows)), shape=(len(inputs), len(refGenome)))
sparseMatrix = sparse.coo_matrix((data, (columns, rows)), shape=(lenInputs, len(refGenome)))
#dfMatrix = pd.DataFrame(sparseMatrix)
#dfMatrix.to_csv("sparseMatrix")
#np.savetxt("sparseMatrix", (sparseMatrix.col, sparseMatrix.row), fmt="%d")

## PLOTS
# Plot of num SNPS per genome
myDict = {}
myDate = ""
filesInDict = []
for i in range(len(sparseMatrix.row)):
    toSplit = fileDatesSorted[sparseMatrix.row[i]].split("|")
    myDate = toSplit[2]
    if fileDatesSorted[sparseMatrix.row[i]] in myDict:
        myDict[fileDatesSorted[sparseMatrix.row[i]]] = myDict[fileDatesSorted[sparseMatrix.row[i]]] + 1
    else:
        myDict[fileDatesSorted[sparseMatrix.row[i]]] = 1

genomes = []
counter = 2
for file in myDict.keys():
    toSplit = file.split("|")
    myDate = toSplit[2]
    if myDate in genomes:
        genomes.append(myDate + "." + str(counter))
        counter += 1
    else:
        genomes.append(myDate)
frequencies = list(myDict.values())

fig, ax = plt.subplots()
plt.bar(genomes, frequencies)
ax.xaxis.set_major_locator(plt.AutoLocator())
plt.xlabel("Genome (by date)")
plt.ylabel("Number of SNPS")
plt.xticks(rotation=60)
plt.grid()
plt.savefig("NumSnpsPerGenome.png", bbox_inches="tight")

# Plot of num genomes with that certain SNP
myDict2 = {}
for i in range(len(sparseMatrix.col)):
    if sparseMatrix.col[i] in myDict2:
        myDict2[sparseMatrix.col[i]] = myDict2[sparseMatrix.col[i]] + 1
    else:
        myDict2[sparseMatrix.col[i]] = 1

snpPos = list(myDict2.keys())
freqOfSNP = list(myDict2.values())
plt.bar(snpPos, freqOfSNP, width=50)
plt.xlabel("Position of SNP")
plt.xticks([0,5000,10000,15000, 20000, 25000, 30000], ["0", "5000", "10000", "15000", "20000", "25000", "30000"], rotation=0)
plt.ylabel("Number of sequences with SNP")
plt.savefig("FreqOfSNPs.png")