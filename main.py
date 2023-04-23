import pandas as pd
import numpy as np

from ReadingFasta import Sequence
from LoadPrimer import PrimerOperation
from Bio.Blast import NCBIWWW
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
import math

class Genome():
    def __init__(self):
        self.country = ""
        self.date = ""
        self.sequence = ""
        self.length = 0
        self.id = ""
        self.virus_type = ""
        self.lineage = ""


def getIdentity(seqA,seqB):
    misMatch = 0
    for base1, base2 in zip(seqA, seqB):
        if (base1 != base2):
            misMatch += 1
    identity =1 - (misMatch / len(seqA))
    return identity
def getBindingAlign(sequence1,sequence2):
    first_string = ""
    middle_string = ""
    second_string = ""
    for seq1,seq2 in zip(sequence1,sequence2):
        if (seq1==seq2):
            middle_string+=seq1
        else:
            middle_string+="|"
        first_string += seq1
        second_string +=seq2
    result = first_string+"\n"+middle_string + "\n"+second_string
    return result
def doAlignAlongTheSequence(genome,primer):
    results = []
    bestIdentity = 0
    for position in range(0, genome.length):
        genomeSequence = genome.sequence[position:position + primer["length"]]
        identity =  getIdentity(genomeSequence,primer["sequence"])
        if (identity>=identify_thresould):
            bindingAlignment = getBindingAlign(genomeSequence,primer["sequence"])
            temp = {"HitStart":position,
                    "HitEnd": position + primer["length"],
                    "Identity":identity,
                    "PrimerID" : primer["id"],
                    "BindingAlignment" : bindingAlignment }
            results.append(temp)
        if identity > bestIdentity:
            bestIdentity = identity
    return results, bestIdentity
def CheckAmplicant(minusList,plusList, probeList):
    results = []
    founded = False
    for plus in plusList:
        for minus in minusList:
            amplicant = abs(minus["HitStart"]-plus["HitStart"])
            if amplicant<amplicant_length:
                for probe in probeList:
                    cond1 = probe["HitStart"]>plus["HitStart"] and probe["HitStart"]<minus["HitStart"]
                    if cond1 == True:
                        temp = {"Plus":plus,"Probe":probe,"Minus": minus,"Amplicant": amplicant}
                        results.append(temp)
                        founded = True
    return founded,results
def saveResult(resultDict,type=None):
    file_handler = open("result/result.txt","a")
    file_handler.write("HitStart    HitEnd      Identity    BindingAlignment    Type\n")
    for item in resultDict:
        file_handler.write("{}\t{}\t{}\t{}\t{}".format(
            str(item["HitStart"]).rjust(12),
            item["HitEnd"],
            item["Identity"],
            item["BindingAlignment"],
            type
        ))
        print(item)
    file_handler.close()
def saveResultInExcel(results):
    import xlsxwriter  # install xlsxwriter

    workbook = xlsxwriter.Workbook("result/result.xlsx")
    worksheet = workbook.add_worksheet("Primers")
    def writeHeader(col):
        worksheet.write(0, col, "ID")
        worksheet.write(0, col+1, "HitStart")
        worksheet.write(0, col+2, "HitEnd")
        worksheet.write(0, col+3, "Identity")
        worksheet.write(0, col+4, "Applicant")
        worksheet.write(0, col+5, "BindingAlignment")
    writeHeader(0)
    writeHeader(6)
    writeHeader(12)

    row = 1
    col = 0
    for item in results:
        col = 0
        worksheet.write(row,col , item["Plus"]["PrimerID"])
        worksheet.write(row,col+1 , item["Plus"]["HitStart"])
        worksheet.write(row,col+2 , item["Plus"]["HitEnd"])
        worksheet.write(row,col+3 , item["Plus"]["Identity"])
        worksheet.write(row,col+4 , item["Amplicant"])
        worksheet.write(row,col+5 , item["Plus"]["BindingAlignment"])

        col=6
        worksheet.write(row, col, item["Probe"]["PrimerID"])
        worksheet.write(row,col+1 , item["Probe"]["HitStart"])
        worksheet.write(row,col+2 , item["Probe"]["HitEnd"])
        worksheet.write(row,col+3 , item["Probe"]["Identity"])
        worksheet.write(row, col + 4, item["Amplicant"])
        worksheet.write(row,col+5 , item["Probe"]["BindingAlignment"])

        col = 12
        worksheet.write(row,col , item["Minus"]["PrimerID"])
        worksheet.write(row,col+1 , item["Minus"]["HitStart"])
        worksheet.write(row,col+2 , item["Minus"]["HitEnd"])
        worksheet.write(row,col+3 , item["Minus"]["Identity"])
        worksheet.write(row, col + 4, item["Amplicant"])
        worksheet.write(row,col+5 , item["Minus"]["BindingAlignment"])
        row+=1
    workbook.close()



def doAlignment(genome):
    for i,key in enumerate(primer_groups.keys()):
        primer = primer_groups[key]
        plus_Result,  bestIdentity = doAlignAlongTheSequence(genome, primer["PlusStrand"])
        probe_Result, bestIdentity = doAlignAlongTheSequence(genome, primer["ProbeStrand"])
        minus_Result, bestIdentity = doAlignAlongTheSequence(genome, primer["MinusStrand"])

        match,results = CheckAmplicant(minusList=minus_Result,plusList=plus_Result,probeList=probe_Result)
        if (match == True):
            print ("Save result in excel")
            saveResultInExcel(results)


def startAligning(filename):
    seq =Sequence()
    with open(filename,"r") as file:
        file.readline()
        while(True==True):
            line = file.readline()
            if not line:
                print("There is no sequence in file")
                break
            data = line.split(",")
            genome = Genome()
            genome.sequence = data[13]
            genome.length= len(genome.sequence)
            genome.id = data[2]

            match = doAlignment(genome)

            #print("{},{},{}", genome.id, data[5],match)
            print("{},{},{}".format(genome.id, data[5],match))

identify_thresould = 0.8
amplicant_length = 1000

primerOperation = PrimerOperation()
primer_groups,primers_id = primerOperation.loadPrimerFromFile("primer/primers.xlsx")
startAligning("clean/test_0.csv")


