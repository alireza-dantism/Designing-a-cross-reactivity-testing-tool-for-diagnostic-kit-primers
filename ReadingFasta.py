import pandas as pd
class Sequence():
    validLength = 29903
    validLengthRange = 500
    validBases = ['T','A','G','C']
    seperatorChar = "/"
    def __int__(self):
        self.id = ""
        self.header= ""
        self.epi_id = "NA"
        self.date = "NA"
        self.sequence =""
        self.country = "NA"
        self.continent = "NA"
        self.length = 0
        self.age = "NA"
        self.sex = "NA"
        self.host = "NA"
        self.submeteddate = "NA"
        self.virustype = "NA"
        self.lineage = "NA"
        self.error = ""
    def isLengthValid(self):
        beginRange = self.validLength-self.validLengthRange
        endRange = self.validLength + self.validLengthRange
        if ( self.length>= beginRange and self.length<=endRange):
            return True
        else:
            return False
    def isHuman(self):
        if self.host=="Human":
            return True
        return False
    def containN(self):
        for nocleotid in self.sequence:
            if (nocleotid not in self.validBases):
                return True
        return False
    def splitLine(self,line):
        sequence = Sequence()
        try:
            data = line.split(self.seperatorChar)
            sequence.country = data[0]
            sequence.sequence = data[1]
            sequence.length = data[2]
            sequence.id = data[3]
            sequence.epi_id = data[4]
            sequence.date = data[5]
            sequence.continent = data[6]
            sequence.sex = data[7]
            sequence.age = data[8]
            sequence.host = data[9]
            sequence.submeteddate = "25/12/156"
            sequence.error = False
        except:
            sequence.error = True
        return sequence
class Covid_FastA():
    def __init__(self,fileName,metaDataFileName=None , lineageFiles = None , outdir = "clean", splitRecordCount=5000000):
        self.splitRecordCount = splitRecordCount
        self.fileName = fileName
        self.outDir = outdir
        self.seperatorChar = " "
        if (metaDataFileName!=None):
            self.metaDataDF = pd.read_csv(metaDataFileName,sep='\t',index_col='strain')
        else :
            self.metaDataDF = None
        if (lineageFiles != None):
            self.lineageDF = pd.read_excel(lineageFiles)

            #self.metaDataDF.set_index('strain')
    def getRecordLength(self):
        with open(self.fileName,"r") as file:
            line1 = file.readline()
            line2 = file.readline()
            file.close()
        return len(line1+line2)
    def saveToFile(self,sequence, fileHandeler):
        fileHandeler.write(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(sequence.country,
                                                     sequence.length,
                                                     sequence.id,
                                                     sequence.epi_id,
                                                     sequence.date,
                                                     sequence.submeteddate,
                                                     sequence.continent,
                                                     sequence.sex,
                                                     sequence.age,
                                                     sequence.host,
                                                     sequence.virustype,
                                                     sequence.lineage,
                                                     sequence.header,
                                                     sequence.sequence,
                                                     ))
    def __isValidSequence(self,sequence):
        hasError = sequence.error
        isValidLength = sequence.isLengthValid()
        isHuman = sequence.isHuman()
        if(not hasError and isValidLength and isHuman):
            return True
        else:
            return False
    def splitCleanCovidFile(self,outputFileName,statisticalFile):

        totalCount = 0
        validCount = 0
        splitCount = 0
        fileCounter = 0

        fileWriter = open("{}/{}_{}.csv".format(self.outDir,outputFileName,fileCounter), "w")

        with open(self.fileName,"r",buffering=1) as fileReader:
            seq = Sequence()

            while (True==True):
                line = fileReader.readline()
                if (not line):
                    break
                line = fileReader.readline()
                sequence = seq.splitLine(line)

                totalCount+=1

                if (self.__isCleanSequence(sequence)):
                    validCount+=1
                    self.saveToFile(sequence,fileWriter)
                    splitCount += 1
                    if (splitCount>=self.splitRecordCount):
                        fileWriter.close()
                        splitCount=0
                        fileCounter+=1
                        fileWriter = open("{}/{}_{}.csv".format(self.outDir,outputFileName,fileCounter),"w")

        fileReader.close()
        fileWriter.close()

        file = open("{}_csv".format(statisticalFile),"w")
        file.write("Total : {}\n Valid :{} \n".format(totalCount,validCount))
    def readSequence(self,fileReader):
        line = fileReader.readline()
        line = line.replace("\n","")
        sequence = ""
        while(line.find(">")<0 ):
            sequence+= line.replace("\n","")
            line = fileReader.readline()
            if (not line):
                break
        return sequence, line
    def __splitHeader(self,header):

        sequence = Sequence()
        try:
            sequence.header = header[1:].replace("\n", "").replace (",","")
            head = sequence.header.split(self.seperatorChar)
            sequence.id = head[0][1:]
            sequence.epi_id = head[0][1:]
            sequence.country = "NA"
            sequence.date = "NA"
            sequence.continent = "NA"
            sequence.sex = "NA"
            sequence.age = "NA"
            sequence.host = "NA"
            sequence.submeteddate = "NA"
            sequence.virustype = "NA"
            sequence.lineage = "NA"
            sequence.error = False
        except:
            sequence.error = True
        return sequence
    def splitHeader(self,line):
        if (line.find(">")>=0):
            header = self.__splitHeader(line)
            return header
        return  None
    def find_Lineage(self,lineage):
        lineages = self.lineageDF.to_numpy()
        for item in lineages:

            if item[1] == lineage:
                return item[0]
        return "Other"

    def findMetaInfo(self,sequenceInfo):
        sequence = Sequence()
        sequence = sequenceInfo
        if (self.metaDataDF!=None):
            iloc = self.metaDataDF.index.get_loc(sequenceInfo.header)
            findedRow = self.metaDataDF.iloc[iloc]
            if (len(findedRow.values.shape)>1):
                findedRow = findedRow.iloc[0]
            sequence.sex = findedRow['sex']
            sequence.age = findedRow['age']
            sequence.country = findedRow['country']
            sequence.date = findedRow['date']
            sequence.submeteddate = findedRow['date_submitted']
            sequence.host = findedRow['host']
            sequence.continent = findedRow['region']
            sequence.epi_id = findedRow['gisaid_epi_isl']
            sequence.virustype = findedRow['virus']
            sequence.lineage = self.find_Lineage(findedRow['pangolin_lineage'])
        else:
            sequence.sex = ""
            sequence.age = ""
            sequence.country = ""
            sequence.date = ""
            sequence.submeteddate = ""
            sequence.host = ""
            sequence.continent = ""
            sequence.epi_id = ""
            sequence.virustype = ""
            sequence.lineage = ""

        return sequence

    def splitCleanCovidGenome(self, outputFileName, statisticalFile):
        totalCount = 0
        validCount = 0
        splitCount = 0
        fileCounter = 0

        fileWriter = open("{}/{}_{}.csv".format(self.outDir, outputFileName, fileCounter), "w")
        fileWriter.write("Country,Length,Id,EpiID,Date,SubmitDate,Continent,Sex,Age,Host,VirusType,Lineage,Header,Sequence\n")

        with open(self.fileName, "r", buffering=1) as fileReader:
            line = fileReader.readline()
            while (True == True):
                if (not line):
                    break
                sequneceInfo = self.splitHeader(line)
                sequence,lastline = self.readSequence(fileReader=fileReader)
                line = lastline
                sequneceInfo.sequence = sequence
                sequneceInfo.length = len(sequence)
                totalCount += 1
                sequneceInfo = self.findMetaInfo(sequneceInfo)

                try:
                    validCount += 1
                    self.saveToFile(sequneceInfo, fileWriter)
                    splitCount += 1
                    if (splitCount >= self.splitRecordCount):
                        fileWriter.close()
                        splitCount = 0
                        fileCounter += 1
                        fileWriter = open("{}/{}_{}.csv".format(self.outDir, outputFileName, fileCounter), "w")
                        fileWriter.write("Country,Sequence,Length,Id,EpiID,Date,SubmitDate,Continent,Sex,Age,Host,VirusType,Lineage,Header\n")
                except Exception as e:
                    print (sequneceInfo.epi_id)
                    continue
        fileReader.close()
        fileWriter.close()

        file = open("{}/{}.csv".format(self.outDir,statisticalFile), "w")
        file.write("Total : {}\n Valid :{} \n".format(totalCount, validCount))

sequence_file = "data/test.fasta"
if (__name__ == '__main__'):
    covid = Covid_FastA(sequence_file ,
                        metaDataFileName =None,
                        lineageFiles = None,
                        outdir= "clean")
    covid.splitCleanCovidGenome("test","statistic")
