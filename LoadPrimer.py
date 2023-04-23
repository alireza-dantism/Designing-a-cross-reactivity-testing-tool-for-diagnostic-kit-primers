import pandas as pd
import numpy as np

class Primer():
    def __int__(self):
        self.id = ""
        self.primers = list()
        self.lengths = list()
        self.startPosition = 0
        self.nsp = ""
        self.group = ""
        self.subgroup = ""
        self.type = ""
        self.temperature = ""


        self.alignScore=0
        self.exactMatch =0
        self.partialMatch =0
        self.totalSequence = 0

class PrimerOperation():
    def __init__(self):
        self.primer_dict = dict()
    def loadPrimerFromFile(self,filename):
        df = pd.read_excel(filename)

        primer_dict = dict()
        groups = df.groupby(["Group"])

        for group in groups:
            primer_info ={}
            subgroups = group[1].groupby("SubGroup")

            for item in subgroups:

                primers = item[1].to_numpy()

                temp_primer = Primer()
                sub_group = primers[0,4]
                temp_dict = {
                                "sequence":primers[0,5],
                                "length":len(primers[0,5]),
                                "position":0,
                                "id":primers[0,1],
                                "gen": "fullGenome",
                             }

                primer_info[sub_group] = temp_dict
            primer_dict[group[0]] = primer_info
        primers = df.to_numpy()[:,0]
        return primer_dict,primers


