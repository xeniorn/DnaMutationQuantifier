

myFile = "RP272_cdna_ko.ab1"

class sequencing_locus:
    def __init__(self, a, c, g, t):
        val_sum = a + c + g + t
        self.a = a / val_sum
        self.c = c / val_sum
        self.g = g / val_sum
        self.t = t / val_sum
        
    def __getitem__(self, key):
        if key == "A":
            return self.a
        elif key == "C":
            return self.c
        elif key == "G":
            return self.g
        elif key == "T":
            return self.t

class normalized_trace:
    def __init__(self, abi_record):
        
        from Bio import SeqIO
        
        record = SeqIO.read(myFile, "abi")
        
        list(record.annotations.keys())
        list(record.annotations["abif_raw"].keys())
        
        
        from collections import defaultdict
        
        order = record.annotations['abif_raw']['FWO_1']
        
        base_to_data = {
        order[0:1].decode(): "DATA9",
        order[1:2].decode(): "DATA10",
        order[2:3].decode(): "DATA11",
        order[3:4].decode(): "DATA12"
        }
        
        channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
        trace = defaultdict(list)
        for c in channels:
            trace[c] = record.annotations["abif_raw"][c]
        
        indices = record.annotations["abif_raw"]["PLOC2"]
        primary_basecall = record.annotations["abif_raw"]["PBAS2"]	
        
        self.sequence = primary_basecall.decode()
        
        seqlen = len(self.sequence)
        
        complex_trace = []
        
        for index in range(seqlen):
            peak_index = indices[index]
            print(str(index) + " " + str(peak_index))
            vA = trace[base_to_data["A"]][peak_index]
            vC = trace[base_to_data["C"]][peak_index]
            vG = trace[base_to_data["G"]][peak_index]
            vT = trace[base_to_data["T"]][peak_index]
            seqloc = sequencing_locus(vA,vC,vG,vT)
            complex_trace.append(seqloc)
        
        self.complex_trace = complex_trace
        
    def __getitem__(self, key):
        return self.complex_trace[key]
    
    def calc_match(self, search_seq, template_start_index):
        seqlen = len(search_seq)
        
        if template_start_index + len(search_seq) > len(self.sequence):
            return 0
        
        score=0
        for probe_index in range(seqlen):
            template_index = template_start_index + probe_index
            score+=self[template_index][search_seq[probe_index]]
        return score
    
    def calc_matches(self, search_seq):
        
        seqlen = len(self.complex_trace)
        scores=[self.calc_match(search_seq, index) for index in range(seqlen)]
        
        return scores
    
    def get_top_match(self, search_seq):
        
        seqlen = len(self.complex_trace)
        maxindex=0
        maxscore=0
        for index in range(1 + seqlen - len(search_seq)):
            score = self.calc_match(search_seq, index)
            if score > maxscore:
                maxscore=score
                maxindex=index
        return maxindex

nt = normalized_trace("a")

nt[0]["G"]


import matplotlib.pyplot as plt
import numpy as np



la = np.array([i.a for i in complex_trace])
lc = np.array([i.c for i in complex_trace])
lt = np.array([i.t for i in complex_trace])
lg = np.array([i.g for i in complex_trace])
    
plt.bar(range(seqlen), lc, color="blue")
plt.bar(range(seqlen), la, bottom=lc      , color="red")
plt.bar(range(seqlen), lt, bottom=lc+la   , color="green")
plt.bar(range(seqlen), lg, bottom=lc+la+lt, color="black")
plt.show()




ax=plt.gca()
ax.set_xlim(0,7000)
	
plt.plot(trace[base_to_data["C"]], color="blue", linewidth=0.5)
plt.plot(trace[base_to_data["A"]], color="red", linewidth=0.5)
plt.plot(trace[base_to_data["T"]], color="green", linewidth=0.5)
plt.plot(trace[base_to_data["G"]], color="black", linewidth=0.5)
plt.show()
	
	