
# coding: utf-8

# In[1]:


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import random, sys
import numpy as np

#Set the frameshifts and stop-codons presence treatment (default = true)
try: 
    fr = int(sys.argv[4])
except:
    fr = 1
try:
    stop = int(sys.argv[3])
except:
    stop = 1
#Open fasta with V, D and J genes
handle = open('./V.fasta') 
records = list(SeqIO.parse(handle, "fasta"))
V = [str(records[x].seq) for x in range(len(records)) ]
handle.close()

handle = open('./D.fas')
records = list(SeqIO.parse(handle, "fasta"))
D = [str(records[x].seq) for x in range(len(records)) ]
handle.close()

handle = open('./J.fas')
records = list(SeqIO.parse(handle, "fasta"))
J = [str(records[x].seq) for x in range(len(records)) ]
handle.close()

nucleotides = 'agtc'
leng = pd.DataFrame(columns =["SeqID", "N1 length"]) #dataframe for lengths
annotation = pd.DataFrame(columns = ["SeqID", "V-REGION", "P-REGION V 3'", "N1-REGION","P-REGION D5'", "D-REGION", "P-REGION D3'",
                  "N2-REGION", "P-REGION J 5'","J-REGION"]) #dataframe for annotation
#Set initial parameters of the N1-zone formation
pp = 0.5 #p-nucleotides
pn = 0.1 #n-nucleotides
pe = 0.9 #exonuclease activity


def pal_3 (sequence): #3' palindrom formation
    pal = sequence[-np.random.geometric(pp):]
    my_seq = str(Seq(pal).reverse_complement())
    return my_seq

def pal_5 (sequence): #5' palindrom formation
    pal = sequence[:np.random.geometric(pp)]
    my_seq = str(Seq(pal).reverse_complement())
    return my_seq

def n_region (): #n region formation
    n_reg = ''
    for n in range(np.random.geometric(pn)):
        n_reg += random.choice(nucleotides)
    return n_reg

def exonuclease(seq): #exonucleases
    new_seq = ''
    for n in range(len(seq)):
        if np.random.binomial(n = 1, p = pe) == 1:
            new_seq += seq[n]
        else:
            continue
    return(new_seq)
            
        
for n in range (int(sys.argv[1])): #generation of antibodies
    V1, D1, J1  = random.choice(V), random.choice(D), random.choice(J)
    p1, n1, p2, p3, n2, p4 = exonuclease(pal_3(V1)), exonuclease(n_region()), exonuclease(pal_5(D1)), exonuclease(pal_3(D1)), exonuclease(n_region()), exonuclease(pal_5(J1))
    V_dom = V1 + p1 + n1 + p2 + D1 + p3 + n2 + p4 + J1
    len_n2 = len(n2)
    if fr == 1:
        if len(V_dom) % 3 == 1: #frameshift treatment
            V_dom = V_dom.replace(V_dom[len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + 1:len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + len_n2],
                                  V_dom[len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + 1:len(V1) + len(p1) + len(n1) + + len(p2) + len(D1) + len(p3) + len_n2] + random.choice(nucleotides) + random.choice(nucleotides))
            len_n2 +=2
        elif len(V_dom) % 3 == 2:
            V_dom = V_dom.replace(V_dom[len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + 1:len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + len_n2],
                                  V_dom[len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + 1:len(V1) + len(p1) + len(n1) + + len(p2) + len(D1) + len(p3) + len_n2] + random.choice(nucleotides))
            len_n2 +=1
        else:
            pass
    annotation.append(pd.DataFrame([["seq" + str(n), len(V1),len(V1) + len(p1), len(V1) + len(p1) + len(n1),
                                    len(V1) + len(p1) + len(n1) + len(p2), len(V1) + len(p1) + len(n1) + len(p2) + len(D1), len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3),
                                    len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + len_n2, len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + len_n2 + len(p4), len(V1) + len(p1) + len(n1) + len(p2) + len(D1) + len(p3) + len_n2 + len(p4) + len(J1)]], 
                                   columns = ["SeqID", "V-REGION", "P-REGION V 3'", "N1-REGION","P-REGION D5'", "D-REGION", "P-REGION D3'","N2-REGION", "P-REGION J 5'","J-REGION"]))

                
    leng.append(pd.DataFrame([["seq"+str(n), len(p1) + len(n1) + len(p2)]], columns =["SeqID", "N1 length"]))
    if stop == 1:
        V_dom = V_dom.replace("tag", "tac")
        V_dom = V_dom.replace("taa", "tac")
        V_dom = V_dom.replace("tga", "tgc")
    var = SeqRecord(Seq(V_dom), id = str(n) + "seq")
    output_handle = open(sys.argv[2] +".fasta", "a")
    SeqIO.write(var, output_handle, "fasta"),
    output_handle.close()
    n+=1
 
annotation.to_csv("./Annotations_" +  str(sys.argv[2]) + ".csv", encoding='utf-8', index=True)
leng.to_csv("./length_N1_" + str(sys.argv[2]) +".csv", encoding='utf-8', index=True) 

