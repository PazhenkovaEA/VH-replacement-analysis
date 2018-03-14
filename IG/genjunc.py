from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os, random, sys

try:
    fr = int(sys.argv[4])
except:
    fr = 1
try:
    stop = int(sys.argv[3])
except:
    stop = 1
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
if not os.path.exists('Annotations'):
    os.makedirs('Annotations')

nucleotides = 'agtc'

def pal_3 (sequence):
    pal = sequence[-random.randint(1, 3):]
    my_seq = str(Seq(pal).reverse_complement())
    return my_seq

def pal_5 (sequence):
    pal = sequence[:random.randint(1, 3)]
    my_seq = str(Seq(pal).reverse_complement())
    return my_seq

def n_region ():
    n_reg = ''
    for n in range(random.randint(5, 20)):
        n_reg += random.choice(nucleotides)
    return n_reg
for n in range (int(sys.argv[1])):
    V1, D1, J1  = random.choice(V), random.choice(D), random.choice(J)
    p1, n1, p2, p3, n2, p4 = pal_3(V1), n_region(), pal_5(D1), pal_3(D1), n_region(), pal_5(J1)
    V_dom = V1 + p1 + n1 + p2 + D1 + p3 + n2 + p4 + J1
    len_n1 = len(n1)
    if fr == 1:
        if len(V_dom) % 3 == 1: #лечим сдвиг рамки считывания
            V_dom = V_dom.replace(V_dom[len(V1) + len(p1) + 1:len(V1) + len(p1) + len(n1)],
                                  V_dom[len(V1) + len(p1) + 1:len(V1) + len(p1) + len(n1)] + random.choice(nucleotides) + random.choice(nucleotides))
            len_n1 +=2
        elif len(V_dom) % 3 == 2:
            V_dom = V_dom.replace(V_dom[len(V1) + len(p1) + 1:len(V1) + len(p1) + len(n1)],
                                  V_dom[len(V1) + len(p1) + 1:len(V1) + len(p1) + len(n1)] + random.choice(nucleotides))
            len_n1 +=1
        else:
            pass
    annotation = {"V-REGION" : [1, len(V1)], "P-REGION V 3'" : [len(V1)+1,len(V1) + len(p1)], "N1-REGION" : [len(V1) + len(p1) +1 , len(V1) + len(p1) + len_n1],
                  "P-REGION D5'": [len(V1) + len(p1) + len_n1 + 1, len(V1) + len(p1) + len_n1 + len(p2)], "D-REGION": [len(V1) + len(p1) + len_n1 + len(p2) + 1, len(V1) + len(p1) + len_n1 + len(p2) + len(D1)],
                  "P-REGION D3'" : [len(V1) + len(p1) + len_n1 + len(p2) + len(D1) +1, len(V1) + len(p1) + len_n1 + len(p2) + len(D1) + len(p3)],
                  "N2-REGION": [len(V1) + len(p1) + len_n1 + len(p2) + len(D1) + len(p3) +1, len(V1) + len(p1) + len_n1 + len(p2) + len(D1) + len(p3) + len(n2)],
                  "P-REGION J 5'": [len(V1) + len(p1) + len_n1 + len(p2) + len(D1) + len(p3) + len(n2) + 1, len(V1) + len(p1) + len_n1 + len(p2) + len(D1) + len(p3) + len(n2) + len(p4)],
                  "J-REGION": [len(V1) + len(p1) + len_n1 + len(p2) + len(D1) + len(p3) + len(n2) + len(p4), len(V1) + len(p1) + len_n1 + len(p2) + len(D1) + len(p3) + len(n2) + len(p4) + len(J1)]}
    df = pd.DataFrame.from_dict(annotation, orient = 'index') #это всё ради того, чтобы знать исходную разметку.
    df.to_csv("./Annotations/" + str(n) + "ann.csv", encoding='utf-8', index=True)
    if stop == 1:
        V_dom = V_dom.replace("tag", "tac") #тк слишком часто образуются стоп-кодоны, то мы грубо избавимся от них.
        V_dom = V_dom.replace("taa", "tac")
        V_dom = V_dom.replace("tga", "tgc")
    var = SeqRecord(Seq(V_dom), id = str(n) + "seq")
    output_handle = open(sys.argv[2] +".fasta", "a")
    SeqIO.write(var, output_handle, "fasta"),
    output_handle.close()
    n+=1
