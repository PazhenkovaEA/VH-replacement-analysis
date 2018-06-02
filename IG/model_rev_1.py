
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import re
from matplotlib import pyplot as plt
import os

os.chdir("/Users/irina/Desktop/project/1-162")
list_ids = [] #list of ids from Partis output, the first id from each clonal family
for n in range(1,162):
    try:
        aa= pd.read_csv(str(n) + '-cluster-annotations.csv')
        for _,row in aa.iterrows():
            list_ids.append(row["unique_ids"].split(":")[0])
    except:
        continue
        
os.chdir("/Users/irina/PycharmProjects/bioinf_2sem/small_phenotypes")
def cut_num(paper_name):
    x = re.search('[0-9]+\)\s', paper_name)    
    start_pos = x.end()
    return paper_name[start_pos:]

def make_ID_dict(fname):
    ID_dict = {}
    counter = 0
    with open(fname,'r') as f:
        for line in f:
            counter += 1
            if counter%3 == 1:
                ID = line.strip()
            elif counter%3 == 2:
                paper = line.strip()
                ID_dict[ID] = paper
    return ID_dict
    
def pnpe(n, pp, pn, pe):
    if n > 0:
        prob = ((pe**2*pn*pp**2*((n-1)*(pe*(1-pp))**(n-2)))/((pn-pp)*(pe+pp-pe*pp)**(n+2))) + ((pe*pn*pp**2*(pe*(1-pn))**n)/((pn-pp)**2*(pe-pe*pn)*(pe+pn-pe*pn)**(n+1))) - ((pe*pn*pp**2 * (pe*(1-pp))**n * (pe-2*pn + 3*pp + 2* pe*pn - 3 *pe*pp))/((pn - pp)**2*(pe-pe*pp)*(pe + pp - pe*pp)**(n+2)))
        return prob
    else:
        prob = ((pe**2*pn*pp**2*((1/(pe**2*(pp-1)**2))+(((n-1)*((-pe*(pp-1))/(pe+pp-pe*pp))**(n-2))/((pe+pp-pe*pp)**2))))/((pn-pp)*(pe-pp-pe*pp)**2)) - ((pn*pp**2*(pe-1)**3)/((pe+pn - pe*pn)*(pe+pp - pe*pp)**2)) +((pe*pn*pp**2*(((-pe*(pn-1))/(pe+pn-pe*pn))**n-1))/((pn-pp)**2*(pe-pe*pn)*(pe+pn-pe*pn))) - ((pe*pn*pp**2*(((-pe*(pp-1))/(pe+pp-pe*pp))**n-1)*(pe-2*pn+3*pp+2*pe*pn - 3*pe*pp))/((pn - pp)**2*(pe-pe*pp)*(pe+pp-pe*pp)**2))
        return prob

def likelihood(params):
    pp, pn, pe = params[0], params[1], params[2]
    return -np.sum([(y_exp_norm[x0])*np.log(pnpe(x0, pp, pn, pe)) for x0 in x_exp])
#split up dataset on phenotypes
handle = pd.ExcelFile('table_phenotypes.xlsx')
data = handle.parse(sheetname="Sheet1",skiprows=1)
phenotypes = set(data.columns)
phen_exclude = set(['Unnamed: 23', '//T helper cells', 'Unnamed: 31', '//SHM and Absence of AID','//Primary cold agglutinin disease','//Anti-platelets Ab','//Human Ig combinatorial library from genomic V segments and synthetic CDR3 fragments','//In vitro assembly and Sterile DJH rearrangements','//Women with myasthenia: evidence for immunization by fetal antigen','//Xeroderma pigmentosum','//IgD-only B cells','//T helper cells'])
phenotypes.difference_update(phen_exclude)

paper_dict = {cut_num(paper):phen for phen in phenotypes for paper in data[phen] if (type(paper)==str and ')' in paper)}
ID_dict = make_ID_dict('id_title-productive_seqs.txt')
    
ID_phen_dict = {ID:paper_dict[paper] for ID,paper in ID_dict.items() if paper in paper_dict}
    
len_data = pd.read_csv('Nt_seq.csv') #V-quest output
    
phen_hist_dict = {phenotype:np.zeros((1000)) for phenotype in phenotypes}
for _,row in len_data.iterrows():
    current_ID = row['Sequence ID']
    if current_ID in list_ids: #if we don't need filter sequences from one clonal family, just delete or comment this line
        if current_ID in ID_phen_dict:
            phen_hist_dict[ID_phen_dict[current_ID]][int(row['lens_n1'])] += 1

        
        
sortednames=sorted(phen_hist_dict.keys(), key=lambda x:x.lower()) 
phen_hist_dict = {sor:phen_hist_dict[sor] for sor in sortednames}
parametrs = []
all_data = []
#find the parameters for each phenotype
for phen in phen_hist_dict:
    if sum(phen_hist_dict[phen]) >= 80: #only for more than 80 samples
        x_exp = list(np.arange(50))
        y_exp = list(phen_hist_dict[phen][:50])
        y_exp_norm = [y / sum(y_exp) for y in y_exp]
        initial_guess = np.array([0.51, 0.49, 0.5]) #initial parameters
        result = minimize(likelihood, x0=initial_guess, bounds=((0.2, 0.999), (0.01, 0.999), (0.001, 0.999))) #minimization of ML function
        coeff = result.x
        y = [pnpe(z, coeff[0], coeff[1], coeff[2]) for z in x_exp]
        parametrs.append([coeff, phen])
        all_data.append([x_exp, y_exp_norm, y, phen])
    else:
        continue
dd = []
ss = []
for i in range(len(all_data)):
    ss = []
    for j in range(3):
        ss.append(parametrs[i][0][j])
    dd.append([ss + [str(parametrs[i][1])]])

data = [tuple(i[0]) for i in dd]
labels = ['pp', 'pn', 'pe', 'phenotype']
par_df = pd.DataFrame.from_records(parametrs)
par_df.to_csv("Phenotypes_no_clones.csv") #Save parameters for each phenotype to csv.

fig = plt.figure(figsize=(25, 25)) #visualization
columns = 4
rows = 4
for i in range(1, 17):
    fig.add_subplot(rows, columns, i)
    plt.plot(all_data[i - 1][0], all_data[i - 1][1], color="red")
    plt.plot(all_data[i - 1][0], all_data[i - 1][2], color="blue")
    plt.xlabel("Length of N1")
    plt.ylabel("Prob")
    plt.title(str(all_data[i - 1][3])[2:30])
plt.savefig('pheno.png', bbox_inches='tight')

