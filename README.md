# Analysis of VH-replacement in human antibodies

## Students

Darya Krytskaya, Elena Pazhenkova

## Supervisors 

Evgeniy Bakin, Oksana Stanevich

## Goal and tasks
Create a statistical model for the formation of the N1 zone of human IgH.

Tasks:
* to study details of sequence analysis in [IMGT-Vquest](https://www.imgt.org/IMGT_vquest/vquest)
* to define biological stages of antibody formation
* to propose mathematical models, describing N1-region formation
* to estimate parameters of the model for synthetic and real datasets
* to find the best model and evaluate its biological relevance

## Methods
### Test data

The data was collected from [NCBI](https://www.ncbi.nlm.nih.gov/nucleotide/). 
The dataset contains 32000 sequences, among them 18000 are productive.
The dataset will be divided by source of Ig (bone marrow or periphery) and by phenotype (healthy or not).

### Synthetic dataset
  
  We generate synthetic dataset with our script [genejunc.py](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/IG/genjunc.py). The idea of the script based on biological stages of antibody formation in human according to Murphy & Casey, 2017. Briefly, the random V, D, J regions are chosen from germline dataset (we use productive genes only). On 3' of V-gene, 3' and 5' of D-gene and 5' of J-gene added 1-4 nucleotides, which are complimentary to 1-4 first (in case of 5'end) or last (in case of 3'end) nucleotides of corresponding region. Thus, random size (2-8) palindroms are formed on the borders of all three genes. To produce only productive sequences (optional), the last symbol in stop-codons is changed to 'C' and 1-2 nucleotides are added in N1-region to avoid frameshifts. Between palindroms random number (up to 20) of random n-nucleotides is added. Such biological processes as nucleases activity and VH-replacement will be modelled.  
  
  The script can be runned from terminal with following command:
  
  ```
  python3 genjunc.py <seq_number> <output> <stop codons> <frameshifts> <VH-replacement> <nucleases>
  ```
  
* seq_number - amount of generated sequences
* output - output fasta name w/o extention
* stop codons - remove stop codons if 1 (default 1)
* frameshifts - treat frameshifts if 1 (default 1)
* VH-replacement - imitate VH-replacement if 1 (default 0) # under development
* nucleases - remove random nucleotides if 1 (default 0)# under development

For example, following command returns "test.fasta" file, contains 50 productive (without frameshifts and stop-codons inside) sequences:
```
python3 genjunc.py 50 test
```
Also, the folder "Annotations" with .csv file for each sequence is created. The tables contains the informations about borders of each region.

And this command returns "test1.fasta" file contains 100 sequences without stop-codons, but with possible frameshifts:
```
python3 genjunc.py 100 test 1 0
```

Files V.fasta, D.fasta and J.fasta contain germline alleles for human IgH from [IMGT-GENE] (http://www.imgt.org/genedb/) database and they should be placed in the working directory.

Python 3.x, Biopython (Cock et al. 2009) and pandas library are required. 

### IMGT-Vquest results

картинку сделаю

## Citation
* Brochet, X., Lefranc, M. P., & Giudicelli, V. (2008). IMGT/V-QUEST: the highly customized and integrated system for IG and TR standardized V-J and V-D-J sequence analysis. Nucleic Acids Research, 36(Web Server issue). https://doi.org/10.1093/nar/gkn316
* Brochet, X., Lefranc, M. P., & Giudicelli, V. (2008). IMGT/V-QUEST: the highly customized and integrated system for IG and TR standardized V-J and V-D-J sequence analysis. Nucleic Acids Research, 36(Web Server issue). https://doi.org/10.1093/nar/gkn316
* Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … De Hoon, M. J. L. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
* Murphy, K., Casey, W. (2017). Janeway's immunobiology. 9th edition. New York: Garland Science, London: Taylor & Francis Group
