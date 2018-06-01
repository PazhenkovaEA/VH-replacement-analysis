# Model of the N1 zone formation in human antibodies’ gene sequences

## Project description

  The N1-zone is a variable region of human antibodies DNA, formed as a result of VDJ-recombination and providing diversity of antigen binding regions. N1-zone generation is a complicated process including formation of palindroms on 5' and 3' ends and addind random nucleotides to 5' and 3' ends with following non-homologous end joining. The most poorly understood factor is exonuclease trimming, which also contributes to the BCR and TCR repertoires (Jackson et al. 2013). Consequently, we made an assumption, that the length of N1-zone depends of several random events. According to this approval, we created a statistical model of N1-zone formation and estimating its parameters for synthetical and experimental datasets using Maximum Likelihood method.

## Students

Elena Pazhenkova, Darya Krytskaya 
## Supervisors 

Evgeniy Bakin, Oksana Stanevich

## Goal and tasks
Create a statistical model for the formation of the N1 zone of human IgH.

Tasks:
* to study details of sequence analysis in [IMGT/V-quest](https://www.imgt.org/IMGT_vquest/vquest)
* to define biological stages of antibody formation
* to propose mathematical models, describing N1-region formation
* to estimate parameters of the model for synthetic datasets
* to find the best model for experimental data and evaluate its biological relevance

## Material and methods
### The data

The data was collected from [NCBI](https://www.ncbi.nlm.nih.gov/nucleotide/). 
The dataset contained 32301 sequences, than it was splitted up on different phenotypes. We divided our dataset on clonal families using [Partis](https://github.com/psathyrella/partis/blob/master/manual.md) and filtered clones to avoid presence of non-random samples. After filtering, we got following phenotypes:
* Healthy people: N = 2742
* Systemic lupus erythematosus: N = 1257
* Multiple sclerosis: N = 202
* X-linked hyper-IgM syndrome: N = 994
* Wegners granuloma: N = 160
* Chronic lymphocytic leukemia: N = 741
* Hodgkin lymphoma: N = 98
* Non-Hodgkin lymphomas: N = 328
* Hepatitis C: N = 317
* HIV-1: N = 80
* Infectious mononucleosis: N = 130
* Acute viral and bacterial infections: N = 249
* Pneumococcal vaccine: N = 97 
Parameters of the model were estimated with Maximum likelihood method and real and modelled distriburions were visualized before and after filtering of dependencies in Partis. 

## Results

### The model
Under general independency assumption it follows that:
* palindrome length is a geometrically distributed random variable:
<a href="https://www.codecogs.com/eqnedit.php?latex=Lp&space;\sim&space;Geom(pp)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Lp&space;\sim&space;Geom(pp)" title="Lp \sim Geom(pp)" /></a>

* N-nucleotides length is a geometrically distributed random variable:
<a href="https://www.codecogs.com/eqnedit.php?latex=Ln&space;\sim&space;Geom(pn)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Ln&space;\sim&space;Geom(pn)" title="Ln \sim Geom(pn)" /></a>

* Resulting N1-zone length is a binomial random variable:
<a href="https://www.codecogs.com/eqnedit.php?latex=L&space;\sim&space;Bin(Lp1&plus;Lp2&plus;Ln,&space;pe)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L&space;\sim&space;Bin(Lp1&plus;Lp2&plus;Ln,&space;pe)" title="L \sim Bin(Lp1+Lp2+Ln, pe)" /></a>

The model has three parameters: 
<a href="https://www.codecogs.com/eqnedit.php?latex=pp,&space;pn,&space;pe" target="_blank"><img src="https://latex.codecogs.com/gif.latex?pp,&space;pn,&space;pe" title="pp, pn, pe" /></a>

The final formula of model is ![model](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/Figures/distribution.png)

### Synthetic dataset
  
  We generate synthetic dataset with our script [genejunc.py](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/IG/genjunc.py). The idea of the script based on biological stages of antibody formation in human according to Murphy & Casey, 2017. Briefly, the random V, D, J regions are chosen from germline dataset (we use productive genes only). Palindrome length depends on position, where Artemis/DNA-PK complex make a single-strand brake, which is determined by depth of Ku70-Ku80 touchdown (Figure 1). 
 
 **Figure 1**
  ![fig1](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/Figures/Bez_imeni-1.png)
  
  Thus, we can expect, than palindrome length is a geometrically distributed random variable. Palindromes are added on 3' of V-gene, 3' and 5' of D-gene and 5' of J-gene, and they are complimentary to the first (in case of 5'end) or last (in case of 3'end) several nucleotides of corresponding region. Thus, geometrically distributed random size palindromes are formed on the borders of all three genes. Between palindroms random number of random n-nucleotides is added. We assume, that N-nucleotides length is a geometrically distributed random variable. Exonuclease trimming evaluated for each nucleotides separately with Bernoulli distribution probability, which is a simplification. To produce only productive sequences (optional), the last symbol in stop-codons is changed to 'C' and 1-2 nucleotides are added in the N2-region to avoid frameshifts. 
  
  The script can be runned from terminal with following command:
  
  ```
  python3 genjunc.py <seq_number> <output> <stop codons> <frameshifts> 
  ```
  
* seq_number - amount of generated sequences
* output - output fasta name w/o extention
* stop codons - remove stop codons if 1 (default 1)
* frameshifts - treat frameshifts if 1 (default 1)


For example, following command returns "test.fasta" file, contains 50 productive (without frameshifts and stop-codons inside) sequences:
```
python3 genjunc.py 50 test
```
And this command returns "test1.fasta" file contains 100 sequences without stop-codons, but with possible frameshifts:
```
python3 genjunc.py 100 test 1 0
```
The file "Annotations_filename.csv" will be created. The table contains the informations about borders of each region.
Also the table "length_filename.csv" with N1-zone length for each sequence will be created.

Files V.fasta, D.fasta and J.fasta contain germline alleles for human IgH from [IMGT-GENE](http://www.imgt.org/genedb/) database and they should be placed in the working directory.

**Prerequisites**: Python 3.x, Biopython (Cock et al. 2009) and pandas library are required. 

### IMGT/V-quest results

We submited synthetic sequences to IMGT V-quest and compared result of analysis to our annotation. On Fig. 2 an example of the whole V-D-J region is represented. Borders of frameworks (Fr1-4) and complementarity-determining regions (CDR1-3) are marked by IMGT.  

  **Figure 2**
![Bigpicture](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/Figures/Fig1.png)

On Fig. 3 the Junction (from 3'V to 5'J) of the same sequence as on Fig. 2 is shown. The main difference between our and V-quest' annotation is that V-quest often includes 3'P-nucleotides in N-region. Considering this, the junction between V and D will be reffered to as N1-zone, and the junction between D and J will be referred to as N2-zone. This terminology also accepted in Meng et al., 2014. 


   **Figure 3**
![smallpicture](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/Figures/Fig2.png)

The length of N1 zones were equal in 99% of 150000 generated sequences according to IMGT and initial annotation (Fig. 4).
  **Figure 4**

![heatmap](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/Figures/heatmap.png)

Initial parameters chosen for synthetic dataset generation and parameters of model, estimated with Maximum likelihood method before and after annotation in V-quest differed not more than 10%. Thus, both model and V-quest are suitable for using in real datasets. 

|   | pp | pn | pe |
|--|--|--|--|
| Initial | 0.5 | 0.1 | 0.9 |
| Before V-quest | 0.5 | 0.1 | 0.9 |
| After V-quest | 0.43 | 0.09 | 0.8 |

### Testing model on experimental datasets

The model was tested on data, collected from NCBI and splitted up on phenotypes. 
For each phenotypes parameters were estimated with Maximum likelyhood method and distributions were plotted (red - experimental data, blue - model distribution)
![Figure 5](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/Figures/1_graph.png)

In cases of healthy phenotypes and some non-healthy (chronic lymphomic leukemia, infectious mononucleosis), model well describes real data N1-length distributions. However, in HIV-1 infection, non-Hodgkin lymohomas, systemic lupus erythematos, multiple sclerosis (and etc) phenotypes some lengths of N1 represented by larger number of copies, than other. This can be caused by presence of clones in our samples, which means, that our N1 lenght sampling is not random. To avoid such limitations, we filtered clones, detected by Partis, and fit the model with filtered data.

![Figure 6](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/Figures/2_graph.png)

However, plots for some unhealthy phenotypes still shows discrepancy in model and real N1 length distributions. This can be explained by presence of additional factors playing role in N1-zone formation. The N1-zone sometimes contains so-called footprints, appeared as a result of VH-replacement and recent studies showed that the length of CDR3 (including V3', N1, D, N2 and J5') is correlated with number of footprints (Meng et al., 2014). Futher investigations are nessesary to test this hypothesis. 


[The script](https://github.com/PazhenkovaEA/VH-replacement-analysis/blob/master/IG/model_rev_1.py) was used to choose the parameters for the model and to visualize the resulted distribution for each phenotype.

**Prerequisites**: Python 3.x, libraries matplotlib, pandas, numpy, scipy

## Further plans
* Analyze model adequacy via Pearson test
* Perform experiments with high-quality datasets
* Loosen independence assumptions in the model


## Citation
* Brochet, X., Lefranc, M. P., & Giudicelli, V. (2008). IMGT/V-QUEST: the highly customized and integrated system for IG and TR standardized V-J and V-D-J sequence analysis. Nucleic Acids Research, 36(Web Server issue). https://doi.org/10.1093/nar/gkn316
* Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … De Hoon, M. J. L. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
* Meng, W., Jayaraman, S., Zhang, B., Schwartz, G. W., Daber, R. D., Hershberg, U., … Luning Prak, E. T. (2014). Trials and tribulations with VH replacement. Frontiers in Immunology, 5(JAN). https://doi.org/10.3389/fimmu.2014.00010
* Murphy, K., Casey, W. (2017). Janeway's immunobiology. 9th edition. New York: Garland Science, London: Taylor & Francis Group
* Partis https://github.com/psathyrella/partis/blob/master/manual.md


P.S. “Essentially, all models are wrong, but some are useful.” George Box
