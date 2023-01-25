# NPOmix v1.0

## Publication:

![logo](https://github.com/tiagolbiotech/NPOmix/blob/main/Screen_Shot_2021-08-12_at_7.18.10_PM.png)

```
NPOmix: A machine learning classifier to connect mass spectrometry fragmentation data to biosynthetic gene clusters 

Tiago F Leão,  Mingxun Wang,  Ricardo da Silva,  Alexey Gurevich,  Anelize Bauermeister, Paulo Wender P Gomes,  
Asker Brejnrod,  Evgenia Glukhov,  Allegra T Aron, Joris J R Louwen,  Hyun Woo Kim,  Raphael Reher,  Marli F Fiore, 
Justin J J van der Hooft,  Lena Gerwick,  William H Gerwick,  Nuno Bandeira, Pieter C Dorrestein

PNAS Nexus, Volume 1, Issue 5, November 2022, pgac257.
DOI: https://doi.org/10.1093/pnasnexus/pgac257
```
[NPOmix PNAS Nexus link](https://academic.oup.com/pnasnexus/article/1/5/pgac257/6847575)

## More information about NPOmix at (including workshops): 

https://www.tfleao.com/npomix1

## Quick tutorial:

### Instructions:

**1) Clone the GitHub repository:**
```
cd path/to/root/folder
git clone https://github.com/tiagolbiotech/NPOmix_python.git
```

**2) Download the training BGCs (ideally in the `test/` folder):** 
```
wget -O antismash_bgcs.zip https://zenodo.org/record/6637083/files/antismash_only_gbk.zip?download=1
unzip antismash_bgcs.zip && rm antismash_bgcs.zip
```

**3) Clone the GNPS job with the metabolomes for the strains in the training set and add your files (ideally to G3 group):**
https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=cc4bd4de2f184baf9bf6defef4492927

**4) Download GNPS output. Run and download the antiSMASH file. Place the antiSMASH file in the `antismash_gbk_only` folder.**

**5) Run BiG-SCAPE for the antiSMASH BGCs (including from your samples) and concatenate outputs for all classes.** 

**6) Run the dereplication or not dereplication notebooks. Make sure to adjust the paths for the antiSMASH, GNPS and BiG-SCAPE files.**

----------------------------------------------------------
For detailed instructions, please check the video below:

https://user-images.githubusercontent.com/12102722/214508959-6eb80a78-a63b-41cb-a8f1-81452f108dc5.mp4

## Submit your samples

If you have difficulties running the NPOmix tool for your samples or have questions, please submit your samples at the link below.

https://www.tfleao.com/general-8

## Overview of the methodology

To use the NPOmix approach (Fig. 1, schematic example for the approach used in only four samples), you will need a dataset of so-called paired genome-MS/MS samples, samples that contain both a genome (or metagenome) and a group of MS/MS spectra obtained via untargeted LC-MS/MS. Many paired datasets are available at the Paired omics Data Platform [(PoDP)](https://pairedomicsdata.bioinformatics.nl), one of the first initiatives to gather paired genome-MS/MS samples. By applying [BiG-SCAPE](https://bigscape-corason.secondarymetabolites.org), each biosynthetic gene cluster (BGC) in the genomes will go through a pairwise similarity comparison (Fig. 1A) to every other BGC in the same set of genomes to compute similarity scores (1 minus BiG-SCAPE raw distance) and to assign BGCs to Gene Cluster Families (GCFs), if possible. In order to create a BGC fingerprint (Fig. 1C), we identify the maximum similarity of the query BGC to one of the many BGCs in each genome (which can be considered a pool of BGCs) in the dataset. Therefore, the BGC fingerprints can be represented as a row of values (a vector with the maximum similarity scores), and each column is a different genome from the selected dataset. Similarity scores range from 0.0 to 1.0; for instance, identical BGCs have a perfect similarity score of 1.0, a score of 0.8 would represent that a homologous BGC is present in the genome, and the score of 0.0 (or below the similarity cutoff of 0.7) represents that the queried BGC is absent in the genome. A similar process happens to create MS/MS fingerprints (Fig. 1B); however, genomes are replaced by MS/MS spectra, and a query BGC is replaced by a query MS/MS spectrum; either a reference spectrum from [GNPS](https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp) or a cryptic MS/MS spectrum from a new sample (that contains a genome and experimental MS/MS spectra). In the case of MS/MS fingerprints (Fig. 1D), we used GNPS to calculate the pairwise modified cosine score and then identify the maximum similarity of the query MS/MS spectrum to one of the many MS/MS spectra in each experimental sample. Of particular note, we did not use the full classical GNPS molecular networking capabilities; we only used the functions required to calculate a modified cosine similarity score between a pair of MS/MS spectra. The BGC fingerprints create a training matrix (Fig. 1E) where rows are the maximum similarity scores, normally thousands of rows (e.g., for our first release, round 4, we have used 5,421 BGCs in 1,040 genomes/metagenomes), where each genome is a column. This matrix can be fed to the k-nearest neighbor (KNN) algorithm in order to train it with the genomic data. Additionally, one extra column is required in this genomic data matrix, a column that labels each BGC fingerprint with a GCF so the KNN algorithm will know the fingerprint patterns that belong together. The KNN algorithm plots the BGC fingerprints in the KNN space (in Fig. 1G, the KNN space is exemplified by only 2 dimensions). Next, the MS/MS fingerprints also form a testing matrix (Fig. 1F); in this case, this matrix also contains 1,040 columns due to the 1,040 sets of paired experimental MS/MS spectra. For example, for our first release, this testing matrix contained 15 reference MS/MS fingerprints (rows) for MS/MS reference spectra from the GNPS database. Each query MS/MS fingerprint (a row in the testing metabolomic matrix and columns are the experimental MS/MS spectra per sample) will also be plotted into the same KNN space (Fig. 1G) so the algorithm can obtain the GCF labels for the k-nearest neighbors to the query MS/MS fingerprint (e.g., for three most similar BGC neighbors, k = 3). We note that GCF labels can be present more than once in the returned list if two or more BGC nearest neighbors belong to the same gene family. This repetition of the GCF classification is a common behavior of the KNN approach. Our approach is suitable for bacterial, fungal, algal, and plant genomes and MS/MS spectra. Metagenomes and metagenome-assembled genomes can also be used instead of genomes; however, complete genomes are preferred. This KNN approach also supports LC-MS/MS from fractions or from different culture conditions; multiple LC-MS/MS files for the same genome were merged together into a single set of experimental MS/MS spectra.

![Fig1_part1](https://github.com/tiagolbiotech/NPOmix/blob/main/Screen%20Shot%202021-06-23%20at%201.35.17%20PM.png)

![Fig1_part2](https://github.com/tiagolbiotech/NPOmix/blob/main/Screen%20Shot%202021-06-23%20at%201.35.53%20PM.png)

## Video overview:

https://user-images.githubusercontent.com/12102722/131203358-9576574d-4526-414a-b6a5-b94125ea8b18.mp4

### References:
1) 3% of the biosynthetic potential: Gavriilidou, A., Kautsar, S. A., Zaburannyi, N., Krug, D, Muller, R., Medema, M. H. & Ziemert, N. A global survey of specialized metabolic diversity encoded in bacterial genomes. bioRxiv (2021).
2) Piared omics Data Platform: Schorn, M. A. et al. A community resource for paired genomic and metabolomic data mining. Nat. Chem. Biol. 17, 363–368 (2021).
3) New approved drugs from 1981 to 2014: Newman, D. J. & Cragg, G. M. Natural Products as Sources of New Drugs from 1981 to 2014. J. Nat. Prod. 79, 629–61 (2016).

## Video details on the methodology:

https://user-images.githubusercontent.com/12102722/133908907-b1ee838c-1180-40bb-aa27-b50bbecbd329.mp4
