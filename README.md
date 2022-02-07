# Litertaure analysis of epilepsy-associated genes

This repository contains scripts used for the generation of Figure 1 and 2 from the manuscript "From gene discovery in human epilepsies to variant interpretation", in review.

## 1) The history of epilepsy-associated genes and number of references since discovery

#### Script: epilepsy_gene_year_analysis.R

A) All NCBI PubMed abstracts that mention ‘epilepsy’ along at least one clinical and one genetic term were downloaded and processed (12/05/21).

```
PubMed query: "epilepsy AND (\"patient\" OR \"patients\" OR \"family\" OR \"proband\" OR \"case\" OR \"cases\") AND (\"mutation\"  OR \"deletion\" OR \"duplication\" OR \"amplification\" OR \"variant\" OR \"CNV\" OR \"fusion\" OR \"microdeletion\" OR \"missense\" OR \"triplication\" OR \"pathogenic variant\" OR \"truncation\" OR \"stop\" OR \"frameshift\") NOT systematic review[Filter] NOT review[Filter]"
```

B) From those abstracts all genes, along with the first year of discovery and the overall count of references, were extracted using PubTator.

C) Manual quality control of all genes that were on average mentioned at least once a year after the initial discovery or that were defined as epilepsy-associated genes by Heyne et al., Lindy et al, or ClinGen (version: 01/10/22).

D) Plotting of timeline of all remaining genes and their year of first reference using ggplot2. (X-axis= year, Y-axis= number of genes)

![Figure 1](/data_gene_year_analysis/screenshot_fig1.png)

## B) Epilepsy-associated genes with Loss-of-Function and/or Gain-of-Function effect variants

#### Script: epilepsy_gof_lof_analysis.R

A) All  NCBI PubMed abstracts that mention ‘epilepsy’ along at least one Gain-of-function (GoF) or one Loss-of-Function (LoF) term were downloaded and processed (12/05/21).

```
PubMed query: epilepsy AND ("GoF" OR "gain of function" OR "gain-of-function" OR "LoF" OR "loss of function" OR "loss-of-function") NOT systematic review[Filter] NOT review[Filter]
```

B) From those abstracts all genes were extracted using PubTator.

C) Abstracts of the 50 most frequent genes were annotated with the information whether the abstract is mentioning a GoF term and/or a LoF term and manually revised. 

D) Plotting (stacked barplot) of all genes with at least five abstracts and their association with LoF or GoF variants, or variants of both categories.

![Figure 2](/data_gof_lof_analysis/screenshot_fig2.png)

## References

1. 	Lindy AS, Stosser MB, Butler E, et al. Diagnostic outcomes for genetic testing of 70 genes in 8565 patients with epilepsy and neurodevelopmental disorders. Epilepsia. 2018;59(5):1062-1071. doi:10.1111/epi.14074

2. 	Rehm HL, Berg JS, Brooks LD, et al. ClinGen — The Clinical Genome Resource. http://dx.doi.org/10.1056/NEJMsr1406261. doi:10.1056/NEJMsr1406261

3. 	Heyne HO, Singh T, Stamberger H, et al. De novo variants in neurodevelopmental disorders with epilepsy. Nature Genetics. 2018;50(7):1048-1053. doi:10.1038/s41588-018-0143-7






