### script to identify gof/lof genes in the literature

# sections:
# - STEP 1: Downlaod all abstracts (csv. and .txt format) from PubMed mentioning epilepsy along one gof / lof term
# - STEP 2: extract all genes from the abstracts using PubTator
# - STEP 3: Generate list of most frequent genes and their abstracts and annotate whether the genes is mentioned with gof or lof term
# - STEP 4: Manual control of lof/gof status of all abstracts from the top50 genes 
# - STEP 5: Barplot of genes with >5 abstracts gof/lof terms


library(readr)
library(stringr)
library(stringi)
library(RCurl)
library(easyPubMed)
library(tidyr)
library(ggrepel)
library(ggpubr)
library(readxl)
library(tidyselect)
library(tidyr)
library(dplyr)

###load data
gene_df <- read_delim("data_gof_lof_analysis/Homo_sapiens.gene_info","\t", escape_double = FALSE, trim_ws = TRUE)

#### STEP 1: Downlaod all abstracts (csv. and .txt format) from PubMed mentioning epilepsy along one gof / lof term ####

#search query used on PubMed website to download csv and PubMed file of all entries: 
#epilepsy AND ("GoF" OR "gain of function" OR "gain-of-function" OR "LoF" OR "loss of function" OR "loss-of-function") NOT systematic review[Filter] NOT review[Filter]

#load csv file from pubmed
PM_df <- read_csv("gof_lof_analysis/csv-epilepsyAN-set.csv")

#load abstracts from pubmed file to annotate them to the pmid-csvfile
pubmed_data = read_lines(
  "gof_lof_analysis/pubmed-epilepsyAN-set.txt",
  skip = 0,
  skip_empty_rows = FALSE,
  n_max = Inf
)

pubmed_data = tibble(pubmed_data)

#tag all abstract lines
for(i in grep("^(AB  -)", pubmed_data$pubmed_data)){
  cat(i)
  n=1
  while(grepl("^(   )",pubmed_data[i+n,1])){
    pubmed_data[i+n,1] = paste0("AB  -", gsub("      "," ",pubmed_data[i+n,1]))
    n = n+1
  }
}

pubmed_data = pubmed_data[grepl("^(AB  -|PMID-)", pubmed_data$pubmed_data),]
saveRDS(pubmed_data, "gof_lof_analysis/pubmed_abstract_data")

#extract abstract of each PMID and add it to pubmed dataframe
for(i in 1:nrow(PM_df)) {
  
  abstract_index = grep(paste0("PMID- ",PM_df$PMID[i]), pubmed_data$pubmed_data)
  
  n=1
  tmp_abstract = c()
  
  if(grepl("^(AB  -)",pubmed_data[abstract_index+n,1])) {
    while(grepl("^(AB  -)",pubmed_data[abstract_index+n,1])){
      tmp_abstract = paste0(tmp_abstract, gsub("AB  -","",pubmed_data[abstract_index+n,1]))
      n = n+1
    }
    PM_df[i, "abstract"] = tmp_abstract
  }
  
}
saveRDS(PM_df, "gof_lof_analysis/pubmed_csv_with_abstract_data")
PM_df <- readRDS("data_gof_lof_analysis/pubmed_csv_with_abstract_data")

#### STEP 2: extract all genes from the abstracts using PubTator ####

#extract genes from PMIDs using PubTator
pmids <- PM_df$PMID
n <- 0
gene_table <- data.frame()

for (pmid_split in split(pmids, ceiling(seq_along(pmids)/99))){
  out.data <- NULL
  n <- n+1
  
  out.data <- tryCatch({    
    getURL(paste("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=", 
                 paste(pmid_split, collapse=","), sep = ""))
  }, error = function(e) {
    print(e)
    next
  }, finally = {
    Sys.sleep(5)
  })
  
  if(!is.null(out.data)){
    
    cat(n)
    
    out.data <- unlist(strsplit(out.data, "\n", fixed = T))
    
    for (i in 3:length(out.data)) {
      temps <- unlist(strsplit(out.data[i], "\t", fixed = T))
      
      if (length(temps) == 5) {
        temps <- c(temps, NA)}
      
      if (length(temps) == 6) {
        gene_table <- rbind(gene_table, temps)
        gene_table <- gene_table[as.character(gene_table[,5]) == "Gene",]
        gene_table <- gene_table[!duplicated(gene_table[,c(1,6)]),] #count genes once per abstract
      }
    }
  }}

colnames= gsub("X.", "", names(gene_table))
colnames= gsub("\\.*$", "", colnames)
gene_table <- rbind(gene_table, colnames) #column names were first row of data
gene_table <- gene_table[,c(1,4,6)]
names(gene_table) <- c("pmid","genename","geneID")

gene_table <- separate_rows(gene_table, geneID, sep=";") #when there are several geneIDs split row in two (eg. CDKL5A/B -> two rows as two different IDs)

saveRDS(gene_table, "gof_lof_analysis/genetable_epilepsy")


#### STEP 3: Generate list of 50 most frequent genes and their abstracts and annotate whether the genes is mentioned with gof or lof term ####

###get counts of genes (table of frequencies of geneIDs)
genecounts_df <- as.data.frame(table(gene_table[,3])) 
genecounts_df <- data.frame(geneID = as.character(genecounts_df[,1]),freq = as.numeric(genecounts_df[,2]))
genecounts_df <- genecounts_df[order(genecounts_df$freq,decreasing=T),]

###annotate abstract information to genes
gene_table <- merge(gene_table, PM_df, by.x = "pmid", by.y = "PMID")
gene_table <- gene_table %>%
  separate(`Create Date`, c("year", "month", "day"), "-")

saveRDS(gene_table, "data_gof_lof_analysis/genetable_epilepsy_with_annotated_abstracts")
gene_table = readRDS("data_gof_lof_analysis/genetable_epilepsy_with_annotated_abstracts")

PM_df_gof_lof = PM_df
LOFGOF_genes = merge(PM_df_gof_lof, gene_table, by.x="PMID", by.y="pmid", all.x = T)

LOFGOF_genes = LOFGOF_genes[!is.na(match(LOFGOF_genes$geneID, gene_df$GeneID)),] 
LOFGOF_genes <- merge(LOFGOF_genes, gene_df, by.x = "geneID", by.y = "GeneID")
top50_lof_gof_genes = names(sort(table(LOFGOF_genes$geneID), decreasing=T)[1:50])
LOFGOF_genes = LOFGOF_genes[LOFGOF_genes$geneID  %in% top50_lof_gof_genes, ]

##only keep genes that are in the epilepsy gene -year plot
plot_data_geneplot <- readxl::read_excel("data_gene_year_analysis/epilepsy_genes_manual_checked.xlsx")
plot_data_geneplot = plot_data_geneplot[plot_data_geneplot$`EXCLUDE?` == "no", ]
LOFGOF_genes$gene_included = LOFGOF_genes$Symbol %in% plot_data_geneplot$Symbol
LOFGOF_genes = LOFGOF_genes[LOFGOF_genes$gene_included == T,]

#annotate whether gof or lof term found in abstract
LOFGOF_genes$gof = grepl("(GOF|gof|GoF|gain-of-function|gain of function|Gain-of-function|Gain of function|gain-of function|Gain-of function|gain of  function|Gain-of-Function|gain  of function|gain  of  function|gain of channel function)", LOFGOF_genes$abstract.x) | grepl("(GOF|gof|GoF|gain-of-function|gain of function|Gain-of-function|Gain of function|gain-of function|Gain-of function|gain of  function|Gain-of-Function|gain  of function|gain  of  function|gain of channel function)", LOFGOF_genes$Title.x)
LOFGOF_genes$lof = grepl("(LOF|lof|LoF|loss of function|loss-of-function|Loss-of-function|Loss of function|loss-of function|Loss-of function|loss of  function|Loss-of-Function|loss  of function|loss  of  function|loss of channel function)", LOFGOF_genes$abstract.x) | grepl("(LOF|lof|LoF|loss of function|loss-of-function|Loss-of-function|Loss of function|loss-of function|Loss-of function|loss of  function|Loss-of-Function|loss  of function|loss  of  function|loss of channel function)", LOFGOF_genes$Title.x)

#delete abstracts without gof/lof term
LOFGOF_genes=LOFGOF_genes[-which(LOFGOF_genes$gof + LOFGOF_genes$lof == 0),]

#### STEP 4: Manual control of all abstracts from the top50 genes ####

writexl::write_xlsx(LOFGOF_genes, "data_gof_lof_analysis/gof_lof_abstract_manual_control.xlsx")

#### STEP 5: Barplot of most frequent genes with gof/lof terms ####

##load excel file of manual quality control and plot data 
data_checked <- readxl::read_excel("data_gof_lof_analysis/gof_lof_abstract_manual_controlled.xlsx")
data_checked <- data_checked[data_checked$gof_verif == "T"|data_checked$gof_verif == "F", ]
data_checked$gof_verif <- ifelse(data_checked$gof_verif == "T", TRUE, FALSE)
data_checked$lof_verif <- ifelse(data_checked$lof_verif == "T", TRUE, FALSE)

plot_data = data.frame(Symbol= unique(data_checked$Symbol), GeneID = unique(data_checked$geneID), freq= as.numeric(table(data_checked$geneID)))

plot_data$lof = 0
plot_data$gof = 0
plot_data$'lof+gof' = 0

for(i in 1:nrow(plot_data)){
  
  tmp.data = data_checked[data_checked$geneID == plot_data$GeneID[i], ]
  tmp.data$'lof+gof_verif' = tmp.data$lof_verif + tmp.data$gof_verif
  
  for(r in 1:nrow(tmp.data)){
    if(tmp.data$'lof+gof_verif'[r] == 2){
      plot_data[i,'lof+gof'] = plot_data[i,'lof+gof'] + 1
    } else if (tmp.data$gof_verif[r]) {
      plot_data[i,"gof"] =  plot_data[i,"gof"] + 1
    } else{
      plot_data[i,"lof"] = plot_data[i,"lof"] + 1
    }
  }
  
}

plot_data$total = plot_data$lof+plot_data$gof+plot_data$'lof+gof'

plot_data = plot_data[plot_data$freq > 5,]

plot_data = gather(plot_data, key = "gof_lof", value = "count", gof, lof, 'lof+gof')

plot = ggplot(plot_data, aes(fill=gof_lof, y=count, x=reorder(Symbol, -count))) + 
  geom_bar(position="stack", stat="identity")+
  ylab("Number of references") +
  xlab("") +
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="bottom", legend.box = "horizontal", legend.title = element_blank())+
  scale_fill_manual(values = c("#52854C", "#4E84C4", "#293352"))


pdf("gof_lof_analysis/gof_lof_epilepsy_genes.pdf",width= 7,height= 4)
plot
dev.off()







