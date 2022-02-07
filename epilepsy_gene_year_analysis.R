### script to identify epilepsy-associated genes in the literature

# sections:
#   - STEP 1: Downlaod all abstracts (csv. and .txt format) from PubMed mentioning epilepsy, one clinical and one genetic term
#   - STEP 2: extract all genes from the abstracts using PubTator
#   - STEP 3: Generate list of most frequent genes and their first year of reference
#   - STEP 4: Manual control of epilepsy associated genes and their year of first reference
#   - STEP 5: Plot timeline of epilepsy-associated genes

# load packages
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
epilepsy_genes <- read_delim("data_gene_year_analysis/epilepsy_genes_10jan22","\t", escape_double = FALSE, trim_ws = TRUE)
gene_df <- read_delim("data_gene_year_analysis/Homo_sapiens.gene_info","\t", escape_double = FALSE, trim_ws = TRUE)

#### STEP1: Downlaod all abstracts (csv. and .txt format) from PubMed mentioning epilepsy and one clinical and one genetic term using the following query: ####
# "epilepsy AND (\"patient\" OR \"patients\" OR \"family\" OR \"proband\" OR \"case\" OR \"cases\") AND (\"mutation\"  OR \"deletion\" OR \"duplication\" OR \"amplification\" OR \"variant\" OR \"CNV\" OR \"fusion\" OR \"microdeletion\" OR \"missense\" OR \"triplication\" OR \"pathogenic variant\" OR \"truncation\" OR \"stop\" OR \"frameshift\") NOT systematic review[Filter] NOT review[Filter]"
# and process the data to a table with all PMIDs and their abstract text

#load csv file from pubmed
PM_df <- read_csv("data_gene_year_analysis/csv-epilepsyAN-set.csv")

#load abstracts from pubmed file to annotate them to the pmid-csvfile
pubmed_data = read_lines(
  "data_gene_year_analysis/pubmed-epilepsyAN-set.txt",
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
saveRDS(pubmed_data, "data_gene_year_analysis/pubmed_abstract_data")

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
saveRDS(PM_df, "data_gene_year_analysis/pubmed_csv_with_abstract_data")



#### STEP2: extract all genes from the abstracts using PubTator ####

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

saveRDS(gene_table, "data_gene_year_analysis/genetable_epilepsy")


#### STEP3: Generate list of most frequent genes and their first year of reference  ####

###get counts of genes (table of frequencies of geneIDs)
genecounts_df <- as.data.frame(table(gene_table[,3])) 
genecounts_df <- data.frame(geneID = as.character(genecounts_df[,1]),freq = as.numeric(genecounts_df[,2]))
genecounts_df <- genecounts_df[order(genecounts_df$freq,decreasing=T),]

###annotate abstract information to genes
gene_table <- merge(gene_table, PM_df, by.x = "pmid", by.y = "PMID")
gene_table <- gene_table %>%
  separate(`Create Date`, c("year", "month", "day"), "-")

saveRDS(gene_table, "data_gene_year_analysis/genetable_epilepsy_with_annotated_abstracts")

###only keep earliest entry of each gene
unique_gene_table <- data.frame()

for(i in unique(gene_table$geneID)){
  
  tmp.df <- gene_table[gene_table$geneID== i,]
  tmp.df <- tmp.df[as.numeric(tmp.df$year) == min(as.numeric(tmp.df$year)),]
  
  if(nrow(tmp.df)>1){
    tmp.df <- tmp.df[as.numeric(tmp.df$month) == min(as.numeric(tmp.df$month)),]
    tmp.df <- tmp.df[as.numeric(tmp.df$day) == min(as.numeric(tmp.df$day)),]
  }
  
  unique_gene_table <- rbind(unique_gene_table,tmp.df)
  unique_gene_table <- unique_gene_table[!duplicated(unique_gene_table[,c("pmid","geneID")]),]
  
}

#remove duplicated genes (sometimes a journal published about a gene in two publications on the same day)
unique_gene_table <- unique_gene_table[!duplicated(unique_gene_table$geneID),] 

#add total gene count to table 
unique_gene_table <- merge(unique_gene_table, genecounts_df, by.x = "geneID", by.y = "geneID")

#only keep human genes and translate geneIDs to HUGO symbol
unique_gene_table <- unique_gene_table[!is.na(match(unique_gene_table$geneID, gene_df$GeneID)),] 
unique_gene_table <- merge(unique_gene_table, gene_df, by.x = "geneID", by.y = "GeneID")

#calculate average of publications per year
unique_gene_table$year <- as.numeric(unique_gene_table$year)
unique_gene_table$abstracts_year <- unique_gene_table$freq / ((2021 - unique_gene_table$year)+1) 
unique_gene_table <- unique_gene_table[order(unique_gene_table$freq, decreasing=T),]
colnames(unique_gene_table) <- c("ID","PMID","Gene","doi","title","year","month","day","jabbrv","journal","freq","Symbol", "abstracts_year")

#annotate established epilepsy genes
unique_gene_table$epilepsygene <- unique_gene_table$Symbol %in% as.character(epilepsy_genes$gene)
unique_gene_table$epilepsygene <- ifelse(unique_gene_table$epilepsygene == TRUE, "epilepsygene (Heyne et al., Lindy et al., ClinGen)", "no epilepsygene")

#filter to genes to plot that have at least one publication per year on average
plot_data <- as.data.frame(unique_gene_table[unique_gene_table$abstracts_year > 1| unique_gene_table$epilepsygene == "epilepsygene (Heyne et al., Lindy et al., ClinGen)",])


###write top genes and known epilepsy genes to excel file for manual check (check genes for PubTator errors and year of first reference)
writexl::write_xlsx(plot_data, "data_gene_year_analysis/epilepsy_genes_manual_check.xlsx")


#### STEP 4: Manual control of epilepsy associated genes and their year of first reference ####

#### STEP 5: Plot timeline of associated genes ####

###read excel file to plot data
plot_data <- readxl::read_excel("data_gene_year_analysis/epilepsy_genes_manual_checked.xlsx")

#modify based on manual control
plot_data = plot_data[plot_data$`EXCLUDE?` == "no", ]
plot_data$year <- ifelse(is.na(plot_data$`NEW YEAR`), plot_data$year, plot_data$`NEW YEAR`)

#update epilepsy gene (update with clingen version 10th of January 22)
plot_data$epilepsygene <- plot_data$Symbol %in% as.character(epilepsy_genes$gene)
plot_data$epilepsygene <- ifelse(plot_data$epilepsygene == TRUE, "epilepsygene (Heyne et al., Lindy et al., ClinGen)", "no epilepsygene")

#y-axis value for plot
plot_data <- plot_data[order(plot_data$year),]
plot_data$cumulative <- seq(1,nrow(plot_data))

plot <- ggplot(plot_data, aes(x= year, y= cumulative))+
  geom_line(size=1)+
  geom_label_repel(aes(fill=freq, label=Symbol, col=epilepsygene), size = 3.25, max.overlaps = 200)+
  ylab("Number of epilepsy-associated genes") +
  xlab("Year of first discovery") +
  theme_light() +
  labs(fill = "Number of references")+
  theme(text = element_text(size=18))+
  scale_x_continuous(limits=c(1990,2022),breaks=c(1990,1995,2000,2005,2010,2015,2020))+
  scale_y_continuous(limits=c(0,185),breaks=c(0,25,50,75,100,125,150,175))+
  scale_fill_gradientn(colours = c("red","yellow"),values=c(1.0, 0.1, 0),na.value = NA)+
  scale_color_manual(values=c("#000000", "#999999"))+
  guides(col=FALSE)+
  theme(
    legend.position = c(.95, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )


pdf("data_gene_year_analysis/epilepsy_genes_timeline.pdf",width= 15,height= 9)
plot
dev.off()
