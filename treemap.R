#treemap
library(readxl)
library(treemap)

epilepsy_gene_grouping <- read_excel("data_treeplot/epilepsy_gene_grouping.xlsx",  sheet = "TabelleA")

epilepsy_genes_manual_checked <- read_excel("data_gene_year_analysis/epilepsy_genes_manual_checked.xlsx")
epilepsy_genes_manual_checked <- epilepsy_genes_manual_checked[!epilepsy_genes_manual_checked$`EXCLUDE?` == "yes",]
epilepsy_genes_manual_checked <- epilepsy_genes_manual_checked[, c("Symbol", "freq")]

#merge datasets to annotate frequency of genes
epilepsy_gene_grouping <- merge(epilepsy_gene_grouping, epilepsy_genes_manual_checked, by.x="gene", by.y="Symbol", all=T)

epilepsy_gene_grouping$subgroup[is.na(epilepsy_gene_grouping$subgroup)] = "other"


pdf("data_treeplot/treeplot.pdf",width= 15,height= 8)
treemap(epilepsy_gene_grouping,
        index=c("subgroup","subgroup2", "gene"),
        vSize="freq",
        type="index",
        fontsize.labels=c(14,11,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","black", "black"),   # Color of labels
        fontface.labels=c(1,1,3),                   # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        #bg.labels=c("transparent"),                # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("left", "top"),
          c("right", "bottom")
        ),                                          # Where to place labels in the rectangle?
        overlap.labels=1,                           # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        palette= c("#8BE064" ,
                   "#D974B6" ,
                   "#E1D65C" ,"#D5D8AD", "#71DAAC" ,"#D7B9CD", "#94CED7", "#B34DDC" ,"#D67863",
                   "#828BDB", "#FFD92F", "#4DAF4A", "#7FC97F" ,"#7570B3", "#CCCCCC", "#FDC086", "#B3E2CD", "#FFD92F" ,"#CBD5E8" ,"#E31A1C", "#999999"),
        inflate.labels= F) 
dev.off()








