library(dplyr)
library(stringr)

setwd("~/git/PD_coloc")
coloc_top <- read.table("coloc_top_snps_per_gene.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
braineac <- read.table("braineac_association_results_coloc_top_snps.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

for (i in 1:length(coloc_top$Gene)){
  gene_p <- braineac %>%
    filter(geneSymbol == coloc_top$Gene[i], rsid ==  coloc_top$Coloc.Top.SNP[i], str_detect(exprID, '^t'))%>%
    select(get(coloc_top$Brain.Region[i]))%>%
    min()
  exon_top_p <- braineac %>%
    filter(geneSymbol == coloc_top$Gene[i], rsid ==  coloc_top$Coloc.Top.SNP[i], !str_detect(exprID, '^t'))%>%
    select(get(coloc_top$Brain.Region[i])) %>%
    min()
  
  coloc_top$Gene_p_value[i] <- gene_p
  coloc_top$Top_exon_p_value[i] <- exon_top_p
}
write.table(coloc_top, "coloc_top_snp_gene_vs_exon.tab", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
