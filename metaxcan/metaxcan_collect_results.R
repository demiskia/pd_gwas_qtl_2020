library(stringr)
setwd("metaxcan_se_results/")

tests <- list.files()
for (i in 1:length(tests)){
  assign(tests[i], read.table(tests[i], sep=",", header=TRUE, stringsAsFactors=FALSE))
}
results <- setNames(lapply(tests, function(x) get(x)), tests)

for(j in 1:length(results))
{
  temp <- results[[j]]
  temp$Pvalue <- as.numeric(as.character(temp$pvalue))
  dat <- subset(temp, !is.na(temp$Pvalue))
  dat$Pvalue_FDR <- p.adjust(dat$Pvalue, method = c("fdr"))
  dat$Pvalue_bonferroni <- p.adjust(dat$Pvalue, method = c("bonferroni"))
  dat <- dat[order(dat$Pvalue),]
  dat$Region <- names(results)[j]
  dat$Region <- str_replace_all(dat$Region, "_0.5_1KG_results_se.txt", "")
  dat$Region <- str_replace_all(dat$Region, "TW_", "")
  dat_fdr <- subset(dat, Pvalue_FDR <= 0.05)
  dat_bonf <- subset(dat, Pvalue_bonferroni <= 0.05)
  write.table(dat_fdr, file = paste(names(results)[j],"_FDR_HITS.tsv",sep = ""), quote = F, sep = "\t", row.names = F, col.names=TRUE)
  write.table(dat_bonf, file = paste(names(results)[j],"_BONF_HITS.tsv",sep = ""), quote = F, sep = "\t", row.names = F, col.names=TRUE)
}

fdrs <- list.files(pattern = "FDR_HITS")
for (i in 1:length(fdrs)){
  assign(fdrs[i], read.table(fdrs[i], sep="\t", header=TRUE, stringsAsFactors=FALSE))
}
fdr_results <- setNames(lapply(fdrs, function(x) get(x)), fdrs)

full_fdr_results <- NULL
for(j in 1:length(fdr_results))
{
  full_fdr_results <- rbind(full_fdr_results, fdr_results[[j]])
}
full_fdr_results <- full_fdr_results[order(full_fdr_results$pvalue),]
write.table(full_fdr_results, file = "metaxcan_merged_fdr_results.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)


bonferronis <- list.files(pattern = "BONF_HITS")
for (i in 1:length(bonferronis)){
  assign(bonferronis[i], read.table(bonferronis[i], sep="\t", header=TRUE, stringsAsFactors=FALSE))
}
bonf_results <- setNames(lapply(bonferronis, function(x) get(x)), bonferronis)
full_bonf_results <- NULL
for(j in 1:length(bonf_results))
{
  full_bonf_results <- rbind(full_bonf_results, bonf_results[[j]])
}
full_bonf_results <- full_bonf_results[order(full_bonf_results$pvalue),]
write.table(full_bonf_results, file = "metaxcan_merged_bonf_results.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)
