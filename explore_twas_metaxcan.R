twas_fdr <- read.table("twas/twas_results/twas_merged_fdr_results.tab", sep="\t", header=TRUE, stringsAsFactors = FALSE)
twas_bonf <- read.table("twas/twas_results/twas_merged_bonf_results.tab", sep="\t", header=TRUE, stringsAsFactors = FALSE)

metaxcan_fdr <- read.table("metaxcan/metaxcan_se_results/metaxcan_merged_fdr_results.tab", sep="\t", header=TRUE, stringsAsFactors = FALSE)
metaxcan_bonf <- read.table("metaxcan/metaxcan_se_results/metaxcan_merged_bonf_results.tab", sep="\t", header=TRUE, stringsAsFactors = FALSE)

twas_fdr_unique_genes <- unique(twas_fdr$ID)
twas_bonf_unique_genes <- unique(twas_bonf$ID)

metaxcan_fdr_unique_genes <- unique(metaxcan_fdr$gene_name)
metaxcan_bonf_unique_genes <- unique(metaxcan_bonf$gene_name)

overlap_fdr <- twas_fdr_unique_genes[twas_fdr_unique_genes%in%metaxcan_fdr_unique_genes]
overlap_bonf <- twas_bonf_unique_genes[twas_bonf_unique_genes%in%metaxcan_bonf_unique_genes]

"VAMP4"%in%overlap_fdr
"VAMP4"%in%overlap_bonf

"LSM7"%in%twas_fdr_unique_genes
"LSM7"%in%twas_bonf_unique_genes
"LSM7"%in%metaxcan_fdr_unique_genes
"LSM7"%in%metaxcan_bonf_unique_genes

"SPPL2B"%in%twas_fdr_unique_genes
"SPPL2B"%in%twas_bonf_unique_genes
"SPPL2B"%in%metaxcan_fdr_unique_genes
"SPPL2B"%in%metaxcan_bonf_unique_genes

"PM20D1"%in%twas_fdr_unique_genes
"PM20D1"%in%twas_bonf_unique_genes
"PM20D1"%in%metaxcan_fdr_unique_genes
"PM20D1"%in%metaxcan_bonf_unique_genes

### Replicating with BONF in both: CD38, GPNMB, NUPL2, RAB7L1, VAMP4
### Replicating with FDR in both: PM20D1
### Replicating only in MetaXcan (BONF & FDR): LSM7
### Not replicating at all: SPPL2B