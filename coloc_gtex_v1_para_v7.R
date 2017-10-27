library(coloc)
library(tidyverse)
library(stringr)

START = Sys.getenv("START")
END = Sys.getenv("END")

# Functions -------------------------------------------------------------------------------------------

generate_parameters_df <- function(p1_or_p2, prob_shared, PD_cc_ratio){
  
  # parameters_df <- data_frame()
  
  priors_df <- 
    data_frame(p1 = rep(p1_or_p2, length(prob_shared)) %>% sort(decreasing = T), 
               prob_shared = rep(prob_shared, length(p1_or_p2))) %>% 
    mutate(p2 = p1, 
           p12 = p1 * prob_shared, 
           priors = str_c(p1, "_", p12), 
           index = 1:length(p2)) %>% 
    select(index, priors, p1, p2, p12, prob_shared)
  
  # coloc.abf does not require cc ratio as input
  # for(i in seq_along(PD_cc_ratio)){
  #   
  #   parameters_df <- bind_rows(parameters_df, priors_df)
  #   
  # }
  # 
  # parameters_df <- 
  # parameters_df %>% 
  #   mutate(PD_cc_ratio = rep(PD_cc_ratio, nrow(parameters_df)/length(PD_cc_ratio)))
  
  # return(parameters_df)
  
  return(priors_df)
  
}

sort_GTEx_brain_regions_eQTL_paths <- function(GTEx_brain_regions_eQTL_paths, coloc_PD_genes_HGNC_symbols_ensembl_38_37_uniq){
  
  GTEx_brain_regions_eQTL_df <- data_frame(path = GTEx_brain_regions_eQTL_paths, 
                                           tissue = str_replace_all(path, "/.*/|_Analysis.*","") %>% str_replace_all("_", " "), 
                                           total_PD_genes_num = coloc_PD_genes_HGNC_symbols_ensembl_38_37_uniq$hgnc_symbol %>% unique() %>% length(), 
                                           genes_present_num = NA,
                                           genes_absent_num = NA,
                                           genes_absent_symbols = NA)
  
  for(i in seq_along(GTEx_brain_regions_eQTL_df$path)){
    
    path <- GTEx_brain_regions_eQTL_df$path[i]
    GTEx_brain_region_eQTL <- read_delim(path, delim = "\t")
    
    genes_present <- 
      coloc_PD_genes_HGNC_symbols_ensembl_38_37_uniq %>% 
      filter(ensembl_gene_id %in% GTEx_brain_region_eQTL$gene_id_clean) %>% 
      .[["hgnc_symbol"]] %>% 
      unique()
    
    genes_absent <- 
      coloc_PD_genes_HGNC_symbols_ensembl_38_37_uniq %>% 
      filter(!hgnc_symbol %in% genes_present) %>% 
      .[["hgnc_symbol"]] %>% 
      unique()
    
    GTEx_brain_regions_eQTL_df$genes_present_num[i] <- length(genes_present)

    if(length(genes_absent) > 0){
      
      GTEx_brain_regions_eQTL_df$genes_absent_num[i] <- length(genes_absent)
      GTEx_brain_regions_eQTL_df$genes_absent_symbols[i] <- str_c(genes_absent, collapse = ", ")
      
    }else{
      
      GTEx_brain_regions_eQTL_df$genes_absent_num[i] <- 0
      
    }
    
  }
  
  return(GTEx_brain_regions_eQTL_df)
  
}

format_PD_GWAS <- function(PD_GWAS_v2){
  
  PD_GWAS_v2_tidy <- 
    PD_GWAS_v2 %>% 
    select(variant_id, BETA = Effect_flip, SE = StdErr, Freq1_flip, Allele1_flip, Allele2_flip) %>% 
    mutate(BETA = as.double(BETA), 
           SE = as.double(SE), 
           varbeta = SE ^ 2,
           Freq1_flip = as.double(Freq1_flip), 
           MAF = ifelse(Freq1_flip > 0.5, 1 - Freq1_flip, Freq1_flip), 
           SE = NULL, 
           Freq1_flip = NULL)
    
  return(PD_GWAS_v2_tidy)
  
}

make_results_dir <- function(results_path, tissue_formatted){
  
  dir_path <- str_c(results_path, "/", tissue_formatted)
  
  if(!dir.exists(dir_path)){
    
    dir.create(dir_path)
    
  }else {
    
    print(str_c(dir_path, " directory already exists.."))
    
  }
  
}

import_GTEx_eQTL <- function(path){
  
  GTEx_PD_genes_eQTL_data <- read_delim(path, delim = "\t", col_types = cols(.default = col_character()))
  
  return(GTEx_PD_genes_eQTL_data)
  
}

format_GTEx_eQTL <- function(GTEx_PD_genes_eQTL_data){
  
  GTEx_PD_genes_eQTL_data_tidy <- 
    GTEx_PD_genes_eQTL_data %>% 
    select(gene_id_clean, variant_id, BETA = slope, SE = slope_se) %>% 
    mutate(BETA = as.double(BETA),
           SE = as.double(SE), 
           varbeta = SE ^ 2, 
           t.stat = NULL,
           SE = NULL)
    
  return(GTEx_PD_genes_eQTL_data_tidy)
  
}

# courtesy of moloc package
match_alleles <- function(data, A1.ref="A1.ref", A2.ref="A2.ref",  A1.data = "A1", A2.data = "A2", BETA.data="BETA", flip = TRUE) {
  
  match_correct = data[,A1.ref] == data[,A1.data] & data[,A2.ref]== data[,A2.data]
  match_flip = data[,A1.ref] == data[,A2.data] & data[,A2.ref] == data[,A1.data]
  match_comp_one = data[,A1.ref] == complement_snp(data[,A1.data]) & data[,A2.ref]== complement_snp(data[,A2.data])
  match_comp_two = data[,A1.ref] == complement_snp(data[,A2.data]) & data[,A2.ref] == complement_snp(data[,A2.data])
  snp_allele_match = match_flip | match_correct | match_comp_one | match_comp_two
  message(sum(snp_allele_match), " SNPs out of ", length(snp_allele_match), " had the correct alleles, discarding SNPs without the correct alleles")
  if (flip) {
    if (any(which(match_flip)>0)) {
      data[match_flip, A1.data]=data[match_flip, A1.ref]
      data[match_flip, A2.data]=data[match_flip, A2.ref]
      data[match_flip, BETA.data]=-data[match_flip, BETA.data]
    }
  }
  removed = data[!snp_allele_match,]
  data = data[snp_allele_match,]
  return(list(data, removed))
}

complement_snp <- function(x){
  
  as = x =="A"
  ts = x == "T"
  gs = x == "G"
  cs = x == "C"
  ins = x == "I"
  dels = x == "D"
  x[as] = "T"
  x[ts] = "A"
  x[gs] = "C"
  x[cs] = "G"
  x[ins] = "NA"
  x[dels] = "NA"
  
  return(x)
}

change_indels <- function(data) {
  
  data$A2[nchar(data$A1)>1] = "D"
  data$A1[nchar(data$A1)>1] = "I"
  data$A1[nchar(data$A2)>1] = "D"
  data$A2[nchar(data$A2)>1] = "I"
  
  return(data)
}

run_coloc_analysis <- function(PD_GWAS_v2_tidy, GTEx_PD_genes_eQTL_data_tidy_gene_filtered, priors_df, 
                               tissue_formatted, gene_to_filter){
  
  
  PD_GWAS_v2_coloc <- 
    PD_GWAS_v2_tidy %>% 
    filter(variant_id %in% GTEx_PD_genes_eQTL_data_tidy_gene_filtered$variant_id) %>% 
    arrange(variant_id)
  
  GTEx_PD_genes_eQTL_data_coloc <- 
    GTEx_PD_genes_eQTL_data_tidy_gene_filtered %>% 
    filter(variant_id %in% PD_GWAS_v2_tidy$variant_id, !duplicated(variant_id)) %>% 
    arrange(variant_id)
  
  print(str_c(tissue_formatted, " - ", gene_to_filter, " - filtered ;)" ))
  
  stopifnot(identical(PD_GWAS_v2_coloc$variant_id, GTEx_PD_genes_eQTL_data_coloc$variant_id))
  
  GTEx_PD_genes_eQTL_data_coloc <- 
  GTEx_PD_genes_eQTL_data_coloc %>% 
    mutate(REF_PD_Al1 = PD_GWAS_v2_coloc$Allele1_flip, 
           REF_PD_Al2 = PD_GWAS_v2_coloc$Allele2_flip)
  
  if(length(GTEx_PD_genes_eQTL_data_coloc$variant_id) == 0){
    
    print(str_c(tissue_formatted, " - ", gene_to_filter, " - no overlap :(" ))
    
    return(NULL)
    
  }
  
  GTEx_PD_genes_eQTL_data_coloc_match_alleles_list <- 
  match_alleles(GTEx_PD_genes_eQTL_data_coloc, A1.ref="REF_PD_Al1", A2.ref="REF_PD_Al2",  
                A1.data = "Al1", A2.data = "Al2", BETA.data="BETA", flip = TRUE)
  
  write_delim(GTEx_PD_genes_eQTL_data_coloc_match_alleles_list[[2]], 
              str_c("/data/kronos/dzhang/coloc/results/GTEx/GTEx_raw_results_full_v7", "/", tissue_formatted, "/", tissue_formatted, "_", 
                    gene_to_filter, "_removed_SNPS.csv"), delim = ",")
  
  GTEx_PD_genes_eQTL_data_coloc_harmonised <- 
    GTEx_PD_genes_eQTL_data_coloc_match_alleles_list[[1]]
  
  PD_GWAS_v2_coloc <- 
    PD_GWAS_v2_coloc %>% 
    filter(variant_id %in% GTEx_PD_genes_eQTL_data_coloc_harmonised$variant_id) %>% 
    arrange(variant_id)
  
  GTEx_PD_genes_eQTL_data_coloc_harmonised <- 
    GTEx_PD_genes_eQTL_data_coloc_harmonised %>% 
    filter(variant_id %in% PD_GWAS_v2_coloc$variant_id)%>% 
    arrange(variant_id) 
  
  print(str_c(tissue_formatted, " - ", gene_to_filter, " - harmonised ;)" ))
  
  stopifnot(identical(PD_GWAS_v2_coloc$Allele1_flip, GTEx_PD_genes_eQTL_data_coloc_harmonised$Al1))
  stopifnot(identical(PD_GWAS_v2_coloc$Allele2_flip, GTEx_PD_genes_eQTL_data_coloc_harmonised$Al2))
  
  coloc_results_gene <- format_coloc_priors(PD_GWAS_v2_coloc, GTEx_PD_genes_eQTL_data_coloc_harmonised, priors_df)

  return(coloc_results_gene)
  
}

format_coloc_priors <- function(PD_GWAS_v2_coloc, GTEx_PD_genes_eQTL_data_coloc_harmonised, priors_df){
  
  coloc_results_gene <- list()
  
  for(i in seq_along(priors_df$index)){
    
    p1 <- priors_df$p1[i]
    p2 <- priors_df$p2[i]
    p12 <- priors_df$p12[i]
    priors <- priors_df$priors[i]
    # PD_cc_ratio <- parameters_df$PD_cc_ratio[i]
    
    coloc_results_prior <- 
      coloc.abf(dataset1 = list( beta = GTEx_PD_genes_eQTL_data_coloc_harmonised$BETA, 
                                 varbeta = GTEx_PD_genes_eQTL_data_coloc_harmonised$varbeta, 
                                 N = 135, 
                                 type = "quant", 
                                 snp = GTEx_PD_genes_eQTL_data_coloc_harmonised$variant_id),
                                 # MAF = GTEx_PD_genes_eQTL_data_coloc_harmonised$MAF), 
                dataset2 = list( beta = PD_GWAS_v2_coloc$BETA,
                                 varbeta = PD_GWAS_v2_coloc$varbeta,
                                 type = "cc", 
                                 # s = PD_cc_ratio,
                                 snp = PD_GWAS_v2_coloc$variant_id),
                                 MAF = PD_GWAS_v2_coloc$MAF, 
                p1 = p1, p2 = p2, p12 = p12)
   
    coloc_results_prior[["summary"]]["p1"] <- p1
    coloc_results_prior[["summary"]]["p2"] <- p2
    coloc_results_prior[["summary"]]["p12"] <- p12
    
    names(coloc_results_prior) <- str_c(names(coloc_results_prior), "_", priors)
    
    coloc_results_gene[[i]] <- coloc_results_prior
     
  }
  
  return(coloc_results_gene)
  
}

save_coloc_results <- function(coloc_results_gene, results_path, tissue_formatted, PD_gene){
  
  results_path_gene_tissue_exprID <- str_c(results_path, "/", tissue_formatted, "/", tissue_formatted, "_", PD_gene)

  for(i in seq_along(coloc_results_gene)){
    
    if(i == 1){
      
      coloc_results_gene_summary <- coloc_results_gene[[i]][[1]] %>% t() %>% as_tibble()
      
    }else{
      
      coloc_results_gene_summary <- 
        bind_rows(coloc_results_gene_summary,
                  coloc_results_gene[[i]][[1]] %>% t() %>% as_tibble())
      
    }
    
    write_delim(coloc_results_gene[[i]][[2]],
                str_c(results_path_gene_tissue_exprID, "_coloc_", names(coloc_results_gene[[i]][2]), ".csv"), 
                delim = ",")
    
  }
  
  coloc_results_gene_summary_complete <- 
    coloc_results_gene_summary %>% 
    mutate(gene = PD_gene, 
           tissue = tissue_formatted) %>% 
    select(gene, tissue, everything())
  
  write_delim(coloc_results_gene_summary_complete, 
              str_c(results_path_gene_tissue_exprID, "_coloc_summary.csv"), 
              delim = ",")
  
}

# Main ------------------------------------------------------------------------------------------------

GTEx_dir <- "/array/dkia/eqtl/gtex/gtex_v7/Brain_regions_gtex_clean_id"

GTEx_brain_regions_eQTL_paths <- list.files(GTEx_dir, full.names = T)

 coloc_PD_genes_HGNC_symbols_ensembl_38_37 <- 
   read_delim("/data/kronos/dzhang/coloc/raw_data/GTEx/coloc_PD_genes_HGNC_symbols_ensembl_38_37.csv", delim = ",", 
              col_types = cols(.default = col_character()))

 coloc_PD_genes_HGNC_symbols_ensembl_38_37_uniq <- 
   coloc_PD_genes_HGNC_symbols_ensembl_38_37 %>% 
     mutate(HGNC_ensembl = (str_c(hgnc_symbol, "_", ensembl_gene_id))) %>%
     filter(!duplicated(HGNC_ensembl)) %>%
     as_tibble()

if(!file.exists("/data/kronos/dzhang/coloc/results/GTEx/GTEx_brain_regions_eQTL_df_v7.csv")){

  GTEx_brain_regions_eQTL_df <- sort_GTEx_brain_regions_eQTL_paths(GTEx_brain_regions_eQTL_paths, coloc_PD_genes_HGNC_symbols_ensembl_38_37)

  write_delim(GTEx_brain_regions_eQTL_df, "/data/kronos/dzhang/coloc/results/GTEx/GTEx_brain_regions_eQTL_df_v7.csv", delim = ",")

}else{

  GTEx_brain_regions_eQTL_df <- read_delim("/data/kronos/dzhang/coloc/results/GTEx/GTEx_brain_regions_eQTL_df_v7.csv", delim = ",")

}

priors_df <- generate_parameters_df(p1_or_p2 = c(1e-04, 1e-05), prob_shared = c(0.01, 0.02, 0.1), 
                              PD_cc_ratio = c(average_PD_cc_ratio, max_PD_cc_ratio, min_PD_cc_ratio, meta_ratio))

PD_GWAS_v2_tidy <- read_delim("/data/kronos/dzhang/coloc/raw_data/pd_newmeta_refallele.txt", delim = "\t",
                         col_types = cols(.default = col_character()))

PD_GWAS_v2_tidy <- format_PD_GWAS(PD_GWAS_v2_tidy)

gc(verbose=TRUE)

for(i in seq_along(GTEx_brain_regions_eQTL_df$path)){
  
  path <- GTEx_brain_regions_eQTL_df$path[i]
  tissue_formatted <- GTEx_brain_regions_eQTL_df$tissue[i] %>% str_replace_all(" ", "_")
  
  make_results_dir(results_path = "/data/kronos/dzhang/coloc/results/GTEx/GTEx_raw_results_full_v7", tissue_formatted)
  
  GTEx_PD_genes_eQTL_data_tidy <- import_GTEx_eQTL(path)
  GTEx_PD_genes_eQTL_data_tidy <- format_GTEx_eQTL(GTEx_PD_genes_eQTL_data_tidy)
  
  gc(verbose=TRUE)
  
  # PD_genes_present <- 
  #   coloc_PD_genes_HGNC_symbols_ensembl_38_37_uniq %>% 
  #   filter(ensembl_gene_id %in% GTEx_PD_genes_eQTL_data_tidy$gene_id_clean) 
  
  GTEx_ensembl_IDs_uniq <-
    GTEx_PD_genes_eQTL_data_tidy$gene_id_clean %>% unique()
  
  for(index in START:END){
    
    if(index %in% seq_along(GTEx_ensembl_IDs_uniq)){
      
      gene_to_filter <- GTEx_ensembl_IDs_uniq[index]
      
      print(str_c(tissue_formatted, " - ", index, " - ", gene_to_filter))
      
      GTEx_PD_genes_eQTL_data_tidy_gene_filtered <- 
        GTEx_PD_genes_eQTL_data_tidy %>% 
        filter(gene_id_clean == gene_to_filter) %>% 
        mutate(variant_id_tmp = variant_id) %>% 
        separate(variant_id_tmp, into = c("chr", "pos", "Al1", "Al2", "gchr"))
      # mutate(SNP = str_c("chr", chr, ":", pos), 
      #        chr = NULL, 
      #        pos = NULL, 
      #        gchr = NULL) %>% 
      # filter(!str_detect(Al1, ".."), !str_detect(Al2, ".."))
      
      coloc_results_gene <- 
        run_coloc_analysis(PD_GWAS_v2_tidy, GTEx_PD_genes_eQTL_data_tidy_gene_filtered, priors_df, 
                           tissue_formatted, gene_to_filter)
      
      if(is.null(coloc_results_gene)){
        
        print("hello")
        next
        
      }
      print(str_c(tissue_formatted, " - ", index, " - ", gene_to_filter, " finished"))
      
      save_coloc_results(coloc_results_gene, results_path = "/data/kronos/dzhang/coloc/results/GTEx/GTEx_raw_results_full_v7", tissue_formatted, gene_to_filter)
      
    }else{
      
      print(str_c(tissue_formatted, " - ", index, " - ", gene_to_filter, " - break"))
      
      break
    
    }
    
    # PD_gene_to_filter <- PD_genes_present$ensembl_gene_id[i]
    # PD_gene <- PD_genes_present$hgnc_symbol[i]
  
   
  
  }
  
  print(str_c("Mission complete: ", tissue_formatted))

gc(verbose=TRUE)  
}
