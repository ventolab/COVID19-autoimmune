# load required packages

library(data.table)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(plyr)
library(purrr)
library(viridis)
library(ggpubr)
library(tibble)

setwd("/disco1/Project_scCOVID/")

# load joint DEG lists

deg_list <- list(
  'ctrl' = fread("lists/deg/Joint_DE_table_Control_vs_rest_filtered.csv"), 
  'MS' = fread("lists/deg/Joint_DE_table_MS_vs_rest_filtered.csv"),
  'PS' = fread("lists/deg/Joint_DE_table_PS_vs_rest_filtered.csv"),
  'RA' = fread("lists/deg/Joint_DE_table_RA_vs_rest_filtered.csv"))

# run GSEA enrichment on the different gmt gene set files

cells <- c("CD14 mono", 'CD4 naïve T', 'CD4 memory T')
gsea_results <- list()
gsea_list <- list()
gsea_list_dis <- list()

gene_sets <- list(
  'hallmark' = read.gmt(gmtfile = "lists/GSEA_genesets/h.all.v7.4.symbols.gmt"),
  'curated' = read.gmt(gmtfile = "lists/GSEA_genesets/c2.all.v7.4.symbols.gmt"),
  'GO' = read.gmt(gmtfile = "lists/GSEA_genesets/c5.all.v7.4.symbols.gmt"),
  'immune' = read.gmt(gmtfile = "lists/GSEA_genesets/c7.all.v7.4.symbols.gmt"))

saveRDS(gene_sets, file = "lists/GSEA_genesets/gene_sets.RDS")

for (i in 1:length(cells)) {
  for (e in 1:length(deg_list)) {
    
  genelist <- deg_list[[e]] %>% 
    dplyr::filter(cluster == cells[i]) %>% 
    dplyr::select(Gene, logFC) %>%
    arrange(desc(logFC)) %>%
    deframe()
  
  gene <- names(genelist)[abs(genelist) > 0]

  for (j in 1:length(gene_sets)) {
    gsea <- GSEA(genelist, TERM2GENE=gene_sets[[j]], verbose=FALSE, pvalueCutoff = 1, nPermSimple = 10000) 
    if (nrow(gsea) == 0) {
      gsea_list[[j]] <- data.frame(
        ID = NA, NES = NA, p.adjust = NA, 
        cell = cells[i], disease = names(deg_list)[e], gene_set = names(gene_sets)[j]) 
      } else {
          gsea_list[[j]] <- gsea %>%
            as.data.frame() %>%
            dplyr::select(ID, NES, p.adjust) %>% 
            mutate(cell = cells[i], disease = names(deg_list)[e], gene_set = names(gene_sets)[j])
      }
    names(gsea_list)[j] <- names(gene_sets)[j]
      }
  gsea_list_dis[[e]] <- gsea_list
  names(gsea_list_dis)[e] <- names(deg_list)[e]
  }
  gsea_results[[i]] <- gsea_list_dis
  names(gsea_results)[i] <- cells[i]
}

# format results

gsea_table <- data.frame()

for (i in 1:length(gsea_results)){
for (j in 1:length(gsea_results[[1]])) {
for (k in 1:length(gsea_results[[1]][[1]])){
  gsea_table <- rbind(gsea_table, gsea_results[[i]][[j]][[k]])
    }
  }
}

fwrite(gsea_table, "lists/gsea_table_results.txt", sep = '\t')
gsea_table <- fread("lists/gsea_table_results.txt")

# plot results

plot <- gsea_table %>%
  mutate(logFDR = -log10(p.adjust), 
         cell = factor(cell, levels = c("CD14 mono", 'CD4 naïve T', 'CD4 memory T'))) %>% 
  dplyr::filter(gene_set %in% c('hallmark','GO')) %>%
  mutate(NES_bimodal = ifelse(NES > 0, 1, -1), NES_abs = abs(NES), 
         logFDR = ifelse(logFDR < -log10(0.05), 0.5, logFDR)) %>%
  na.omit() 

categories <- c(
  'REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT',
  "GOBP_INTERLEUKIN_1_BETA_PRODUCTION", 
  'GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION', 
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  'HALLMARK_HYPOXIA',
  'GOBP_RESPONSE_TO_VIRUS',
  'GOBP_RESPONSE_TO_TYPE_I_INTERFERON',
  'GOBP_RESPONSE_TO_INTERFERON_GAMMA',
  'GOMF_MHC_PROTEIN_COMPLEX_BINDING',
  'REACTOME_INTERLEUKIN_17_SIGNALING',
  'GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE',
  'GOBP_POSITIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE', 
  'GOBP_RESPONSE_TO_INTERLEUKIN_12')

ggscatter(plot %>% dplyr::filter(ID %in% categories), size = 'logFDR',
          'disease', 'ID', facet.by = 'cell', color = 'NES',
          order = c('ctrl', 'RA', 'MS', 'PS'), ylab = '', xlab = '') %>%
  ggpar(legend.title = list(size = expression("-log"[10]*'FDR'))) +
  scale_size(limits = c(0,8)) +
  scale_color_distiller(palette = 'RdYlGn', na.value = 'gray') +
  theme_bw()

ggplot(plot %>% dplyr::filter(ID %in% categories), aes(disease, ID, fill = NES)) +
  geom_tile() +
  scale_fill_viridis(option = 'E') +
  theme_bw() +
  facet_grid(~cell)

ggsave("figures/gsea_bubble_square.pdf", height = 3.5, width = 9)
ggsave("figures/gsea_bubble_RdYlGn.pdf", height = 3.5, width = 9)