library(biomaRt)
library(tidyverse)
library(patchwork)
library(edgeR)

setwd("~/Documents/_Altman/_project/COLLAB/SH2B3")

ensembl <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart='ensembl')

key <- biomaRt::getBM(attributes=c("ensembl_gene_id", "ensembl_exon_id", "mgi_symbol"), 
                          filters = "mgi_symbol",
                          values = "Sh2b3",
                          mart = ensembl)

key.loc <- biomaRt::getBM(attributes=c("ensembl_exon_id", "chromosome_name",
                                       "exon_chrom_start", "exon_chrom_end",
                                       "rank"), 
                          filters = "ensembl_exon_id",
                          values = key$ensembl_exon_id,
                          mart = ensembl)

key.all <- full_join(key, key.loc)

attach("voom_exp1.RData")
attach("voom_exp2.RData")

meta <- voom.exp1$targets %>% 
  dplyr::select(libID, Animal, Genotype, Time, Batch) %>% 
  mutate(Exp = "Experiment 1")
meta <- voom.exp2$targets %>% 
  dplyr::select(libID, Animal, Genotype, Time, Batch) %>% 
  mutate(Exp = "Experiment 2") %>% 
  bind_rows(meta)

count <- read_tsv("combined_feature_counts_exon.tsv") 

cpm <- count %>% 
  column_to_rownames("Geneid") %>% 
  edgeR::cpm(log=FALSE) %>% 
  as.data.frame() %>% 
  rownames_to_column("Geneid") %>% 
  pivot_longer(-Geneid) %>% 
  inner_join(key.all, by=c("Geneid"="ensembl_exon_id")) %>% 
  # group_by(name, mgi_symbol, rank) %>% 
  # summarise(total_cpm = sum(value, na.rm=TRUE)) %>% 
  # ungroup() %>% 
  inner_join(meta, by=c("name"="libID"))

non0 <- cpm %>% 
  group_by(Geneid) %>% 
  summarise(max_cpm=max(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(max_cpm>0)
  
cpm_format <- cpm %>% 
  filter(Geneid %in% non0$Geneid) %>% 
  mutate(Time = gsub("\\+","\n\\+", Time),
         x = paste(Genotype, Time, sep="\n"),
         x = factor(x, levels=c(
           "WT\n0hr","WT\n2hr","WT\n6hr","WT\n24hr","WT\n24hr\n+antiIL2",
           "KO\n0hr","KO\n2hr","KO\n6hr","KO\n24hr","KO\n24hr\n+antiIL2")),
         Geneid = paste0("Chr", chromosome_name, "\n",
                         exon_chrom_start, "-", exon_chrom_end, "\n",
                         Geneid),
         Geneid = fct_rev(Geneid))

p1 <- cpm_format %>% 
  mutate(Genotype = factor(Genotype, levels=c("WT","KO"))) %>% 
  ggplot() +
  aes(x=Genotype, y=value) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Exp, shape=Exp), width=0.2, height=0) +
  facet_wrap(~Geneid, scales="free", ncol=4) +
  theme_bw() +
  labs(x="", y="Counts per million (CPM)", shape="", color="") +
  scale_color_manual(values=c("#E69F00","#56B4E9"))

# p1

ggsave(p1, file = "SH2B3_exon_expression.pdf", width=9, height=5)
