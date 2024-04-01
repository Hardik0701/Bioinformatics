#A follow along attempt to learn Bioinformatic using 'R' 
#Credits : YouTube : Bioinformagician 
#Playlist name : Bioinformatics 101

#Description:
# Form the NCBI GEO breast cancer datasets the expression of various genes associated with 
# breast cancer is assessed 


#importing library

library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)

# reading data downloaded from NCBI GEO datasets 
dat <- read.csv(file = "C:/Users/Hitesh/Downloads/GSE183947_fpkm.csv/GSE183947_fpkm.csv")

# Getting meta data of cancer dataset 
Gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(Gse[[1]]))

metadata.modified <- metadata %>% select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue:","", tissue)) %>%
  mutate(metastasis = gsub("metastasis:","",metastasis))

#reshaping data
dat.long <- dat %>% 
  rename(gene = X ) %>%
  gather(key = 'samples', value = 'FPKM', -gene)

#joining 2 dataframe(left join)

dat.long <- dat.long %>%
  left_join(., metadata.modified, by =c("samples"="description"))

dat.long %>% 
  filter(gene == 'BRCA1'| gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM))%>%
  arrange(mean_FPKM)

# plotting
  
# 1. Barplot 
#Description : Compares the expression of BRCA1 gene among cancer and normal tissue 

dat.long %>%
  filter( gene =='BRCA1') %>%
  ggplot(., aes(x = samples, y = FPKM, fill=tissue))+
  geom_col()

#2 density
#Description: Density distribution of BRCA1 gene among malignant and normal tissue

dat.long %>%
  filter( gene =='BRCA1') %>%
  ggplot(., aes(x = FPKM, fill=tissue))+
  geom_density(alpha = 0.3)

#3. box plot
#Description: Assessing expression of BRCA1 based on metastasis of tumor

dat.long %>%
  filter( gene =='BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = 'red'))+
  geom_violin()

# scatter plot
#Description : Comparative expression of BRCA1 and BRCA2 base on malignant and normal tissue

dat.long %>%
  filter(gene =='BRCA1'| gene =='BRCA2')%>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x= BRCA1, y =BRCA2, color = tissue))+
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE)

#heat-map (multiple gene)
#Description: Assessing expression of multiple genes across samples (malignant and normal)
gene.of.intrest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

dat.long %>% 
  filter(gene %in% gene.of.intrest)%>%
  ggplot(., aes(x=samples, y =gene, fill= FPKM))+
  geom_tile()+
  scale_fill_gradient(low = "black",high = "red")