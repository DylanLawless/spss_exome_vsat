library(ggplot2)
library(dplyr)
library(gridExtra)

#/////////////////////
# Eigen vector plot
#/////////////////////

# import data

val <- read.table("../../data/joint/pca_output/bcftools_gatk_norm_pca.eigenval", header=F) 
vec <- read.table("../../data/joint/pca_output/bcftools_gatk_norm_pca.eigenvec", header=F) 
pheno <- read.csv("../../data/joint/pca_output/phenotypes.csv", header=T, sep = "," )
# pheno <- phenotypes %>% select(V1, V2, age.days, study.site, gender, ethnicity)
# pheno <- phenotypes %>% select(V1, pheno)

vec_pheno <- merge(x=vec, y=pheno, by="V1", all=TRUE ) 
vec_pheno_raw <- vec_pheno

vec_pheno_raw %>% 
	ggplot(aes(x=V3, y=V4))+ 
	geom_point(aes( fill = pheno ), shape=21, alpha = 0.5, position = position_jitter(height = 0))+
	labs(x = "PC1", y = "PC2") +
	geom_hline(yintercept =.1, linetype="dotted")+
	geom_hline(yintercept =-.1, linetype="dotted")
# vec_pca_pheno_raw.pdf 4x4

vec_pheno <- vec_pheno %>% filter(V3 < .1 ) %>% filter(V3 > -.1)
# test plott
vec_pheno %>% 
  ggplot(aes(x=V3, y=V4))+ 
  geom_point(aes( fill = pheno ), shape=21, alpha = 0.5, position = position_jitter(height = 0))+
  labs(x = "PC1", y = "PC2")
# vec_pca_pheno_PC1PC2.pdf 4x4

# Set plot style
theme_set <- theme(text = element_text(face="bold"), legend.position="none", panel.background = element_rect("#F7F7F7"))
geom_set <- geom_point(aes(fill = pheno), shape=21, alpha = 0.5)

# replace with a loop function here
p1 <- vec_pheno %>% ggplot(aes(x=V3, y=V4))+geom_set + labs(x = "PC1", y = "PC2")+ theme_set
p2 <- vec_pheno %>% ggplot(aes(x=V4, y=V5))+geom_set + labs(x = "PC2", y = "PC3")+ theme_set
p3 <- vec_pheno %>% ggplot(aes(x=V5, y=V6))+geom_set + labs(x = "PC3", y = "PC4")+ theme_set
p4 <- vec_pheno %>% ggplot(aes(x=V6, y=V7))+geom_set + labs(x = "PC4", y = "PC5")+ theme_set
p5 <- vec_pheno %>% ggplot(aes(x=V7, y=V8))+geom_set + labs(x = "PC5", y = "PC6")+ theme_set
p6 <- vec_pheno %>% ggplot(aes(x=V8, y=V9))+geom_set + labs(x = "PC6", y = "PC7")+ theme_set
p7 <- vec_pheno %>% ggplot(aes(x=V9, y=V10))+geom_set + labs(x = "PC7", y = "PC8")+ theme_set
p8 <- vec_pheno %>% ggplot(aes(x=V10, y=V11))+geom_set + labs(x = "PC8", y = "PC9")+ theme_set
p9 <- vec_pheno %>% ggplot(aes(x=V11, y=V12))+geom_set + labs(x = "PC9", y = "PC10")+ theme(text = element_text(face="bold"), legend.position="right", panel.background = element_rect("#F7F7F7"))

grid.arrange (p1,p2,p3,p4,p5,p6,p7,p8,p9, ncol=4, top = "Cohort Principal component analysis", left = "")
# vec_pca_pheno_grid.pdf 12x12

#/////////////////////
# Eigen value plot
#/////////////////////

# Give row name column to plot
val <- val %>% mutate(V2 = rownames(val))

# Plot
val %>% 
  ggplot(aes(x = as.numeric(V2), y = V1))+  
  geom_point()
# vec_pca_val.pdf 4x4


