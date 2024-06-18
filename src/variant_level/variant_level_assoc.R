# Load the necessary library
library(dplyr)
library(ggplot2)

plink.assoc <- read.table("../../data/variant_level/plink.assoc", header = TRUE)
head(plink.assoc)

plink.assoc$p_val <- plink.assoc$P
plink.assoc$row <- rownames(plink.assoc)
plink.assoc$row <- as.numeric(plink.assoc$row)

# psig threshold
psig <- .05/length(plink.assoc$row) 

# Filter rows where F_A and F_U are not 0
# plink.assoc.filtered <- plink.assoc |> filter(F_A != 0, F_U != 0)

p1 <- plink.assoc |>
	# sample_n(10000) |>
	ggplot( aes(x = row, y = -log10(p_val))) +
  geom_hline(yintercept = -log10(psig),
             linetype="dotted", 
             color="red") + # Bonferroni correction threshold
	geom_point() +
	theme_minimal()

# Show the plot
p1
# ggsave("../../data/variant_level/plink.assoc.pdf", p1)
ggsave("../../data/variant_level/plink_assoc.png", p1)

# install.packages("qqman",repos="http://cran.xl-mirror.nl/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
# library("qqman",lib.loc="~")  

library(qqman)
library(dplyr)

# mlma manhattan plot
# Chr	SNP	bp	A1	A2	Freq	b	se	p
# results_log <- read.table("./SPSS_SHCS_merge.filt.mlma.EU_only.loco.mlma", head=TRUE)  %>% filter(p > 0) 
plink.assoc$chr <- as.numeric(plink.assoc$CHR)
plink.assoc$bp <- as.numeric(plink.assoc$BP)
plink.assoc$p <- as.numeric(plink.assoc$P)

plink.assoc <- plink.assoc[!is.na(plink.assoc$p) & !is.infinite(plink.assoc$p), ]

manhattan(plink.assoc)

png("../../data/variant_level/plink_assoc_ppman.png", width = 5, height = 5, units = 'in', res = 600)
manhattan(plink.assoc,chr="CHR", bp="BP",p="P", 
          suggestiveline = F,
          genomewideline =  -log10(psig)
          )
dev.off()

# QQ plot
png("../../data/variant_level/plink_assoc_ppman_qq.png", width = 5, height = 5, units = 'in', res = 600)
qq(plink.assoc$p 
   #main = "Q-Q plot of GWAS p-values : log"
)
dev.off()
