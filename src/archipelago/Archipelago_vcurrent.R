
# VSAT_level ----
source("../post_ppi/vsat_log_result_vcurrent.R")
df1 <- skat_log |> dplyr::select(MCL_ID,p_val)
df1$set_ID <- df1$MCL_ID
df1$P <- df1$p_val
df1$P <- -log10(df1$P)
df1_ready <- df1 |> dplyr::select(set_ID, P)

# custom cleanup ----
# qc list
# 22.log: C1orf162 MNDA STAC CSF1R LY86 TREM1 TREM2 TREML1 TREML2 TREML4 CLEC5A MS4A6A MS4A7 CLEC4A KLRC2 KLRD1 IGSF6 CD300E CD300LB CD300LF KIR2DL1 KIR2DL4 KIR3DL1 KIR3DL2 KIR3DL3 LILRA1 LILRB1 LILRB4 SIGLEC14 SIGLEC7 TYROBP SIGLEC1 SIRPB1 SIRPG 

# 586.log: FCGR2A FCGR2B FCGR3A FCGR3B FCAR LILRA5 LILRA6 LILRB5 

# 836.log:
  # nucleosomal DNA binding
  # histone deacetylase binding
  # RNA polymerase II-specific DNA-binding transcription factor binding
  # ARID4B
  # ETV3
  # GATAD2B
  # REST
  # PRDM9
  # HDAC2
  # BRMS1
  # CHD4
# ARID4A
# CHD3
# SAP30BP
# GATAD2A
# HDAC10

# variant_level -----
source("../variant_level/variant_level_assoc.R")
head(plink.assoc)
df2 <- plink.assoc |> select(BP, P, CHR, F_A, F_U)
df2$SNP <- df2$BP

# One could remove the small subset of extreeme single variants which have an odds ratio of 0 which are benign, and not enriched in the single-variant analysis, but distracting since the are above the VSAT line. 
# I will keep them in version one - if reviewers find that confusing, it would be a reasonable function to add to the archipelago package.
# df2 <- plink.assoc |> filter(!(P < 1e-4 & (OR == 0) ))

# Get the variant CHR, BP, MCL_ID ----
vsat_log_v4_result_snp_mac <- read.table("../../data/post_ppi/vsat_log_v4_result_snp_mac.txt", fill = TRUE, header=TRUE)
# vsat_log_v4_result_snp_mac_bak.txt

# Separate SNP column
vsat_log_v4_result_snp_mac <- vsat_log_v4_result_snp_mac |>
	separate(SNP, into = c("CHR", "variant"), sep = ":", remove = TRUE)  |>
	separate(variant, into = c("BP", "variant"), sep = "_", remove = TRUE)

# Load stringr package
library(stringr)

# Remove "chr" from the "CHR" column
vsat_log_v4_result_snp_mac$CHR <- str_replace(vsat_log_v4_result_snp_mac$CHR, "chr", "")

# Merge CHR, BP, MCL_ID with variant_level result to get MCL_ID
df2_ready <- merge(df2, vsat_log_v4_result_snp_mac, all.x = T)
df2_ready$MCL_ID <- df2_ready$MCL_ID |> as.numeric()
df2_ready$set_ID <- df2_ready$MCL_ID
df2_ready <- df2_ready |> filter(CHR < 23)

# plot ----
archipelago_plot <- function(df1, df2, 
                             plot_title = "Archipelago Plot", 
                             add_title = FALSE,
                             plot_subtitle = "Variant Set Association Test\nwith individual variant contributions\nand contributions for significant VSAT",
                             add_subtitle = FALSE,
                             chr_ticks = TRUE, 
                             show_legend = TRUE,
                             color_theme = NULL,
                             custom_colors = NULL,
                             color_labels = c("Chr: Individual p-val", "Chr: Individual p-val" ,"Individual variant\nfrom enriched VSAT","VSAT p-val"),
                             crit_val_VSAT = NULL,
                             crit_val_single_variant = NULL,
                             point_size = 1,
                             point_size_large = 3,
                             output_path_pdf = "../../data/archipelago/archipelago_plot.pdf",
                             output_path_png = "../../data/archipelago/archipelago_plot.png",
                             output_path_jpg = "../../data/archipelago/archipelago_plot.jpg",
                             output_raw_pdf = "../../data/archipelago/vsat_raw_plot.pdf",
                             output_raw_png = "../../data/archipelago/vsat_raw_plot.png",
                             output_raw_jpg = "../../data/archipelago/vsat_raw_plot.jpg",
                             output_patch_jpg = "../../data/archipelago/archipelago_raw_patch.jpg"
                             ) {

library(ggplot2)
library(dplyr)

# threshold for VSAT
if(is.null(crit_val_VSAT)) {
  set_ID_max <- df1$set_ID %>% unique() %>% length()
  crit_val_VSAT <- .05/set_ID_max
}

# threshold for single-variant test
if(is.null(crit_val_single_variant)) {
  crit_val_single_variant <- .05/nrow(df2)
}



df2$P <- -log10(df2$P)
df2 <- df2[order(df2$CHR, df2$BP), ]
df2$index=NA
ind = 0
for (i in unique(df2$CHR)){
  ind = ind + 1
  df2[df2$CHR==i,]$index = ind
}
df2$index = rep.int(seq_along(unique(df2$CHR)), times = tapply(df2$BP,df2$CHR,length))  

nchr = length(unique(df2$CHR))
df2$pos=NA

# The following code block replaces my old method - to more closely match qqman - to make it easier to integrate in future.
lastbase=0
ticks=NULL
for (i in unique(df2$index)) {
  if (i==1) {
    df2[df2$index==i, ]$pos=df2[df2$index==i, ]$BP
  } else {
    lastbase = lastbase +max(df2[df2$index==(i-1),"BP"])  
    df2[df2$index == i,"BP"] = df2[df2$index == i,"BP"]-min(df2[df2$index==i,"BP"]) +1
    df2[df2$index == i, "pos"] = df2[df2$index == i,"BP"] + lastbase 
    
  }
}
ticks <-tapply(df2$pos,df2$index,quantile,probs=0.5) 
xlabel = 'Chromosome'
labs <- unique(df2$CHR)

# Merging df1 and df2
merged_df <- bind_rows(df1, df2)

# Add a new column to specify variant_set points
merged_df$metric <- ifelse(is.na(merged_df$BP), "variant_set", "variant")

# Calculate the total of locations per set_ID
pos_sum <- aggregate(merged_df$pos, by=list(merged_df$set_ID), FUN=sum, na.rm=TRUE)

# Rank set_ID by the sum of locations
pos_sum$rank <- rank(pos_sum$x) 

# Calculate the average location within each set_ID group from the original dataframe
average_pos <- aggregate(merged_df$pos, by=list(merged_df$set_ID), FUN=mean, na.rm=TRUE)

# Replace NA locations with the average location of their set_ID group
merged_df$pos[is.na(merged_df$pos)] <- average_pos$x[match(merged_df$set_ID[is.na(merged_df$pos)], average_pos$Group.1)]

# Create a new grouping variable
merged_df$grouping <- with(merged_df, ave(pos, set_ID, FUN = function(x) cumsum(!is.na(x))))

# Get a simple ordered distribution of variant_set p value positions for the x-axis
# sort the "variant_set" group by "pos"
variant_set_sorted <- merged_df %>%
  filter(metric == "variant_set") %>%
  arrange(pos)

# test - if we have a subset of data, NA value may cause failure:
# Set default value as a random number between 1 and 22
if(anyNA(merged_df$pos)) {
	# Handle NA values
	# e.g. replace each NA with a different random number between 1 and 22
	merged_df$pos[is.na(merged_df$pos)] <- sample(1:22, sum(is.na(merged_df$pos)), replace = TRUE)
}

# create a sequence of evenly spaced numbers across the range of all "pos" values
even_pos <- seq(from = min(merged_df$pos),
                to = max(merged_df$pos),
                length.out = nrow(variant_set_sorted))

# assign these evenly spaced numbers to the "pos" column of the "variant_set" group
variant_set_sorted$pos <- even_pos

# replace the "pos" values in the original dataframe with these new evenly spaced numbers
merged_df <- merged_df %>%
  mutate(pos = ifelse(metric == "variant_set",
                      variant_set_sorted$pos[match(set_ID, variant_set_sorted$set_ID)],
                      pos))

# Add a new column for color groups
merged_df$color_group <- ifelse(merged_df$metric == "variant_set", "variant_set", 
                                ifelse(merged_df$CHR %% 2 == 0, "Chr_B", "Chr_A"))

# Create a dataframe with the mid-point of each chromosome
chromosome_ticks <- merged_df %>%
  filter(metric == "variant") %>%
  group_by(CHR) %>%
  summarise(mid_point = mean(pos))

# Split the data into variant_set and variant datasets
variant_set_data <- merged_df %>% filter(metric == "variant_set")
variant_data <- merged_df %>% filter(metric == "variant")

# Join the datasets to create a line dataset
line_data <- left_join(variant_data, variant_set_data, by = "set_ID", suffix = c("_variant", "_variant_set"))

# Raw plot ----
raw <- merged_df %>% 
  filter(metric == "variant_set")
raw$rank <- rank(raw$P)

p_raw <- 
  raw %>%
  ggplot() +
  geom_point(aes(x = rank, y = P), color='#27afea', size = point_size) +
  ylab("-log10 (p-value)") +
  geom_hline(linetype="dotted", 
             # color = "#135775",
             yintercept=-log10(crit_val_VSAT)) +
  theme_bw() +
  ggtitle("Raw VSAT p-value") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  guides(color = "none") 

# p_raw

ggsave(p_raw, filename = output_raw_pdf, width = 6, height = 4)
ggsave(p_raw, filename = output_raw_png, width = 6, height = 4)
ggsave(p_raw, filename = output_raw_jpg, width = 6, height = 4)

# Plot 2 ----
# Highlight individual variant contributions 
# Create a new variable 'color_condition' that checks the condition

# Step 1: Identify all set_ID groups where the metric is "variant_set" and P > -log10(crit_val_VSAT)
set_IDs_to_color <- merged_df %>%
	filter(metric == "variant_set", P > -log10(crit_val_VSAT)) %>%
	pull(set_ID)

# Step 2: Assign "condition_met" to "variant" points within those groups
merged_df <- 
	merged_df %>%
	mutate(color_condition = 
			 	ifelse(metric == "variant" & set_ID %in% set_IDs_to_color, 
			 			 "condition_met", 
			 			 color_group))


# merged_df <- 
#   merged_df %>%
#   group_by(set_ID) %>%
#   mutate(color_condition = 
#            ifelse(any(P > -log10(crit_val_VSAT)) & metric == "variant", 
#                   "condition_met", 
#                   color_group))

# Change the alpha variable accordingly
merged_df$alpha <- ifelse(merged_df$color_condition == "condition_met", 1, ifelse(merged_df$metric == "variant_set", 1, 0.5))

# Plot 4 ----
# Clearer condition_met layer
# Separate the points and lines that meet the condition
condition_met_points <- merged_df[merged_df$color_condition == "condition_met", ]
condition_met_lines <- line_data[line_data$set_ID %in% condition_met_points$set_ID, ]

# Figure  ----
# Define color themes
color_themes <- list(
  retro = c("#ffe28a", "#6fcb9f", "#666547", "#fb2e01"),
  metro = c("#f37735", "#ffc425", "#d11141", "#00aedb"),
  summer = c("#ddf098", "#f9d62e", "#ff4e50", "#fc913a"),
  messenger = c("#44bec7", "#ffc300", "#0084ff", "#fa3c4c"),
  sunset = c("#eeaf61", "#fb9062", "#6a0d83", "#ee5d6c"),
  alice = c("#f7c297", "#ffecb8", "#90d2d8", "#f6a6b2"),
  buckley = c("#9abfd8", "#cac1f3", "#371c4b", "#2a5b7f"), 
  romance = c("#ffcad4", "#c08497", "#3a4440", "#cc8562"),
  meme = c("#f4d4b8", "#a9dada", "#474747", "#f17255"),
  saiko = c("#a2b3d8", "#6289b3", "#d0444a", "#005e90"),
  pagliacci = c("#ffdd75", "#59a4d2", "#f7955d", "#49518a"),
  ambush = c("#f5e9be", "#ffe184", "#174c4f", "#207178"),
  sunra = c("#34b1ff", "#2d93d6", "#eaa221", "#1a406f"),
  caliber = c("#cccccc", "#bbbbbb", "#900303", "#000000"),
  yawn = c("#bcbcbc", "#999999", "#000000", "#5b5b5b"),
  lawless = c('#f6d992', '#f6a192', 'black', '#27afea')
)

# Check if custom_colors is provided, otherwise use color_theme or default
if (!missing(custom_colors)) {
  colors <- custom_colors
} else if (!missing(color_theme) && color_theme %in% names(color_themes)) {
  colors <- color_themes[[color_theme]]
} else {
  colors <- color_themes[['lawless']] # default
}


# Create a data frame with only the "variant_set" points
variant_set_points <- merged_df[merged_df$metric == "variant_set", ]

p_arch_leg <- 
	ggplot() +
	# Adding points where color_condition is not "condition_met"
	geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
				  aes(x = pos, y = P, color = color_condition), size = point_size, alpha = 1) +
	# Adding lines
	geom_segment(data = condition_met_lines, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.5) +
	# Adding points where color_condition is "condition_met"
	geom_point(data = condition_met_points, 
				  aes(x = pos, y = P, color = color_condition), size = point_size, alpha = 1) +
	# Adding points where metric is "variant_set", always on top
	geom_point(data = variant_set_points, 
				  aes(x = pos, y = P, color = color_condition), size = point_size_large, alpha = 1) +
  scale_color_manual(values = colors, 
                     labels= color_labels) +
  ylab("-log10 (p-value)") +
  theme_bw() +
  labs(colour = "P value") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10, "pt")) 

# Add the title only if add_title is TRUE
if (add_title) {
  p_arch_leg <- p_arch_leg + ggtitle(plot_title)
}

# Add the subtitle only if add_subtitle is TRUE
if (add_subtitle) {
  p_arch_leg <- p_arch_leg + labs(subtitle = plot_subtitle)
}

# Only show chromosome ticks if chr_ticks is TRUE
if (chr_ticks) {
  p_arch_leg <- p_arch_leg + scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR) + xlab("Chromosome") 
}

# Only show the legend if show_legend is TRUE
if (!show_legend) {
  p_arch_leg <- p_arch_leg + guides(color = "none") 
}

# Set critical p-val line
if (!is.null(crit_val_VSAT)) {
  p_arch_leg <- p_arch_leg + 
  geom_hline(linetype="dotted", color = "#135775",
             yintercept=-log10(crit_val_VSAT)) +
    geom_hline(linetype="dotted", color = "#f6a192",
               yintercept=-log10(crit_val_single_variant))
    
}

# label each p-val line
if (!is.null(crit_val_VSAT)) {
  p_arch_leg <- p_arch_leg + 
    geom_text(aes(x = max(merged_df$pos), y = (-log10(crit_val_VSAT)+0.3), 
                  label = "VSAT\nthreshold"), hjust = 1) +
    geom_text(aes(x = max(merged_df$pos), y = (-log10(crit_val_single_variant)+0.3), 
                  label = "single-variant\nthreshold"), hjust = 1)
}


p_arch_leg
# ggsave(p_arch_leg, filename = output_path, width = 8, height = 4)
ggsave(p_arch_leg, filename = output_path_pdf, width = 14, height = 6)
ggsave(p_arch_leg, filename = output_path_png, width = 14, height = 6)
ggsave(p_arch_leg, filename = output_path_jpg, width = 14, height = 6)

library(patchwork)
# patch1 <- (p_arch_leg / p_raw)
patch1 <- (p_arch_leg | p_raw) + plot_annotation(tag_levels = 'A')

ggsave(patch1, filename = output_patch_jpg, width = 14, height = 9)

# return the final plot
return(p_arch_leg)
}

# archipelago_plot(df1, df2)
archipelago_plot(df1, df2_ready)

 # 13
# 44949840
# A
# 0.0000000
# 0.021920
# T
# 24.190
# 8.727e-07
# 0.00000
# 8.727e-07