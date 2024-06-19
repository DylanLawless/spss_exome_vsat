source("ACMGuru_post_ppi_vcurrent.R")
# hold <- df_report
# df_report <- hold
df_report <- df

# NB !!!! VERSION SETTINGS:
# v4 10, 1e-3 (.001), HIGH-MODERATE
# df_report <- df_report |> filter(gnomAD_AF <= 1e-3)

rm(list=setdiff(ls(), c("df_report", "df", "geneset_MCL_ID", "file_suffix", "file_suffix")))
# file_suffix <- "gene_illustrate_"

# df_report_position <- df_report_position |> filter(ACMG_score > 0 )

library(dplyr)
library(ggplot2)
library(plotly)

# Notes ----
# Once the loop works. 
# Save each individual plot.
# Then use df_uniprot_meta table to hyperlink to page with plot + evidence table.

# import ----
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rtracklayer")
library(rtracklayer)

# Identdf_reportify canonical transcript ----
# Data sources ----
# Ensembl: https://www.ensembl.org/Homo_sapiens/
# Uniprot under format tab: https://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22# 
# "Homo sapiens (Human) [9606]"
# Consists of 204,185 entries

# Data import ----
# start a for loop here for all genes
df_uniprot <- readGFF("../../data/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.gff")
# df_uniprot <- readGFF("../data/uniprot_HUMAN_Q2TBE0_CWF19L2.gff")
df_uniprot_meta <- read.csv("../../data/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.tab", sep="\t")

# rename header to match
# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
colnames(df_uniprot_meta)[colnames(df_uniprot_meta) == 'Entry'] <- 'seqid'

df_uniprot$seqid
df_uniprot_meta$seqid
df_uniprot_meta$Gene.names

# Separate the Gene.names column where names are separated by space
df_uniprot_meta_tidy <- df_uniprot_meta %>% 
	separate_rows(Gene.names, sep = " ")

# Rename Gene.names to SYMBOL
df_uniprot_meta_tidy <- df_uniprot_meta_tidy %>% rename(SYMBOL = Gene.names)

# Create list 
# df_uniprot_meta <- 
  # df_uniprot_meta %>% 
  # filter(grepl("^reviewed$", Status))

# Uniprot plot ----
library(wesanderson)
library(stringr)

# Define a Wes style block of 25 colors
color_block <- c(
	"#00A08A", "#F2AD00", "#F98400", "#5BBCD6", # Darjeeling1 (removed "#FF0000" hard red)
	"#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", # Darjeeling2 (removed #"#000000" black)
	"#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20", #  FantasticFox1
	"#F1BB7B", "#FD6467", "#5B1A18", "#D67236", # GrandBudapest1
	"#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4", # GrandBudapest2
	"#446455", "#FDD262", "#D3DDDC", "#C7B19C" #Chevalier1
)

# Repeat the color block 7 times
wes_pal <- rep(color_block, times = 7)

# library(htmlwidgets) # saveWidget
# library(htmltools) # save_html
# data tidy ----
# if column is all NA, drop. This DT method is more complex than filter, but will be fast on large data.
library(data.table)
dt <- as.data.table(df_uniprot)

# 2 columns of NA. I will keep these for now while testing a bug
# dt_uniprot <- dt[,which(unlist(lapply(dt, function(x)!all(is.na(x))))),with=F]
dt_uniprot <- dt
rm(dt)

# add notes as shapes
dt_uniprot$Note <- as.character(dt_uniprot$Note)

types <- dt_uniprot %>% select(type)
types <- types[order(types$type), ] # sort
types <- types[!duplicated(types), ] # uniq

# test subset ----
# temp <- dt_uniprot
# dt_uniprot <- dt_uniprot[1:100,]
# dt_uniprot <- temp

# info ----
df_report |> 
	group_by(SYMBOL, #seqid, 
				Uniprot_acc, Uniprot_entry) |>
	summarise(v=n()) |> 
	arrange(desc(v))

# fill uniprot seqid ----
# get first entry for Uniprot_acc ID (delim = &)
df_report$seqid <- sapply(strsplit(as.character(df_report$Uniprot_acc), "&"), "[", 1)

df_report <- df_report %>%
	group_by(SYMBOL) %>%
	fill(seqid, .direction = "downup") |>
	ungroup()

# This might fill empty seqid, but we may have NA still. 
# We could fill all Uniprot data for each SYMBOL such that a match based on SYMBOL instead of seqid can be made.

# Display the dataframe
df_report |> select(SYMBOL, seqid)
df_report_missing <- df_report |> filter(is.na(seqid)) |> dplyr::select(SYMBOL) |> unique()
df_report_missing

# fill uniprot SYMBOL ----
# find largest uniprot entry per SYMBOL


# 1. Identify SYMBOLs present in both df_uniprot_meta_tidy_short and df_report
df_uniprot_meta_tidy_short <- df_uniprot_meta_tidy |> 
	dplyr::select(SYMBOL, seqid)

common_SYMBOLs <- intersect(df_uniprot_meta_tidy_short$SYMBOL, df_report$SYMBOL)

# Get seqids that correspond to the common_SYMBOLs
seqids_for_common_SYMBOLs <- df_uniprot_meta_tidy_short %>% 
	filter(SYMBOL %in% common_SYMBOLs) %>%
	pull(seqid)

rm(common_SYMBOLs)

# 2. Filter dt_uniprot for these seqid
dt_uniprot_filtered <- dt_uniprot %>% filter(seqid %in% seqids_for_common_SYMBOLs)

rm(seqids_for_common_SYMBOLs)

# 3. Merge df_uniprot_meta_tidy_short with dt_uniprot
merged_df <- left_join(dt_uniprot_filtered, df_uniprot_meta_tidy_short, by = "seqid")

rm(dt_uniprot_filtered)

# 4.. Group by SYMBOL and seqid, then count the number of rows for each combination
grouped_df <- merged_df %>%
	group_by(SYMBOL, seqid) %>%
	summarise(n = n(), .groups = "drop")

rm(merged_df)

# 5. For each SYMBOL, keep only the rows corresponding to the seqid with the maximum count
grouped_df_max <- grouped_df %>%
	group_by(SYMBOL) %>%
	slice_max(n) |> 
	filter(!SYMBOL=="")  |>
	dplyr::select(-n)

rm(grouped_df)

class(grouped_df_max$SYMBOL)
class(grouped_df_max$seqid)

# 6. overwrite seqid fully
df_report$seqid_old <- df_report$seqid

df_report <- left_join(df_report, grouped_df_max, by = "SYMBOL", suffix = c("", ".y")) |>
	mutate(seqid = seqid.y)|>
	select(-seqid.y) |>
	dplyr::select(SYMBOL, seqid, Protein_position, CDS_position, rownames, everything())

# get first entry for Uniprot_acc ID (delim = &)
# df_report$seqid <- sapply(strsplit(as.character(df_report$Uniprot_acc), "&"), "[", 1)

# subset
unique_uniprot <- unique(df_report$seqid)
filtered_dt_uniprot <- dt_uniprot[dt_uniprot$seqid %in% unique_uniprot, ]

# get variant protein position ----
df_report_position <- df_report |> dplyr::select(SYMBOL, Protein_position, CDS_position, rownames, seqid, ACMG_score) |> unique()

df_report_position$ACMG_score <- as.numeric(df_report_position$ACMG_score)
df_report_position$Protein_position <- as.numeric(df_report_position$Protein_position)
df_report_position$Protein_position
# Remove NA - maybe use NA to detect non-coding variants like splice site
df_report_position <- df_report_position[!is.na(df_report_position$Protein_position), ]

# add SYMBOL to uniprot ----
df_report_seqid_SYMBOL <- df_report |> dplyr::select(SYMBOL, seqid) |> unique()
filtered_dt_uniprot <- merge(df_report_seqid_SYMBOL, filtered_dt_uniprot )
filtered_dt_uniprot_hold <- filtered_dt_uniprot
filtered_dt_uniprot <- filtered_dt_uniprot_hold

# test ----
# PDF single evidence plot ----

# plot ----
create_vlines <- function(data) {
	vlines <- list()
	vlines[[1]] <- geom_vline(data = data %>% filter(ACMG_score == 0), aes(xintercept = Protein_position), color = "#84A98C"
									  )
	vlines[[2]] <- geom_vline(data = data %>% filter(ACMG_score == 1), aes(xintercept = Protein_position),color = "#52796F")
	vlines[[3]] <- geom_vline(data = data %>% filter(ACMG_score == 2), aes(xintercept = Protein_position),color = "#354F52")
	vlines[[4]] <- geom_vline(data = data %>% filter(ACMG_score == 4), aes(xintercept = Protein_position),color = "#2F3E46")
	vlines[[5]] <- geom_vline(data = data %>% filter(ACMG_score == 8), aes(xintercept = Protein_position),color = "black")
	return(vlines)
}

# create_and_save_plots <- function(plot_list, file_suffix, target_SYMBOL, filename) {
# create_plot <- function(plot_list, file_suffix, target_SYMBOL, filename) {
create_plot <- function(filter_acmg_score, target_SYMBOL = NULL) {    
# to print just one chosen plot, set the "target_SYMBOL"
  if (!is.null(target_SYMBOL)) {
    filtered_dt_uniprot <- filtered_dt_uniprot %>% filter(SYMBOL == target_SYMBOL)
  }
  
# create_plot <- function(filter_acmg_score) {
	plot_list <-
		lapply(split(filtered_dt_uniprot, filtered_dt_uniprot$seqid, drop = TRUE),
				 function(x_dt_uniprot) {
				 	current_seqid <- unique(x_dt_uniprot$seqid)
				 	
				 	current_SYMBOL <- unique(x_dt_uniprot$SYMBOL)
				 	
				 	current_protein_positions <- df_report_position |>
				 		filter(seqid == current_seqid) |>
				 		dplyr::select(Protein_position, ACMG_score) |>
				 		unique()
				 	
				 	if(filter_acmg_score) {
				 		current_protein_positions <- current_protein_positions |> 
				 			filter(ACMG_score > 0 )
				 		# If current_protein_positions is empty, skip this iteration since we have no reason to plot the gene
				 		if(nrow(current_protein_positions) == 0) {
				 			return(NULL)
				 		}
				 	}
				 	
				 	# Filter seqid present in df_report_position$seqid
				 	current_seqid <- x_dt_uniprot$seqid[x_dt_uniprot$seqid %in% df_report_position$seqid] |> unique()
				 	
				 	
				 	# Filter seqid present in df_report_position$seqid
				 	current_seqid <- x_dt_uniprot$seqid[x_dt_uniprot$seqid %in% df_report_position$seqid] |> unique()
				 	
				 	# Convert ACMG_score to factor
				 	current_protein_positions$ACMG_score <- as.factor(current_protein_positions$ACMG_score)
				 	
				 	x_dt_uniprot$position_label <-
				 		rowMeans(x_dt_uniprot[, c('start', 'end')], na.rm = TRUE)
				 	
				 	# Ensuring there is a visible difference between start and end if they're identical
				 	x_dt_uniprot$end[x_dt_uniprot$end == x_dt_uniprot$start] <-
				 		x_dt_uniprot$end + 1
				 	
			 	Domain <- c("Chain", "Domain", "Region", "Motif")
			 	Structure <- c("Helix", "Turn", "Beta strand")
			 	types <-
			 		setdiff(unique(x_dt_uniprot$type), c(Domain, Structure))
			 	
			 	# This function maps 'type' to a 'label' according to the specified conditions
			 	x_dt_uniprot <- x_dt_uniprot |>
			 		mutate(
			 			label = case_when(
			 				type %in% Domain ~ "Family & Domain",
			 				type %in% Structure ~ "Structure",
			 				type %in% types ~ "Features"
			 			)
			 		)
			 	
			 	# Ensuring uniqueness within the 'Family & Domain' category
			 	x_dt_uniprot |>
			 		group_by(label) |>
			 		mutate(type = if_else(label == "Family & Domain", 
			 									 make.unique(as.character(type), sep = "_"), 
			 									 type)) -> x_dt_uniprot
			 	
			 	
			 	# Replacing "character(0)" with ""
			 	x_dt_uniprot <- x_dt_uniprot |>
			 		mutate(across(c(Note, Dbxref), str_replace_all, "character\\(0\\)", "")) |>
			 		mutate(across(c(Note, Dbxref, evidence), str_wrap, width = 30))
			 	
			 	# The plot ----
			 	p <- x_dt_uniprot |>
			 		group_by(type, Note) |>
			 		ggplot(aes(
			 			y = type,
			 			x = start,
			 			label = end,
			 			label2 = Note,
			 			label3 = Dbxref,
			 			label4 = evidence
			 		)) +
			 		geom_segment(aes( x = start, xend = end, y = type, yend = type, color = type ), 
			 						 size = 4, show.legend = FALSE) +
			 		facet_grid(vars(label), scales = "free", space = "free") +
			 		geom_text( data = x_dt_uniprot |> filter(label == "Family & Domain"),
			 					  aes(label = Note, x = position_label, y = type), hjust = 0, vjust = 0, size = 2)  +
			 		ylab("") + xlab("") + labs(title = current_SYMBOL) +
			 		theme_bw() +
			 		theme(panel.background = element_blank()) +
			 		scale_color_manual(values = wes_pal) 
			 	
			 	vlines <- create_vlines(current_protein_positions)
			 	p <- p + vlines
			 	
				 })
	# Remove NULL elements from the list
	plot_list <- plot_list[!sapply(plot_list, is.null)]
	
	return(plot_list)
}

# This function calculates the layout dimensions based on the number of plots.
# For every 4 plots, it assigns a width of 24 and a height of 7.
# If the number of plots is less than 4, it adjusts the dimensions accordingly.
# Input: n_plots - number of plots
# Output: A list with the number of columns (ncol), and the dimensions of height and width
get_layout_dims <- function(n_plots) {
	# Handle the case when there are no plots to display
	if(n_plots == 0) {
		return(list(ncol = 1, height = 1, width = 1))
	}
	base_width_per_ncol <- 5 # Adjust this to your preference
	height_per_nrow <- 7
	s <- sqrt(n_plots)
	ncol <- if(s == floor(s)) s else ceiling(n_plots / floor(s))
	nrow <- if(s == floor(s)) s else floor(s)
	height <- nrow * height_per_nrow
	width <- ncol * base_width_per_ncol
	return(list(ncol = ncol, height = height, width = width))
}


# This function calculates the layout dimensions based on the number of plots.
# It tries to make the layout as square as possible (e.g., 2x2, 3x3, 4x4, etc.),
# but if the number of plots is not a perfect square, it adjusts the dimensions accordingly.
# If there are no plots, it returns a layout of 1x1.
# The base width for each column is 6, and the base height for each row is 7.
# Inputs: 
#    n_plots - number of plots
# Outputs: 
#    A list with the number of columns (ncol), and the dimensions of height and width
create_and_save_plots <- function(plot_list, file_suffix, filename) {
	n_plots <- length(plot_list)
	dims <- get_layout_dims(n_plots)
	ncol <- dims$ncol
	nrow <- ceiling(n_plots / ncol) 
	plot_list <- c(plot_list)
	
	# Arrange all plots together
	p.evidence_plots <-  
		annotate_figure(
			ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow),
			left = textGrob("Evidence type", rot = 90, vjust = 1 ),
			bottom = textGrob("Position" )
		)
	
	# Construct the filename with file_suffix
	filename <- paste(file_suffix, filename, sep = "")
	
	# Save the plots
	# ggsave(paste("../../data/ACMGuru_post_ppi/", filename, sep = ""),plot = p.evidence_plots, 
			 # height = dims$height, width = dims$width, limitsize = FALSE)

# Save the plots
ggsave(paste("../../images/ACMGuru_gene_illustrate_post_ppi/AMCGuru_", filename, sep = ""),plot = p.evidence_plots, 
       height = dims$height, width = dims$width, limitsize = FALSE)

}

# Then can call create_plot() with either TRUE or FALSE for filtering
plot_list_with_filter <- create_plot(TRUE)
plot_list_without_filter <- create_plot(FALSE)

# Save the plots
create_and_save_plots(plot_list_with_filter, file_suffix, "evidence_plots_with_filter.pdf")
create_and_save_plots(plot_list_without_filter, file_suffix, "evidence_plots_without_filter.pdf")

# individual plots ----
# Create list of SYMBOLs with ACMG_total_score >= 6
SYMBOL_list <- df_report %>%
  # filter(ACMG_total_score >= 6) %>%
  filter(ACMG_total_score >= 0) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL)

# Iterate over SYMBOL_list and make single individual plot for each
for (target_SYMBOL in SYMBOL_list) {
  # Create plots
  # plot_list_with_filter_target_SYMBOL <- create_plot(TRUE, target_SYMBOL)
  # plot_list_without_filter_target_SYMBOL <- create_plot(FALSE, target_SYMBOL)
  plot_list_with_filter_target_SYMBOL <- create_plot(TRUE)
  plot_list_without_filter_target_SYMBOL <- create_plot(FALSE)
  
  # Save plots
  create_and_save_plots(plot_list_with_filter_target_SYMBOL, file_suffix, paste("evidence_plots_with_filter_", target_SYMBOL, ".pdf"))
  create_and_save_plots(plot_list_without_filter_target_SYMBOL, file_suffix, paste("evidence_plots_without_filter_", target_SYMBOL, ".pdf"))
  
}

# save file for ACMGuru_gene_uniprotr
saveRDS(grouped_df_max, paste0("../../data/ACMGuru_post_ppi/acmguru_gene_illustrate_grouped_df_max", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

saveRDS(df_report, paste0("../../data/ACMGuru_post_ppi/acmguru_gene_illustrate_df_report", paste(geneset_MCL_ID, collapse="_"), ".Rds"))



# Save the plots with specified height and width
# create_and_save_plots(plot_list_with_filter, file_suffix, "evidence_plots_with_filter.pdf", height = 28, width = 24)
# create_and_save_plots(plot_list_without_filter, file_suffix, "evidence_plots_without_filter.pdf", height = 36, width = 24)


# END ----
# END ----
# END ----
