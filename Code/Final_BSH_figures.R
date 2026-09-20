######################################################################################################################################
######################################################################################################################################
#########All input files used in the script are provided in the GitHub linked in the Data Availability section########################
#########FBMN job link with NEW synthesis filtered library: https://gnps2.org/status?task=6af4031f410d4284b964b0d3faf5155b
######################################################################################################################################
######################################################################################################################################

#set working directory and load libraries
setwd("F:/OneDrive - University of California, San Diego Health/BSH_screening_April/Final_assay_BSH_Exploris")
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(tibble)
library(grid)
library(data.table)
library(stringr)
library(purrr)
library(RColorBrewer)

###read the metadata
df1 <- read.csv("F:/OneDrive - University of California, San Diego Health/BSH_screening_April/Samples_MS.csv", header = TRUE)
df2 <- read.csv("F:/OneDrive - University of California, San Diego Health/BSH_screening_April/Sample_metadata_ALL.csv", header = TRUE)
df_sample_meta <- merge(df1, df2, by = "code", all = FALSE)

#replace A0A317T5Z6 and A0A952U252 with the updated phylum code
df_sample_meta <- df_sample_meta %>%
  mutate(
    Phylum = case_when(
      code == "A0A317T5Z6" ~ "Bacteroidota",
      code == "A0A952U252" ~ "Cyanobacteriota",
      TRUE ~ Phylum
    )
  )

##############################################################################
##############################################################################

#read the final list of bile amidates
df_amidates <- fread("NEW_Stage2_BAs_amines_for_heatmap_manual.csv", header = TRUE)

df_amidates <- df_amidates %>%
  filter(!is.na(Code))

# Step 1: Identify rows where Code contains "Negative_ctrl" or "Only_sub"
control_rows <- df_amidates %>%
  filter(str_detect(Code, "Negative_ctrl|Only_sub"))

control_rows <- as.data.frame(control_rows)

# Step 2: Identify numeric columns from column 4 onward
cols_to_check <- names(df_amidates)[4:ncol(df_amidates)]

# Step 3: Find columns where any value exceeds 1000 in control rows
cols_exceeding <- cols_to_check[
  sapply(control_rows[, cols_to_check], function(col) any(as.numeric(col) > 1000, na.rm = TRUE))
]

# Step 4: Remove those columns from the original dataframe
df_amidates_filtered <- df_amidates %>%
  select(-all_of(cols_exceeding))


##plotting heatmap
merge2 <- merge(df_sample_meta, df_amidates_filtered, by.x = "code", by.y = "Code", all.y = TRUE)
merge3 <- merge2[,-c(2,3,6:9)]
merge3 <- merge3 %>%
  mutate(across(6:ncol(.), as.numeric))


#use merge3 for removing and arranging based on phylum tree from April 
merge3_df <- merge3

#remove Melainabacteria, Euryarchaeota, Not Assigned
merge3_filtered <- merge3_df %>%
  dplyr::filter(!Phylum %in% c("Candidatus Melainabacteria", "Chlorobiota") & !is.na(Phylum))

merge3_controls <- merge3_df %>%
  dplyr::filter(is.na(Phylum))

merge3_filtered$Phylum[merge3_filtered$Phylum == "Euryarchaeota"] <- "Methanobacteriota"

merge3_filtered <- merge3_filtered %>%
  filter(code != "BBD47700")

#combine the three replicates
merge3_collapsed <- merge3_filtered %>%
  group_by(code) %>%
  summarise(
    Phylum = first(Phylum),
    Genus = first(Genus),
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

merge3_collapsed <- merge3_collapsed %>%
  as.data.frame()

merge3_collapsed <- merge3_collapsed %>%
  mutate(code = ifelse(code == "Pencillin_amidase",
                       "B1C4T2",
                       code))

##Updated list of BSH sequences to include in the heatmap that were validated by proteomics
df_anna_list <- fread("BSH_final_seq_Anna_06102026.csv", header = TRUE)

#subset for only those IDs that are in proteomics validated list
merge3_collapsed <- merge3_collapsed[merge3_collapsed$code %in% df_anna_list$Final_IDs, ]

# Remove non-numeric columns and set row names
rownames(merge3_collapsed) <- merge3_collapsed$code  
merge3_numeric <- merge3_collapsed %>% select(-c(code, Phylum, Genus))
merge3_matrix <- as.matrix(merge3_numeric)

# Add epsilon (small value) to avoid log(0) issues
epsilon <- 1e-5  # Small constant
merge3_matrix_log <- log10(merge3_matrix + epsilon)

# Flatten the matrix into a vector
values <- as.vector(merge3_matrix_log)

# Filter only positive values (greater than zero)
positive_values <- values[values > 0]

# Define min and max values in the matrix
min_value <- min(merge3_matrix_log, na.rm = TRUE)#1.901
max_value <- max(merge3_matrix_log, na.rm = TRUE)#8.07

# Create column annotation for the phylum of each sample
annotation_col <- data.frame(
  Phylum = merge3_collapsed$Phylum[match(colnames(merge3_matrix_log), merge3_collapsed$code)] # Match genus to sample columns
)

annotation_row <- data.frame(
  Phylum = merge3_collapsed$Phylum[match(rownames(merge3_matrix_log), merge3_collapsed$code)]
)
rownames(annotation_row) <- rownames(merge3_matrix_log)
# Set row names for annotation to match column names of the heatmap
rownames(annotation_col) <- colnames(merge3_matrix_log)

##################################################################################
###add signal peptide length and the border to the heatmap based on CV############
##################################################################################

#to overlay the information about signal peptide
df_SP <-  fread("Signal_peptide_length.csv")

# keep only needed columns and rename
df_SP2 <- df_SP %>%
  dplyr::select(`Tree node label`, SP_length) %>%
  dplyr::rename(code = `Tree node label`) %>%
  dplyr::mutate(
    SP_length = as.numeric(SP_length),
    SP_present = ifelse(!is.na(SP_length) & SP_length > 0, "Yes", "No")
  )

merge3_stats <- merge3_filtered %>%
  group_by(code) %>%
  summarise(
    Phylum = first(Phylum),
    Genus  = first(Genus),
    across(
      where(is.numeric),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  ) %>%
  left_join(df_SP2, by = "code")

# ---- Build mean and SD matrices ----
mean_matrix <- merge3_stats %>%
  select(ends_with("_mean")) %>%
  as.matrix()

sd_matrix <- merge3_stats %>%
  select(ends_with("_sd")) %>%
  as.matrix()

rownames(mean_matrix) <- merge3_stats$code
rownames(sd_matrix)   <- merge3_stats$code

# Clean column names so mean and sd matrices match sample names
colnames(mean_matrix) <- sub("_mean$", "", colnames(mean_matrix))
colnames(sd_matrix)   <- sub("_sd$", "", colnames(sd_matrix))

# ---- Log-transform mean matrix ----
epsilon <- 1e-5
mean_matrix_log <- log10(mean_matrix + epsilon)


# ---- Row annotation ----
annotation_row4 <- data.frame(
  Phylum = merge3_stats$Phylum,
  SP_length = merge3_stats$SP_length
)
rownames(annotation_row4) <- merge3_stats$code

# ---- Annotation colors ----
annotation_colors4 <- list(
  Phylum = c(
    "Actinomycetota" = "#A11D21",
    "Bacillota" = "#EF3C2D",
    "Bacteroidota" = "#88439A",
    "Fusobacteriota" = "#757575",
    "Planctomycetota" = "#541953",
    "Pseudomonadota" = "#185D38",
    "Thermodesulfobacteriota" = "#C3902C",
    "Cyanobacteriota" = "#76C376",
    "Methanobacteriota" = "#0A303F"
  ),
  SP_length = colorRampPalette(c("#EEEEEE", "#6FCF97", "#1F6F5F"))(100)
)

# ---- CV matrix and border logic ----
cv_matrix <- sd_matrix / (mean_matrix + epsilon)

# Cells to highlight with bold border
border_cells <- (cv_matrix < 0.4) & (mean_matrix_log > 1)

# ---- Plot heatmap ----
b4 <- pheatmap(
  mat = mean_matrix_log,
  color = colorRampPalette(c("#FFFFFF", "#618DB5", "#C07B74"))(300),
  annotation_row = annotation_row4,
  annotation_colors = annotation_colors4,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cellwidth = 10,
  cellheight = 3.5,
  fontsize_row = 4,
  fontsize_col = 6,
  display_numbers = FALSE,
  border_color = "grey80",
  main = "Heatmap with Phylum, Signal Peptide, and CV-Based Borders"
)

# ---- Reorder border matrix to match plotted heatmap ----
if (!is.null(b4$tree_row)) {
  border_cells <- border_cells[b4$tree_row$order, , drop = FALSE]
}
if (!is.null(b4$tree_col)) {
  border_cells <- border_cells[, b4$tree_col$order, drop = FALSE]
}

# ---- Add thick borders to selected cells ----
gt <- b4$gtable
mat_i <- which(gt$layout$name == "matrix")
g <- gt$grobs[[mat_i]]

rect_id <- which(vapply(g$children, inherits, logical(1), what = "rect"))[1]
cell_rect <- g$children[[rect_id]]

nr <- nrow(border_cells)
nc <- ncol(border_cells)

border_idx <- which(border_cells, arr.ind = TRUE)

target_x <- (border_idx[, "col"] - 0.5) / nc
target_y <- 1 - (border_idx[, "row"] - 0.5) / nr

rx <- as.numeric(convertX(cell_rect$x, "npc", valueOnly = TRUE))
ry <- as.numeric(convertY(cell_rect$y, "npc", valueOnly = TRUE))

rect_indices <- integer(nrow(border_idx))
for (k in seq_len(nrow(border_idx))) {
  rect_indices[k] <- which.min((rx - target_x[k])^2 + (ry - target_y[k])^2)
}

N <- length(rx)
col_vec <- rep("grey80", N)
lwd_vec <- rep(0.4, N)

col_vec[rect_indices] <- "black"
lwd_vec[rect_indices] <- 1.8

cell_rect$gp <- gpar(
  col = col_vec,
  lwd = lwd_vec,
  fill = cell_rect$gp$fill
)

g$children[[rect_id]] <- cell_rect
gt$grobs[[mat_i]] <- g

grid.newpage()
grid.draw(gt)

pdf("202600627_bile_amides_heatmap_cv_signalpeptide_length.pdf", width = 14, height = 10)
grid.draw(gt)
dev.off()

############################################################################
############################################################################
#########Supplementary figure for Figure 3##################################
############################################################################
############################################################################
# Read data
df_BBAA <- fread("BBAAs_for_Figure.csv")

# Metadata columns that are NOT BBAAs
metadata_cols <- c("code", "Phylum", "Genus", "filename", "Replicate")

# BBAA columns
bbaa_cols <- setdiff(names(df_BBAA), metadata_cols)

# ------------------------------------------------------------
# 1. Determine whether each BBAA was detected for each enzyme
#    Detection = value > 0 in at least one replicate
# ------------------------------------------------------------

enzyme_detection <- df_BBAA[, 
                            lapply(.SD, function(x) any(x > 0, na.rm = TRUE)),
                            by = code,
                            .SDcols = bbaa_cols
]

# ------------------------------------------------------------
# 2. Calculate detection frequency across enzymes
# ------------------------------------------------------------

frequency <- data.table(
  BBAA = bbaa_cols,
  n_detected = sapply(enzyme_detection[, ..bbaa_cols], sum)
)

frequency[, percent_detection := n_detected / uniqueN(df_BBAA$code) * 100]

# ------------------------------------------------------------
# 3. Select top 10
# ------------------------------------------------------------

top10 <- frequency[
  order(-percent_detection)
][1:10]

top10

# Order BBAAs by frequency
top10[, BBAA := factor(
  BBAA,
  levels = rev(BBAA)
)]

plot <- ggplot(top10, aes(x = percent_detection, y = BBAA)) +
  
  # Bar color
  geom_col(
    width = 0.7,
    fill = "#678BB0"
  ) +
  
  # Percentage labels
  geom_text(
    aes(label = sprintf("%.1f%%", percent_detection)),
    hjust = -0.15,
    size = 4
  ) +
  
  # Axis labels
  labs(
    x = "Detection frequency (%)",
    y = "BBAA"
  ) +
  
  # X-axis limits
  scale_x_continuous(
    limits = c(0, max(top10$percent_detection) * 1.15),
    expand = c(0, 0)
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    # Show BOTH x and y axis lines
    axis.line.x = element_line(
      color = "black",
      linewidth = 0.7
    ),
    axis.line.y = element_line(
      color = "black",
      linewidth = 0.7
    ),
    
    # Axis ticks
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.7
    ),
    
    # Text formatting
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 13,
      color = "black"
    )
  )

ggsave(
  "Top10_BBAA_detection_frequency.png",
  plot = plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  "Top10_BBAA_detection_frequency.svg",
  plot = plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
