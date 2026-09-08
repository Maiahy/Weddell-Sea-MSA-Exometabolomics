####Data cleaning and filtering for multivariate analysis####

#This code was written based on the recommendations provided by Broadhurst and al. 
#(DOI: https://doi.org/10.1007/s11306-018-1367-3). Please cite this article if you are using this code.
#And with the help of this website:
#https://www.davidzeleny.net/anadat-r/doku.php/en:pcoa_nmds

#Loading of the packages:
library(ggplot2)
library(tidyverse)
library(readxl)
library(tidyr)
library(janitor)
library(dplyr)
library(stringr)

library(here)
base_path <- here()
base_path

#####1) Data import ####


#feature table
feature_table <- read_tsv("SKYLINE_Peak_areas.txt")
colnames(feature_table)
feature_table$Replicate_Name


#metadata of the samples
metadata_samples <- read.table("Metadata_final_with_operators.txt", sep = "\t") %>% 
  row_to_names(row_number = 1)
colnames(metadata_samples)
unique(metadata_samples$sample_type)

#Creation of a short ID
feature_table <- feature_table %>%
  mutate(sample_id = Replicate_Name)
unique(feature_table$sample_id)
colnames(feature_table)
metadata_samples <- metadata_samples %>%
  mutate(
    sample_id = str_extract(
      file_name,
      "Blank_.*_R\\d+_R\\d+|Blank_SAMEA\\d+_.*_R\\d+|SCQC_.*_R\\d+|Sample_.*_R\\d+"
    )
  )

#Merging
merged_data <- left_join(feature_table,metadata_samples, by = "sample_id")
unique(merged_data$sample_type)
colnames(merged_data)
rm(feature_table)
rm(metadata_samples)

#####2) Exploration of the samples intensities #####
colnames(merged_data)
library(ggplot2)
library(dplyr)
library(stringr)

#Filtering one molecule
samples_only <- merged_data %>%
  filter(sample_type == "Sample") %>%
  filter(Molecule == "Gonyol")  # You can re-run the code with different molecules!


#Cleaning of the names befor plotting
samples_only <- samples_only %>%
  mutate(
    sample_label = str_replace(Replicate_Name, "Sample_(SAMEA\\d+)_R\\d+_R(\\d+)", "\\1 R\\2"),
    sample_label = str_replace_all(sample_label, "_", " "),
    sample_label = str_to_title(sample_label)
  )

#From little intensities to high
samples_only <- samples_only %>%
  arrange(Area) %>%
  mutate(sample_label = factor(sample_label, levels = unique(sample_label)))
samples_only$Area <- as.numeric(samples_only$Area) 


#Plot
ggplot(samples_only, aes(x = sample_label, y = Area)) +
  geom_col(fill = "turquoise4", width = 0.8) +  
  labs(
    x = "Sample ID",
    y = "Peak area"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.box = "vertical"
  )

#PDF save
ggsave(
  filename = "1_SAMPLES_plot_intensities_Gonyol.pdf",
  device = cairo_pdf,
  width = 18, height = 6, units = "in"
)

#TIFF 
ggsave(
  filename = "1_SAMPLES_plot_intensities_Gonyol.tiff",
  device = "tiff",
  dpi = 600,                 # haute résolution
  width = 18, height = 6, units = "in"
)


#Check the ratio between the samples
library(dplyr)
library(ggplot2)
library(stringr)

samples_only$Area <- as.numeric(samples_only$Area)

#Matrix with all ratios
ratio_matrix <- outer(
  samples_only$Area,
  samples_only$Area,
  function(x, y) pmax(x, y) / pmin(x, y)
)
pairs <- which(upper.tri(ratio_matrix), arr.ind = TRUE)

#Results
ratio_table <- data.frame(
  Sample1 = samples_only$sample_label[pairs[,1]],
  Sample2 = samples_only$sample_label[pairs[,2]],
  Fold_change = ratio_matrix[upper.tri(ratio_matrix)]
)

#From little to large
ratio_table <- ratio_table[order(ratio_table$Fold_change), ]

#All info
summary(ratio_table$Fold_change)

quantile(ratio_table$Fold_change,
         probs = c(0.25, 0.5, 0.75, 0.9, 0.95),
         na.rm = T)
min(samples_only$Area)
max(samples_only$Area)

##### 3) QC intensities  #####

# <!> Attention, this code has to be repeated for each single molecule!

colnames(merged_data)
unique(merged_data$Molecule)

#Filtering
unique(merged_data$sample_type)
qc_only <- merged_data %>%
  filter(Molecule == "MSA") %>%
  filter(sample_type == "SCQC") %>%
  filter(Replicate_Name != "SCQC_SCQC_R01_R1") %>%
  mutate(
    replicate_number = str_extract(Replicate_Name, "_R\\d+$"),           
    replicate_number = as.numeric(str_remove(replicate_number, "_R"))   
  ) %>%  
  mutate(
    sample_label = paste0("SCQC", replicate_number)
  )

#Injection order
qc_only <- qc_only %>%
  arrange(replicate_number) %>%
  mutate(sample_label = factor(sample_label, levels = unique(sample_label)))
qc_only$Area <- as.numeric(qc_only$Area)

#Minum and maximum signal
min_area <- min(qc_only$Area, na.rm = TRUE)
max_area <- max(qc_only$Area, na.rm = TRUE)

#Ratio max/min
signal_ratio <- max_area / min_area

cat("Minimum peak area =", signif(min_area, 4), "\n")
cat("Maximum peak area =", signif(max_area, 4), "\n")
cat("Max/Min ratio =", round(signal_ratio, 2), "x\n")

mean_area <- mean(qc_only$Area, na.rm = TRUE)
sd_area <- sd(qc_only$Area, na.rm = TRUE)
rsd_percent <- (sd_area / mean_area) * 100

#Show the results
cat("Relative Standard Deviation (RSD) in the QCs =", round(rsd_percent, 2), "%\n")



median_area <- median(qc_only$Area, na.rm = TRUE)
mad_area <- mad(qc_only$Area, na.rm = TRUE)

rmad_percent <- (mad_area / median_area) * 100


cat("Relative Median Absolute Deviation (RMAD) in the QCs =",
    round(rmad_percent, 2), "%\n")

#Plot
qc_plot <- ggplot(qc_only, aes(x = replicate_number, y = Area)) +
  geom_col(fill = "turquoise4", width = 0.8) +
  labs(
    title = "DMSOP",
    x = "QC replicate",
    y = "Peak area"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("RSD = ", round(rmad_percent, 2), "%"),
    hjust = 1.1, vjust = 1.5,
    size = 4.5,
    fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )


print(qc_plot)

#PDF
ggsave(
  filename = "1_QC_plot_DMSOP_intensities_ordered.pdf",
  plot = qc_plot,
  device = "pdf",
  width = 10, height = 6, units = "in"
)

#TIFF
ggsave(
  filename = "1_QC_plot_DMSOP_intensities_ordered.tiff",
  plot = qc_plot,
  device = "tiff",
  dpi = 600,
  width = 10, height = 6, units = "in"
)

##### 4) Volume normalization ####

#A normalization based n the volume can be performed in 2 steps:
#1) Calculatation of a dilution factore with DF = V sample/ Volume reference. Here we consider the reference volume to be 2L
#2) Normalization of the intensities based on the dilution factor: divide the intensity of the metabolite by the dilution factor.
#Here we have already calculated DF in the excel file we then just have to divide the entire table by the values of DF contained in the metadata.

merged_data$Area <- as.numeric(merged_data$Area)
filtered_data <- merged_data %>%
  mutate(
    Area_per_L = ifelse(
      is.na(Volume_normalization_factor),
      Area,  # si NA : ne pas normaliser, garder l'aire corrigée telle quelle
      Area / as.numeric(Volume_normalization_factor)  # sinon : normaliser
    )
  )

colnames(filtered_data)

write_csv(filtered_data, "1_Table_peak_area_after_volume_normalization.csv")

minidata <- select(
  filtered_data,
  Replicate_Name,
  Area,
  Volume_normalization_factor,
  Area_per_L
)

library(writexl)

write_xlsx(minidata, "minidata_check_volume_normalization.xlsx")


#####5) Blanks ####

######5.1) Exploration #####

library(dplyr)
library(ggplot2)


plot_data <- filtered_data %>%
  filter(
    sample_type %in% c("Sample", "Blank"),
    Molecule == "DMSP"
  ) %>%
  mutate(Area_per_L = as.numeric(Area_per_L))

ggplot(plot_data,
  aes(x = ESSENTIAL_sample_id,
    y = Area_per_L,
    color = sample_type)) +
  geom_point(
    position = position_jitter(width = 0.15),
    size = 2,
    alpha = 0.7) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5) +
  scale_color_manual(
    values = c(
      "Sample" = "black",
      "Blank" = "red")) +
  labs(
    x = "Sample ID",
    y = "Normalized peak area",
    title = "DMSP",
    color = "Sample type") +
   theme_minimal() +
  theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5))

colnames(plot_data)
ggplot(plot_data,
       aes(x = sample_type,
           y = Area_per_L,
           color = ESSENTIAL_sample_id)) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.7,
    size = 2
  ) +
  
  scale_color_manual(
    values = c(
      "SAME7905926" = "red",
      "Sample" = "black"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Normalized peak area",
    title = "DMSA"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none"
  )


###### 5.2) Assessment of global blank contribution ####

# We first average the technical triplicates for each biological sample
# and for the procedural blank.

replicate_means <- filtered_data %>%
  filter(
    sample_type == "Sample" |
      ESSENTIAL_sample_id == "SAME7905926"
  ) %>%
  mutate(
    Area_per_L = as.numeric(Area_per_L)
  ) %>%
  group_by(
    Molecule,
    ESSENTIAL_sample_id,
    sample_type
  ) %>%
  summarise(
    mean_area = mean(Area_per_L, na.rm = TRUE),
    .groups = "drop"
  )

replicate_means


sample_level <- filtered_data %>%
  filter(
    sample_type == "Sample" |
      ESSENTIAL_sample_id == "SAME7905926"
  ) %>%
  group_by(
    Molecule,
    ESSENTIAL_sample_id,
    sample_type
  ) %>%
  summarise(
    mean_area = mean(Area_per_L, na.rm = TRUE),
    .groups = "drop"
  )

blank_statistics <- sample_level %>%
  group_by(Molecule) %>%
  summarise(
    Blank = mean_area[ESSENTIAL_sample_id == "SAME7905926"],
    Sample = mean(mean_area[sample_type == "Sample"], na.rm = TRUE),
    Ratio_blank_sample = Blank / Sample,
    .groups = "drop"
  )


###### 5.3) Blank contribution sample by sample ######

#Filtering
blank_filter_data <- filtered_data %>%
  filter(
    sample_type == "Sample" |
      ESSENTIAL_sample_id == "SAME7905926"
  ) %>%
  mutate(
    Area_per_L = as.numeric(Area_per_L)
  )


#Averaging technical triplicates
sample_means <- blank_filter_data %>%
  group_by(
    Molecule,
    ESSENTIAL_sample_id
  ) %>%
  summarise(
    mean_area = mean(Area_per_L, na.rm = TRUE),
    n = sum(!is.na(Area_per_L)),
    .groups = "drop"
  )

sample_means
blank_means <- sample_means %>%
  filter(ESSENTIAL_sample_id == "SAME7905926") %>%
  select(
    Molecule,
    blank_mean = mean_area
  )

blank_means


sample_blank_ratios <- sample_means %>%

filter(ESSENTIAL_sample_id != "SAME7905926") %>%

left_join(
    blank_means,
    by = "Molecule"
  ) %>%
  

#Blank/Sample ratio
mutate(
    Ratio_blank_sample = blank_mean / mean_area
  )

sample_blank_ratios

sample_blank_ratios %>%
  filter(Molecule == "DMSP") %>%
  arrange(desc(Ratio_blank_sample))


plot_ratio <- sample_blank_ratios %>%
  filter(Molecule == "MSA")


ggplot(
  plot_ratio,
  aes(
    x = reorder(ESSENTIAL_sample_id, Ratio_blank_sample),
    y = Ratio_blank_sample
  )
) +
  
  geom_point(
    size = 2.5
  ) +
  
  geom_hline(
    yintercept = 0.2,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dotted"
  ) +
  
  labs(
    x = "Sample ID",
    y = "Blank / sample ratio",
    title = "MSA"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


###### 5.4) Average technical replicates and blank substraction #####

#Samples: mean of technical triplicates
sample_means <- filtered_data %>%
  filter(sample_type == "Sample") %>%
  group_by(
    ESSENTIAL_sample_id,
    Molecule
  ) %>%
  summarise(
    mean_area_sample = mean(Area_per_L, na.rm = TRUE),
    n_replicates = sum(!is.na(Area_per_L)),
    .groups = "drop"
  )


#Procedural blank: mean of its 3 technical injections
blank_means <- filtered_data %>%
  filter(ESSENTIAL_sample_id == "SAME7905926") %>%
  group_by(Molecule) %>%
  summarise(
    mean_area_blank = mean(Area_per_L, na.rm = TRUE),
    n_blank_replicates = sum(!is.na(Area_per_L)),
    .groups = "drop"
  )

blank_corrected <- sample_means %>%
  left_join(
    blank_means,
    by = "Molecule"
  ) %>%
  mutate(
    Area_blank_corrected_raw =
      mean_area_sample - mean_area_blank,
#Negative values to 0
    Area_blank_corrected =
      pmax(Area_blank_corrected_raw, 0)
  )

#CSV
write_csv(blank_corrected, "2_Table_peak_area_MSA_after_blank_correction.csv")

##### 6) Average onboard technical replicates #####


metadata_samples <- read.table(
  "Metadata_final_with_operators.txt",
  sep = "\t"
) %>%
  row_to_names(row_number = 1)


metadata_onboard <- metadata_samples %>%
  select(
    ESSENTIAL_sample_id,
    ESSENTIAL_station_label,
    log_samples_with_barcode,
    ESSENTIAL_event_datetime_midpoint,
    ESSENTIAL_event_latitude_N_midpoint,
    ESSENTIAL_event_longitude_E_midpoint,
    ESSENTIAL_elevation_nominal_median_m,
    event_label
  ) %>%
  distinct()

blank_corrected_metadata <- blank_corrected %>%
  left_join(
    metadata_onboard,
    by = "ESSENTIAL_sample_id"
  )


onboard_groups <- blank_corrected_metadata %>%
  group_by(
    Molecule,
    ESSENTIAL_station_label,
    log_samples_with_barcode,
    ESSENTIAL_event_datetime_midpoint,
    ESSENTIAL_event_latitude_N_midpoint,
    ESSENTIAL_event_longitude_E_midpoint,
    ESSENTIAL_elevation_nominal_median_m,
    event_label
  ) %>%
  summarise(
    mean_area_onboard = mean(
      Area_blank_corrected,
      na.rm = TRUE
    ),
    
    n_onboard_replicates = n_distinct(ESSENTIAL_sample_id),
    
    .groups = "drop"
  )

print(onboard_groups)

#Export of the data:
write_csv(onboard_groups, "3_Table_peak_area_MSA_after_average_technical_replicates.csv")

check_onboard_replicates <- blank_corrected_metadata %>%
  group_by(
    ESSENTIAL_station_label,
    log_samples_with_barcode,
    ESSENTIAL_event_datetime_midpoint,
    ESSENTIAL_event_latitude_N_midpoint,
    ESSENTIAL_event_longitude_E_midpoint,
    ESSENTIAL_elevation_nominal_median_m,
    event_label
  ) %>%
  summarise(
    n_samples = n_distinct(ESSENTIAL_sample_id),
    sample_ids = paste(
      unique(ESSENTIAL_sample_id),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

View(check_onboard_replicates)

