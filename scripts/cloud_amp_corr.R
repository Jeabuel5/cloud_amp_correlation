# All require Libraries for this analysis
library(dplyr)
library(readxl)
library(tibble)



# This contains image anlaysis using high-content imager
cloud <- read_csv("data/combined_data.csv", show_col_types = FALSE)

# This file contains the running wheel amplitude data from Jenny
amp <- read_csv("data/Running_wheel_cosinor_data_cleaned.csv", show_col_types = FALSE)
amp <- amp %>% rename(Note = `...15`)
amp <- amp %>% select(-Note)
write_csv(amp, "data/Running_wheel_cosinor_data_cleaned_v2.csv")


# Normalize the lncRNA_Pairs
cloud$Meg3_Snhg14_Area_norm <- (cloud$Meg3_Snhg14_Area / cloud$DAPI_Area) * 100
cloud$Meg3_Xist_Area_norm   <- (cloud$Meg3_Xist_Area / cloud$DAPI_Area) * 100
cloud$Snhg14_Xist_Area_norm <- (cloud$Snhg14_Area / cloud$DAPI_Area) * 100



# Fix Genotypes column. 
cloud <- cloud %>%
  #rename(Genotype = Genotypes) %>%
  mutate(Genotype = sub("\\s\\d+$", "", Genotype),
         Genotype = recode(Genotype, "DEL" = "PWS"))

View(cloud)


# Rename "Sex" column in cloud~ instead of Male1, we replace the actual animal number for each samples
id_map_11 <- tibble(
  Sex = c("Male1", "Female1", "Female2", "Female3",
          "Male1", "Male2", "Male3", "Female1", "Female2", "Female3",
          "Male1", "Male2", "Male3", "Female1", "Female2", "Female3"),
  Genotype = c("WT", "WT", "WT", "WT",
               "PWS", "PWS", "PWS", "PWS", "PWS", "PWS",
               "COMP", "COMP", "COMP", "COMP", "COMP", "COMP"),
  Light_Expo = c(rep(11, 16)),
  Animal_ID = c(499, 831, 777, 858,
                609, 627, 636, 836, 837, 854,
                656, 607, 629, 839, 870, 850)
)

id_map_12 <- tibble(
  Sex = c("Female1", "Female2", "Female3", "Male1", "Male2", "Male3",
          "Female1", "Female2", "Female3", "Male1", "Male2", "Male3",
          "Female1", "Female2", "Female3", "Male1", "Male2", "Male3"),
  Genotype = c("WT", "WT", "WT", "WT", "WT", "WT",
               "PWS", "PWS", "PWS", "PWS", "PWS", "PWS",
               "COMP", "COMP", "COMP", "COMP", "COMP", "COMP"),
  Light_Expo = c(rep(12, 18)),
  Animal_ID = c(836, 740, 778, 730, 509, 506,
                753, 881, 825, 579, 503, 647,
                810, 738, 797, 610, 593, 634)
)

id_map <- bind_rows(id_map_11, id_map_12)

cloud <- cloud %>%
  left_join(id_map, by = c("Sex", "Genotype", "Light_Expo"))



# Calculate the average value per Animal_ID
cols_to_average <- c("DAPI_Area", "DAPI_Int.Intensity", "Meg3_Area", "Meg3_Int.Intensity",
                     "Snhg14_Area", "Snhg14_Int.Intensity", "Xist_Area", "Xist_Int.Intensity",
                     "Tri.coloc_Area", "Tri.coloc_Int.Intensity", "Snhg14_Xist_Area",
                     "Snhg14_Xist_Int.Intensity", "Meg3_Snhg14_Area", "Meg3_Snhg14_Int.Intensity",
                     "Meg3_Xist_Area", "Meg3_Xist_Int.Intensity",
                     "Meg3_norm", "Snhg14_norm", "Xist_norm", "Tri.coloc_norm", "Meg3_Xist_Area_norm",
                     "Meg3_Snhg14_Area_norm", "Snhg14_Xist_Area_norm")  # added normalized cols

animal_avg <- cloud %>%
  group_by(Animal_ID, Sex, Genotype, Light_Expo) %>%
  summarise(across(all_of(cols_to_average), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

animal_avg <- animal_avg %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "PWS", "COMP"))) %>%
  arrange(Genotype)

# Rename the columns and some labels
animal_avg <- animal_avg %>% rename(Cohort = ID)
animal_avg <- animal_avg %>%
  mutate(Genotype = recode(Genotype,
                           "WT" = "WT/WT",
                           "PWS" = "HET/WT",
                           "COMP" = "HET/TG"
  ))

library(stringr)
library(tidyverse)

animal_avg <- animal_avg %>%
  mutate(Sex = case_when(
    str_detect(Sex, "^Male") ~ "M",
    str_detect(Sex, "^Female") ~ "F",
    TRUE ~ Sex
  )) %>%
  rename(ID = Animal_ID,
         Cycle = Light_Expo)

View(animal_avg)
write.csv(animal_avg, "animal_averages.csv", row.names = FALSE)

cell_counts <- cloud %>%
  count(Animal_ID, Site.ID)


animal_avg <- animal_avg %>%
  mutate(Cycle = factor(Cycle, levels = c("12", "11")))

# Plot using the average calculated Meg3_norm (normalized to area)
# Meg3_RNA cloud size plot
ggplot(animal_avg, aes(x = Genotype, y = Meg3_norm, fill = Sex, color = Sex)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              shape = 21,
              color = "black",
              size = 3.0,
              stroke = 0.8,
              alpha = 1.0) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3, color = "black") +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  scale_fill_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  scale_color_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  facet_wrap(~ Cycle, labeller = labeller(Cycle = c("11" = "11:11", "12" = "12:12"))) + # updated from Sample -> Cycle (confirm this is right)
  labs(x = "Genotype", y = "Meg3 RNA-cloud size (% nuclear area)", fill = "Sex", color = "Sex") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )



# Snhg14_RNA cloud size plot
ggplot(animal_avg, aes(x = Genotype, y = Snhg14_norm, fill = Sex, color = Sex)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              shape = 21,
              color = "black",
              size = 3.0,
              stroke = 0.8,
              alpha = 1.0) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3, color = "black") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 0.5)) +
  scale_fill_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  scale_color_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  facet_wrap(~ Cycle, labeller = labeller(Cycle = c("11" = "11:11", "12" = "12:12"))) + # updated from Sample -> Cycle (confirm this is right)
  labs(x = "Genotype", y = "Snhg14 RNA-cloud size (% nuclear area)", fill = "Sex", color = "Sex") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )



# Xist_RNA cloud size plot
animal_avg %>%
  filter(Sex == "F") %>%
  ggplot(aes(x = Genotype, y = Xist_norm)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, color = "black", fill = "pink") +
  geom_jitter(width = 0.0, size = 3.0, alpha = 1.0, color = "black") +
  scale_y_continuous(breaks = seq(0, 15, by = 5)) +
  coord_cartesian(ylim = c(NA, 15)) +
  facet_wrap(~ Cycle, labeller = labeller(Cycle = c("11" = "11:11", "12" = "12:12"))) +
  labs(x = "Genotype", y = "Xist RNA-cloud size (% nuclear area)") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )



# ----------------------



# Meg3_Xist_Area_norm COLOCALIZATION
ggplot(animal_avg, aes(x = Genotype, y = Meg3_Xist_Area_norm, fill = Sex, color = Sex)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              shape = 21,
              color = "black",
              size = 3.0,
              stroke = 0.8,
              alpha = 1.0) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3, color = "black") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
  scale_fill_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  scale_color_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  facet_wrap(~ Cycle, labeller = labeller(Cycle = c("11" = "11:11", "12" = "12:12"))) + # updated from Sample -> Cycle (confirm this is right)
  labs(x = "Genotype", y = "Meg3_Xist colocalization (% nuclear area)", fill = "Sex", color = "Sex") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )



# Meg3_Snhg14_Area_norm COLOCALIZATION
ggplot(animal_avg, aes(x = Genotype, y = Meg3_Snhg14_Area_norm, fill = Sex, color = Sex)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              shape = 21,
              color = "black",
              size = 3.0,
              stroke = 0.8,
              alpha = 1.0) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3, color = "black") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  scale_fill_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  scale_color_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  facet_wrap(~ Cycle, labeller = labeller(Cycle = c("11" = "11:11", "12" = "12:12"))) + # updated from Sample -> Cycle (confirm this is right)
  labs(x = "Genotype", y = "Meg3--Snhg14 colocalization (% nuclear area)", fill = "Sex", color = "Sex") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )



# Snhg14_Xist_Area_norm COLOCALIZATION
ggplot(animal_avg, aes(x = Genotype, y = Snhg14_Xist_Area_norm, fill = Sex, color = Sex)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              shape = 21,
              color = "black",
              size = 3.0,
              stroke = 0.8,
              alpha = 1.0) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3, color = "black") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  scale_fill_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  scale_color_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  facet_wrap(~ Cycle, labeller = labeller(Cycle = c("11" = "11:11", "12" = "12:12"))) + # updated from Sample -> Cycle (confirm this is right)
  labs(x = "Genotype", y = "Snhg14_Xist colocalization (% nuclear area)", fill = "Sex", color = "Sex") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )



# Tri.coloc_norm COLOCALIZATION
ggplot(animal_avg, aes(x = Genotype, y = Tri.coloc_norm, fill = Sex, color = Sex)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              shape = 21,
              color = "black",
              size = 3.0,
              stroke = 0.8,
              alpha = 1.0) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3, color = "black") +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.5)) +
  scale_fill_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  scale_color_manual(values = c("F" = "pink", "M" = "#56B4E9")) +
  facet_wrap(~ Cycle, labeller = labeller(Cycle = c("11" = "11:11", "12" = "12:12"))) + # updated from Sample -> Cycle (confirm this is right)
  labs(x = "Genotype", y = "All probes colocalization (% nuclear area)", fill = "Sex", color = "Sex") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )




# ------------

# packages
library(dplyr)
library(purrr)
library(broom)
library(writexl)

# make sure these are factors
animal_avg$Sex <- factor(animal_avg$Sex)
animal_avg$Cycle <- factor(animal_avg$Cycle)
animal_avg$Genotype <- factor(animal_avg$Genotype)

# list of measurements
measures <- c("Meg3_norm", "Snhg14_norm", "Xist_norm",
              "Meg3_Snhg14_Area_norm", "Meg3_Xist_Area_norm", "Snhg14_Xist_Area_norm")

# function to run ANOVA and tidy results
run_anova <- function(response, formula_text, anova_type) {
  fit <- aov(as.formula(paste(response, formula_text)), data = animal_avg)
  tidy(fit) %>%
    mutate(
      measure = response,
      model = anova_type
    )
}

# 2-way ANOVA: Sex + Cycle
anova_2way_df <- map_dfr(measures, ~ run_anova(.x, "~ Sex * Genotype", "2-way"))

# 3-way ANOVA: Sex + Cycle + Genotype
# anova_3way_df <- map_dfr(measures, ~ run_anova(.x, "~ Sex * Cycle * Genotype", "3-way"))

# combine into one table
anova_all <- bind_rows(anova_2way_df)

# write to Excel file
write_xlsx(list(
  ANOVA_2way = anova_2way_df,
  #ANOVA_3way = anova_3way_df,
  #ANOVA_all = anova_all
), path = "anova_results.xlsx")





# ----------
# ONLY 2-WAY ANOVA on 12:12

# packages
library(dplyr)
library(purrr)
library(broom)
library(writexl)

# keep only cycle 12 and make sure factors are set
animal_avg_cycle12 <- animal_avg %>%
  filter(Cycle == 12) %>%
  mutate(
    Sex = factor(Sex),
    Genotype = factor(Genotype)
  )

# list of measurements
measures <- c("Meg3_norm", "Snhg14_norm", "Xist_norm",
              "Meg3_Snhg14_Area_norm", "Meg3_Xist_Area_norm", "Snhg14_Xist_Area_norm")

# function to run ANOVA and tidy results
run_anova <- function(response, formula_text, anova_type, data) {
  fit <- aov(as.formula(paste(response, formula_text)), data = data)
  tidy(fit) %>%
    mutate(
      measure = response,
      model = anova_type
    )
}

# 2-way ANOVA for cycle 12 only: Sex * Genotype
anova_2way_df <- map_dfr(
  measures,
  ~ run_anova(.x, "~ Sex * Genotype", "2-way_cycle12", animal_avg_cycle12)
)

# write to Excel file
write_xlsx(list(
  ANOVA_2way_cycle12 = anova_2way_df
), path = "anova_results_cycle12.xlsx")
















# --------------------------------------------------------------------------------
# Correlation Analysis of Meg3 RNA-cloud with RW amplitude
# Merge amp data with the animal_avg and add the Meg3_norm column to it.

# Rename Cohort back to ID in cloud_subset's source
animal_avg <- animal_avg %>%
  rename(ID = Cohort)

# Rebuild cloud_subset with the correct column
cloud_subset <- animal_avg %>%
  select(ID, Sex, Cycle, Genotype, Meg3_norm) %>%
  mutate(Cycle = as.numeric(as.character(Cycle)))

# Join on ID (+ Sex, Cycle, Genotype as extra safety checks)
merged_data <- amp %>%
  left_join(cloud_subset, by = c("ID", "Sex", "Cycle", "Genotype"))

sum(!is.na(merged_data$Meg3_norm))

write_csv(merged_data, "data/merged_amp_cloud_data.csv")



complete_data <- merged_data %>%
  filter(!is.na(Amplitude_Base) & !is.na(Meg3_norm))

nrow(complete_data)  # should be ~31

cor_test <- cor.test(complete_data$Amplitude_Base, complete_data$Meg3_norm)
cor_test


# This is to check if I need to normalize my data by log transform
hist(complete_data$Meg3_norm)
hist(complete_data$Amplitude_Base)
shapiro.test(complete_data$Meg3_norm)
shapiro.test(complete_data$Amplitude_Base)










library(ggplot2)
library(dplyr)

# ---- Set genotype order ----
geno_order <- c("WT/WT", "HET/WT", "HET/TG")

complete_data <- complete_data %>%
  mutate(Genotype = factor(Genotype, levels = geno_order))

# ---- Per-genotype correlation stats ----
cor_by_genotype <- complete_data %>%
  group_by(Genotype) %>%
  summarise(
    r = round(cor(Amplitude_Base, Meg3_norm, use = "complete.obs"), 2),
    p = signif(cor.test(Amplitude_Base, Meg3_norm)$p.value, 2),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", r, ", p = ", p),
    facet_label = paste0(Genotype, " (n=", n, ")")
  )

# ---- Attach facet labels back to complete_data, preserving genotype order ----
complete_data <- complete_data %>%
  left_join(cor_by_genotype %>% select(Genotype, facet_label), by = "Genotype") %>%
  mutate(facet_label = factor(
    facet_label,
    levels = cor_by_genotype$facet_label[order(match(cor_by_genotype$Genotype, geno_order))]
  ))

# ---- Plot ----
ggplot(complete_data, aes(x = Meg3_norm, y = Amplitude_Base)) +
  geom_point(size = 2.3, alpha = 0.75, color = "#2C5F8A") +
  geom_smooth(method = "lm", se = TRUE, color = "#D9534F", fill = "grey85", linewidth = 0.9) +
  facet_wrap(~ facet_label) +
  labs(
    title = "Meg3 RNA-cloud vs. Baseline Amplitude",
    subtitle = "By genotype",
    x = "Meg3_norm",
    y = "Amplitude_Base"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  geom_text(
    data = complete_data %>% distinct(facet_label, Genotype) %>%
      left_join(cor_by_genotype %>% select(Genotype, label), by = "Genotype"),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.5, size = 3.3, fontface = "bold", color = "black",
    inherit.aes = FALSE
  )






# Recalculate stats specifically for Amplitude_F4
cor_by_genotype_F4 <- complete_data %>%
  group_by(Genotype) %>%
  summarise(
    r = round(cor(Amplitude_F4, Meg3_norm, use = "complete.obs"), 2),
    p = signif(cor.test(Amplitude_F4, Meg3_norm)$p.value, 2),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", r, ", p = ", p))

ggplot(complete_data, aes(x = Meg3_norm, y = Amplitude_F4)) +
  geom_point(size = 2.3, alpha = 0.75, color = "#2C5F8A") +
  geom_smooth(method = "lm", se = TRUE, color = "#D9534F", fill = "grey85", linewidth = 0.9) +
  facet_wrap(~ facet_label) +
  labs(
    title = "Meg3 RNA_cloud vs. First 4 Amplitude",
    subtitle = "By genotype",
    x = "Meg3_norm",
    y = "Amplitude_F4"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  geom_text(
    data = complete_data %>% distinct(facet_label, Genotype) %>%
      left_join(cor_by_genotype_F4 %>% select(Genotype, label), by = "Genotype"),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.5, size = 3.3, fontface = "bold", color = "black",
    inherit.aes = FALSE
  )






# Recalculate correlation stats specifically for Amplitude_L4
cor_by_genotype_L4 <- complete_data %>%
  group_by(Genotype) %>%
  summarise(
    r = round(cor(Amplitude_L4, Meg3_norm, use = "complete.obs"), 2),
    p = signif(cor.test(Amplitude_L4, Meg3_norm)$p.value, 2),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", r, ", p = ", p))

ggplot(complete_data, aes(x = Meg3_norm, y = Amplitude_L4)) +
  geom_point(size = 2.3, alpha = 0.75, color = "#2C5F8A") +
  geom_smooth(method = "lm", se = TRUE, color = "#D9534F", fill = "grey85", linewidth = 0.9) +
  facet_wrap(~ facet_label) +
  labs(
    title = "Meg3 RNA-cloud vs. Last 4 Amplitude",
    subtitle = "By genotype",
    x = "Meg3_norm",
    y = "Amplitude_L4"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  geom_text(
    data = complete_data %>% distinct(facet_label, Genotype) %>%
      left_join(cor_by_genotype_L4 %>% select(Genotype, label), by = "Genotype"),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.5, size = 3.3, fontface = "bold", color = "black",
    inherit.aes = FALSE
  )








bind_rows(
  cor_by_genotype %>% mutate(Phase = "Baseline"),
  cor_by_genotype_F4 %>% mutate(Phase = "First 4"),
  cor_by_genotype_L4 %>% mutate(Phase = "Last 4")
) %>%
  select(Phase, Genotype, r, p, n) %>%
  arrange(Genotype, Phase)


