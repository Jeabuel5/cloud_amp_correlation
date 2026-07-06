# All require Libraries for this analysis
library(dplyr)
library(readxl)
library(tibble)



# This contains image anlaysis using high-content imager
cloud <- read_csv("data/combined_data.csv", show_col_types = FALSE)
# This file contains the running wheel amplitude data from Jenny
amp <- read_excel("data/Running wheel Cohort 1-21 combined cosinor data.xlsx", sheet = 1)


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




