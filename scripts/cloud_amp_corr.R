# All require Libraries for this analysis
library(dplyr)
library(readxl)
library(tibble)




# This contains image anlaysis using high-content imager
cloud <- read_csv("data/combined_data.csv", show_col_types = FALSE)
# This file contains the running wheel amplitude data from Jenny
amp <- read_excel("data/Running wheel Cohort 1-21 combined cosinor data.xlsx", sheet = 1)



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
                     "Meg3_Xist_Area", "Meg3_Xist_Int.Intensity")

animal_avg <- cloud %>%
  group_by(Animal_ID, Sex, Genotype, Light_Expo) %>%
  summarise(across(all_of(cols_to_average), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

animal_avg <- animal_avg %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "PWS", "COMP"))) %>%  # adjust labels to match yours
  arrange(Genotype)

View(animal_avg)
# Write a CSV file containing the calculated average of all cells per animal.
write.csv(animal_avg, "animal_averages.csv", row.names = FALSE)

cell_counts <- cloud %>%
  count(Animal_ID, Site.ID)




# 
library(ggplot2)
ggplot(animal_avg, aes(x = DAPI_Area, fill = Sex)) +
  geom_histogram(bins = 30, color = "black") +
  facet_grid(Sex ~ Genotype) +
  labs(title = "Distribution of DAPI_Area by Genotype and Sex",
       x = "DAPI_Area", y = "Count") +
  theme_minimal()


