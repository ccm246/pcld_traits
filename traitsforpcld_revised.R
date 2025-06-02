
library(dplyr)
library(tidyverse)

#### Initial Set-up ####

# m <- read.csv('martin2019.csv', header=T)
# p <- read.csv('pcld_dump.csv', header=T)
# 
# s <- left_join(p,m,'Genus')
# head(s)
# unique(s$Taxonomic_group)
# 
# 
# traits = unique(trimws(unlist(strsplit(docs$Tags, split=";"))))

# run once to update input
setwd('~/Downloads/')
pcld <- read.csv('pcld_dump.csv', header=T)
names(pcld) <- tolower(names(pcld))
pcld$genus[pcld$genus == 'unknown'] <- NA
pcld <- pcld %>% drop_na(genus)
genuses <- unique(pcld[,71])

setwd('~/Desktop/TEEMS_2025/LivingDatabase/Trait Work/traits/unfilteredbypcld/')
a <- read.csv('arthropodtraits.csv', header=T)
names(a) <- tolower(names(a))
a$genus <- tolower(a$genus)
anew <- a[a$genus %in% genuses,]

p <- read.csv('emerypitfalltrap.csv', header=T)
names(p) <- tolower(names(p))
p$genus <- tolower(p$genus)
pnew <- p[p$genus %in% genuses,]

s <- read.csv('emerystickytrap.csv', header=T)
names(s) <- tolower(names(s))
s$genus <- tolower(s$genus)
snew <- s[s$genus %in% genuses,]

g <- read.csv('globalants.csv', header=T)
names(g) <- tolower(names(g))
g$genus <- tolower(g$genus)
gnew <- g[g$genus %in% genuses,]

j <- read.csv('jeppsson.csv', header=T)
names(j) <- tolower(names(j))
j$genus <- tolower(j$genus)
jnew <- j[j$genus %in% genuses,]

m <- read.csv('martin2019.csv', header=T)
names(m) <- tolower(names(m))
m$genus <- tolower(m$genus)
mnew <- m[m$genus %in% genuses,]

c <- read.csv('sesync.csv', header=T)
names(c) <- tolower(names(c))
c$genus <- tolower(c$genus)
cnew <- c[c$genus %in% genuses,]

t <- read.csv('tamburini.csv', header=T)
names(t) <- tolower(names(t))
t$genus <- tolower(t$genus)
tnew <- t[t$genus %in% genuses,]

ca <- read.csv('carabids_database.csv', header=T)
names(ca) <- tolower(names(ca))
ca$genus <- tolower(ca$genus)
ca$species <- tolower(ca$species)
canew <- ca[ca$genus %in% genuses,]


# write.csv(canew, 'newcarabids.csv')
# write.csv(anew, 'newarthropodtraits.csv')
# write.csv(snew, 'newemerystickytrap.csv')
# write.csv(pnew, 'newemerypitfalltrap.csv')
# write.csv(gnew, 'newglobalants.csv')
# write.csv(jnew, 'newjeppsson.csv')
# write.csv(mnew, 'newmartin2019.csv')
# write.csv(cnew, 'newsesync.csv')
# write.csv(tnew, 'newtamburini.csv')


# #### Try One ####
# 
# setwd('~/Desktop/TEEMS_2025/LivingDatabase/Trait Work/traits/')
# 
# 
# # Load necessary libraries
# library(dplyr)
# library(readxl)
# 
# # Define the path to your data folder
# data_path <- "~/Desktop/TEEMS_2025/LivingDatabase/Trait Work/traits"
# 
# # Step 1: Load all spreadsheets
# files <- list.files(data_path, pattern = ".csv", full.names = TRUE)
# data_list <- lapply(files, read.csv)
# 
# 
# # Step 2: Standardize column names across all dataframes
# data_list <- lapply(data_list, function(df) {
#   names(df) <- tolower(gsub(" ", "_", names(df))) # Lowercase and replace spaces with underscores
#   df
# })
# 
# # make sure all spreadsheets distinct
# data_list <- lapply(data_list, function(df) {
#   df <- distinct(df) 
#   df
# })
# 
# # Step 3: Harmonize taxonomic levels
# data_list <- lapply(data_list, function(df) {
#   df <- df %>%
#     mutate(
#       genus = ifelse(is.na(genus) & !is.na(species), sub(" .*", "", species), genus), # Extract genus from species
#       species = ifelse(is.na(species), paste0(genus, "_sp"), species) # Fill species as genus_sp when missing
#     )
#   df
# })
# 
# # Convert all columns in each data frame to categorical
# data_list <- lapply(data_list, function(df) {
#   df[] <- lapply(df, as.character)  # Convert each column to character
#   return(df)
# })
# str(data_list)
# 
# # # add source info to all trait variable
# # data_list <- lapply(data_list, function(df) {
# #   df %>% 
# #   rename_with(~ paste0(files[], sep='_',.), -c(genus, species))
# # })
# 
# # Step 4: Combine data into a single dataframe
# combined_data <- as_tibble(bind_rows(data_list, .id = "source"))
# 
# # Step 5: Wide-format trait unification
# # Identify all unique traits across datasets
# all_traits <- setdiff(names(combined_data), c("source", "genus", "species"))
# 
# # Reshape to ensure all traits are present in all rows
# combined_data <- combined_data %>%
#   dplyr::select(genus, species, all_of(all_traits)) %>%
#   pivot_longer(cols = all_of(all_traits), names_to = "trait", values_to = "value") %>%
#   group_by(genus, species, trait) %>%
#   #summarise(
#   #  value = ifelse(first(na.omit(value))),
#   #  .groups = "drop"
#   # ) %>%
#   pivot_wider(names_from = trait, values_from = value)
# 
# # Step 6: Resolve genus-level duplicates
# resolved_data <- combined_data %>%
#   group_by(genus) %>%
#   mutate(
#     across(
#       starts_with("trait_"), 
#       ~ ifelse(is.na(.), mean(., na.rm = TRUE), .), 
#       .names = "genus_{col}"
#     )
#   ) %>%
#   ungroup()
# 
# # Step 7: Export the final combined dataset
# df <- apply(resolved_data, 2, as.character)
# df <- as.data.frame(df)
# 
# df[] <- lapply(df, function(x) {
#   if (is.character(x)) gsub(",", "", x) else x
# })
# df[] <- lapply(df, function(x) {
#   if (is.character(x)) gsub("NA ", "", x) else x
# })
# df[] <- lapply(df, function(x) {
#   if (is.character(x)) gsub(")", "", x) else x
# })
# df[] <- lapply(df, function(x) {
#   if (is.character(x)) gsub("c(", "", x) else x
# })
# 
#  #for some reason this messes up the dataframe and so I am not automating clean up of the df
#  # df <- lapply(df, function(x) {
#  #   gsub("NA," , "", x)
#  # })
#  # df <- lapply(df, function(x) {
#  #   gsub("," , "", x)
#  # })
#  # df <- apply(resolved_data, 2, as.character)
# 
# # write.csv(df, file="Combined_Insect_Traits_3.csv")
# 
# df2 <- read.csv('Combined_Insect_Traits_2_new.csv', header=T)
# f <- df2 %>% select(where(~ any(!is.na(.))))
# 
# # Colleen then does some postprocessing in excel to clean up lists and some characters
#   # removes "NA," phrases, so that if there is a value that is reported as NA many times, it only says NA once
#     # n.b., this does not remove lists that read as "1, 5, NA", so if values are varied and NA is included it is still listed. Could be removed though
#   # removes parentheses from values as well so that the dataset is more readable moving forward
#   # removes quotation marks from values as well so that the dataset is more readable moving forward
#   # removes all commas
# 

#### Try Two -- Start here ####

setwd('~/Desktop/TEEMS_2025/LivingDatabase/Trait Work/traits/pcld_traits/')

library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
library(purrr)

# Define the path to your data folder
data_path <- "~/Desktop/TEEMS_2025/LivingDatabase/Trait Work/traits/pcld_traits/"

# Step 1: Load all spreadsheets
files <- list.files(data_path, pattern = ".csv", full.names = TRUE)
names(files) <- tools::file_path_sans_ext(basename(files)) # Name by file base

data_list <- lapply(files, read.csv)

# Step 2: Standardize column names across all dataframes
data_list <- lapply(data_list, function(df) {
  names(df) <- tolower(gsub(" ", "_", names(df)))
  df
})

# Step 3: Ensure rows are unique
data_list <- lapply(data_list, distinct)

# Step 4: Harmonize taxonomic levels
data_list <- lapply(data_list, function(df) {
  df %>%
    mutate(
      genus = ifelse(is.na(genus) & !is.na(species), sub(" .*", "", species), genus),
      species = ifelse(is.na(species), paste0(genus, "_sp"), species)
    )
})

# Step 4.5: Make sure all data frames have genus/species as character
data_list <- lapply(data_list, function(df) {
  df$genus <- as.character(df$genus)
  df$species <- as.character(df$species)
  df
})

# Step 4.6: Drop rows with no usable genus
data_list <- lapply(data_list, function(df) {
  df <- df %>%
    mutate(genus = str_trim(as.character(genus))) %>%
    filter(!is.na(genus) & genus != "")
})

# make sure they are distinct
data_list <- lapply(data_list, function(df) {
  df %>% distinct(genus, species, .keep_all = TRUE)
})

# Optional: Tag traits with filename to avoid collisions
data_list <- Map(function(df, fname) {
  trait_cols <- setdiff(names(df), c("genus", "species"))
  df <- df %>%
    rename_with(~ paste0(fname, "_", .), all_of(trait_cols))
}, data_list, names(files))

# Step 5: Combine into one long-format table
combined_data <- bind_rows(data_list, .id = "source")

# Step 6: Make all columns character (we'll clean later)
combined_data[] <- lapply(combined_data, as.character)

# Step 6.5: Inspect trait names
trait_cols <- setdiff(names(combined_data), c("source", "genus", "species"))
trait_cols[duplicated(trait_cols)]

# Step 7: Pivot longer
trait_cols <- setdiff(names(combined_data), c("source", "genus", "species"))
long_data <- combined_data %>%
  pivot_longer(cols = all_of(trait_cols), names_to = "trait", values_to = "value") %>%
  filter(!is.na(value) & value != "") %>%
  group_by(genus, species, trait) %>%
  summarise(value = paste(unique(value), collapse = "; "), .groups = "drop")  # Combine multi-values with semicolon

# Step 8: Pivot back to wide
wide_data <- long_data %>%
  pivot_wider(names_from = trait, values_from = value)

# Step 8.5: Check for duplicate taxa
dup_taxa <- wide_data %>%
  group_by(genus, species) %>%
  tally() %>%
  filter(n > 1)

print(dup_taxa)

# Step 9: Attempt to clean numeric-like columns
wide_data[] <- lapply(wide_data, function(x) {
  if (all(str_detect(x, "^\\d+(\\.\\d+)?(; \\d+(\\.\\d+)?)*$") | is.na(x))) {
    as.character(x)  # leave as is for now
  } else {
    str_replace_all(x, "[\\(\\),]", "") %>%
      str_replace_all("^NA ", "") %>%
      str_trim()
  }
})

# Step 9.5: Clean up characters again
wide_data[] <- lapply(wide_data, function(x) {
  if (is.character(x)) {
    x <- str_remove_all(x, "c\\(")         # remove 'c('
    x <- str_remove_all(x, "\\)")          # remove ')'
    x <- str_replace_all(x, ",", ";")      # convert comma to semicolon
    x <- str_replace_all(x, "NA; ?", "")   # clean trailing NA;
    x <- str_trim(x)
  }
  return(x)
})

# Step 9.5: Remove columns that are entirely NA or empty
wide_data <- wide_data[, colSums(!is.na(wide_data) & wide_data != "") > 0]

# Step 9.6: Flag genus/species with no trait data
wide_data$no_trait_data <- apply(wide_data[, -c(1, 2)], 1, function(row) {
  all(is.na(row) | row == "")
})

# Step 9.7: Sanity Checks
# Check that each column now has only one kind of trait
head(names(wide_data), 20)

# Look for any cells with multiple values again
any(grepl(";", wide_data))

# How many taxa have no trait data?
sum(wide_data$no_trait_data)

# Step 9.8: Split wide_data into trait type CSVs by column suffix
# Identify the trait columns (excluding genus, species, no_trait_data)
trait_cols <- setdiff(names(wide_data), c("genus", "species", "no_trait_data"))

# Step 10: Option 1 - export files that are organized by trait type
setwd('~/Desktop/TEEMS_2025/LivingDatabase/Trait Work/traits/pcld_traits/CleaningOutput/')

# Extract unique suffixes (like "_source", "_habitat", etc.)
suffixes <- trait_cols %>%
  str_extract("_[^_]+$") %>%
  na.omit() %>%
  unique()

# Loop through each suffix to create and export a subset dataframe
for (suffix in suffixes) {
  subset_df <- wide_data %>%
    select(genus, species, ends_with(suffix))
  
  # Drop if it ends up only with genus/species (i.e., no traits matched)
  if (ncol(subset_df) > 2) {
    file_name <- paste0("TEEMS_traits", suffix, ".csv")
    write.csv(subset_df, file_name, row.names = FALSE)
  }
}

# # Step 10: Export or continue using
# write.csv(wide_data, "TEEMS_combined_traits_cleaned.csv", row.names = FALSE)

