# Load required libraries
library(readxl)
library(dplyr)
library(writexl)


# Constants
total_rxn <- 1659  # Number of reactions with non-zero min and max values and assigned genes after FVA analysis

# Read data from Excel file
excel_file <- "fishersTest.xlsx"
path_rxn_counts <- read_excel(excel_file, sheet = "PathRxnCounts")
path_rxn_counts <- setNames(path_rxn_counts$Number_of_Reactions, path_rxn_counts$KEGG_Pathways)

cond_changed_rxn <- read_excel(excel_file, sheet = "PerturbedRxnCounts")
cond_changed_rxn <- setNames(cond_changed_rxn$Number_of_Perturbed_Reactions, cond_changed_rxn$Condition)

cond_path_rxn_counts <- read_excel(excel_file, sheet = "CondPathRxnCounts")

# Function to create dictionaries for each condition
create_dictionaries <- function(data) {
  dictionaries <- list()
  for (i in seq(1, ncol(data), by = 2)) {
    key_col <- colnames(data)[i]
    val_col <- colnames(data)[i + 1]
    dict_name <- gsub("\\..*", "", key_col)
    
    non_na_rows <- complete.cases(data[, c(key_col, val_col)])
    keys <- data[non_na_rows, key_col] %>% pull() %>% as.character()
    values <- data[non_na_rows, val_col] %>% pull() %>% as.numeric()
    
    dictionaries[[dict_name]] <- setNames(values, keys)
  }
  dictionaries
}

dictionaries <- create_dictionaries(cond_path_rxn_counts)

# Function to perform Fisher's exact test
calculate_p_values <- function(dict, cond_changed_rxn, path_rxn_counts, total_rxn) {
  p_values <- list()
  for (condition in names(dict)) {
    if (condition %in% names(cond_changed_rxn)) {
      for (path in names(path_rxn_counts)) {
        if (path %in% names(dict[[condition]]) && path %in% names(path_rxn_counts)) {
          contingency_table <- matrix(c(
            dict[[condition]][[path]],
            cond_changed_rxn[[condition]] - dict[[condition]][[path]],
            path_rxn_counts[[path]] - dict[[condition]][[path]],
            total_rxn - path_rxn_counts[[path]] - (cond_changed_rxn[[condition]] - dict[[condition]][[path]])
          ), nrow = 2)
          
          p_value <- fisher.test(contingency_table, alternative = "greater")$p.value
          result_identifier <- paste(condition, path, sep = "_")
          p_values[[result_identifier]] <- p_value
        }
      }
    }
  }
  p_values
}

p_values_list <- calculate_p_values(dictionaries, cond_changed_rxn, path_rxn_counts, total_rxn)

# Filter significant conditions (p-value < 0.05)
significant_conditions <- lapply(p_values_list, function(p_value) {
  if (p_value < 0.05) {
    return(list(condition = names(p_values_list)[which(p_values_list == p_value)], p_value = p_value))
  } else {
    return(NULL)
  }
}) %>% Filter(Negate(is.null), .)

# Convert significant conditions to a data frame and write to Excel
df_significant_conditions <- do.call(rbind, lapply(significant_conditions, as.data.frame))
write_xlsx(df_significant_conditions, "significant_conditions_pathways.xlsx")

# Convert p-values list to a data frame and write to Excel
df_pval <- do.call(rbind, lapply(names(p_values_list), function(condition) {
  data.frame(condition = condition, p_value = p_values_list[[condition]])
}))
write_xlsx(df_pval, "all_conditions_pathways.xlsx")
