library(readxl)
library(writexl)
library(tidyverse)


data <- read_excel("RxnsPaths_input.xlsx")

# Create a list to store the new results
new_data_list <- list()

# Identify the columns that contain reaction data
reaction_columns <- grep("^R_", colnames(data), value = TRUE)

# Process each reaction column
for (reaction_column in reaction_columns) {
  
  # Clean and split reaction entries in the column
  reactions <- gsub("^\\s+|\"|c\\(|\\)$", "", unlist(strsplit(as.character(data[[reaction_column]]), ',')))
  
  # Determine the last non-empty reaction entry
  last_filled_row <- max(which(!is.na(reactions) & reactions != "NA"))
  
  # Process each reaction up to the last non-empty entry
  for (j in 1:last_filled_row) {
    reaction <- reactions[j]
    
    if (!is.na(reaction) && reaction != "NA") {
      # Filter the data to find corresponding metabolisms for the reaction
      metabolisms <- data %>% filter(Rxns == reaction) %>% select(Metabolism_of_Model) %>% pull()
      
      if (length(metabolisms) == 0) {
        # If no corresponding metabolism is found, mark it as "NA"
        new_data_list[[length(new_data_list) + 1]] <- data.frame(
          Reaction = reaction,
          Metabolism_of_Model = "NA",
          Reaction_column = reaction_column
        )
      } else {
        # If corresponding metabolisms are found, add each to the new data list
        for (met in metabolisms) {
          new_data_list[[length(new_data_list) + 1]] <- data.frame(
            Reaction = reaction,
            Metabolism_of_Model = met,
            Reaction_column = reaction_column
          )
        }
      }
    } else {
      # If the reaction entry is NA, mark the metabolism as "NA"
      new_data_list[[length(new_data_list) + 1]] <- data.frame(
        Reaction = reaction,
        Metabolism_of_Model = "NA",
        Reaction_column = reaction_column
      )
    }
  }
}

# Combine all the new data into a single data frame
new_data <- do.call(rbind, new_data_list)

# Write the new data to an Excel file
write_xlsx(new_data, path = "RxnsPaths_output.xlsx")
