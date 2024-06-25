# Load required libraries
library(ggplot2)
library(reshape2)
library(readxl)
library(dplyr)

# Read the data from the Excel file
excel_file <- "Supplementary tables.xlsx"
data <- read_excel(excel_file, sheet = "Table S5")

# Rename the first column to "Reporter Metabolite Name"
colnames(data)[1] <- "Reporter Metabolite Name"

# Select the columns from the 2nd column onward and the first 20 rows
data <- data[1:20, 2:ncol(data)]

# Convert the data to long format
data_long <- melt(data, id.vars = "Reporter Metabolite Name", variable.name = "Sample", value.name = "Count")

# Add a new column to identify the group
data_long <- data_long %>%
  mutate(Group = ifelse(Sample %in% colnames(data)[2:6], "Strains", "Antibiotics"))

# Define color palettes
strains_colors <- c("1207552" = "darksalmon", 
                    "1428368" = "darkolivegreen",
                    "1457504" = "firebrick",
                    "34654" = "slateblue1",
                    "478810" = "#7ac5cd")

antibiotics_colors <- c("Amikacin sulfate" = "#008080", 
                        "Ciprofloxacin" = "#DAA520",
                        "Colistin" = "#FF6347",  # Tomato color
                        "Meropenem" = "#800080")  # Purple color

# Create the bubble plot
bubble <- ggplot(data_long, aes(x = Sample, y = `Reporter Metabolite Name`, size = Count, color = Sample)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(0.5, 50), name = "Count") +
  scale_color_manual(values = c(strains_colors, antibiotics_colors), guide = "none") +
  theme_minimal() +
  labs(x = "Strains and Antibiotics",
       y = "Reporter Metabolite Name",
       size = "Count") +
  theme(legend.position = "right",
        axis.title.x = element_text(margin = margin(t = 30), face = "bold", size = 30),
        axis.title.y = element_text(margin = margin(r = 20), face = "bold", size = 30),
        axis.text.x = element_text(size = 30, face = "bold", margin = margin(t = 20, b = 20)),
        axis.text.y = element_text(size = 30, face = "bold"),
        plot.title = element_text(face = "bold", size = 30),
        legend.title = element_text(face = "bold", size = 40),
        legend.text = element_text(face = "bold", size = 40),
        legend.margin = margin(t = 30),
        plot.margin = margin(t = 30, r = 30, b = 30, l = 30)) +
  theme(axis.title.x = element_text(size = 45),
        axis.title.y = element_text(size = 45))

# Save the bubble plot as a PNG file
ggsave("bubble_plot_RepMet.png", plot = bubble, dpi = 300, width = 40, height = 30, units = "in", limitsize = FALSE)
