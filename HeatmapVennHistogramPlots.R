# Load necessary libraries
library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)  # Required for plot.margin
library(VennDiagram)
library(tidyr)
library(gridExtra)

# Read Excel file
data_heatmap <- read_excel("HeatmapVennHistogram.xlsx", sheet = "Heatmap")

# Reshape the data frame
dfm <- melt(data_heatmap, id.vars = "Pathway")

# Create a column with counts of values greater than 1 and sort
dfm <- dfm %>%
  group_by(Pathway) %>%
  mutate(above_two_count = sum(value > 1)) %>%
  ungroup() %>%
  arrange(Pathway)

# Plot the heatmap
heatmap_plot <- ggplot(dfm, aes(variable, Pathway, fill = value)) +
  geom_tile(color = "white", width = 1.2, height = 0.9) +
  scale_fill_gradient(low = "grey100", high = "lightblue3") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 37, face = "bold"),
    axis.text.y = element_text(size = 40, face = "bold"),
    text = element_text(size = 40, face = "bold"),
    axis.title.y = element_blank(),
    plot.margin = unit(c(1, 6, 1, 1), "cm"),
    legend.box.spacing = unit(1, "cm"),
    panel.grid.major.y = element_blank(),
    strip.text = element_blank(),
    legend.text = element_text(size = 25)
  ) +
  geom_text(aes(label = ifelse(value > 1, value, ""), x = variable, y = Pathway), color = "black", size = 10, hjust = 0.5, vjust = 0.5) +
  scale_y_discrete(limits = rev(unique(dfm$Pathway))) +
  facet_wrap(~ substr(variable, 1, 3), scales = "free_x", nrow = 1) +
  labs(fill = "Number of strains") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, label.theme = element_text(size = 35)))

print(heatmap_plot)

# Save the heatmap image
ggsave("Heatmap_significantPathways.png", plot = heatmap_plot, dpi = 300, width = 40, height = 30, units = "in", limitsize = FALSE)

# Read Venn data from Excel
venn_data <- read_excel("HeatmapVennHistogram.xlsx", sheet = "Venn")

# Define column names for the Venn diagram
column_names <- c("Amikacin sulfate", "Ciprofloxacin", "Colistin", "Meropenem")

# Extract data from selected columns
column_data <- lapply(column_names, function(name) {
  column_data <- venn_data[[name]]
  column_data <- column_data[!is.na(column_data)]
  column_data <- as.character(column_data)
  column_data <- gsub("\\s+", "", column_data)
  column_data
})

# Set vector names to column names
names(column_data) <- column_names

# Set circle colors for Venn diagram
circle_colors <- c("Amikacin sulfate" = "darksalmon",
                   "Ciprofloxacin" = "cadetblue3",
                   "Colistin" = "firebrick",
                   "Meropenem" = "orchid")

# Plot Venn diagram
venn_plot <- venn.diagram(column_data, filename = "Venn_Antibiotics.png",
                          fill = circle_colors, alpha = 0.45,
                          cat.col = circle_colors, cat.cex = 0.98,
                          cat.fontface = "bold", cex = 1.4,
                          fontfamily = "sans", cex.lab = 30,
                          cat.default.pos = "outer", cat.dist = c(0.23, 0.23, 0.1, 0.1), margin = 0.05)

# Print common and unique elements
common_all <- Reduce(intersect, column_data)
unique_amikacin <- setdiff(venn_data$`Amikacin sulfate`, union(venn_data$Ciprofloxacin, union(venn_data$Colistin, venn_data$Meropenem)))
unique_ciprofloxacin <- setdiff(venn_data$Ciprofloxacin, union(venn_data$`Amikacin sulfate`, union(venn_data$Colistin, venn_data$Meropenem)))
unique_colistin <- setdiff(venn_data$Colistin, union(venn_data$`Amikacin sulfate`, union(venn_data$Ciprofloxacin, venn_data$Meropenem)))
unique_meropenem <- setdiff(venn_data$Meropenem, union(venn_data$`Amikacin sulfate`, union(venn_data$Ciprofloxacin, venn_data$Colistin)))

# Create data frame for histogram plot
histogram_data <- read_excel("HeatmapVennHistogram.xlsx", sheet = "Histogram")

# Prepare data frame and convert to long format
data_long <- histogram_data %>%
  pivot_longer(-Condition) %>%
  separate(name, into = c("Type", "Antibiotic"), sep = "(?<=RxnNum|KEGGMetNum)(?=.)") %>%
  mutate(Type = ifelse(Type == "RxnNum", "Rxn", "KEGG"))

# Define MIC25 and MIC75 color scales
mic_colors <- c("MIC25" = "darksalmon", "MIC75" = "cadetblue3")

# Create bar plots for reactions and KEGG pathways
p_rxn <- ggplot(data_long %>% filter(Type == "Rxn"), aes(x = value, y = Condition, fill = Antibiotic)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Total number of reactions",
    y = "Condition",
    fill = "Antibiotic"
  ) +
  scale_fill_manual(values = mic_colors) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
    panel.grid.major.x = element_line(color = "white"),
    panel.grid.major.y = element_line(color = "grey88", size = 0.2, linetype = 2),
    axis.title.x = element_text(vjust = 0.5, margin = margin(t = 15, b = 15), face = "bold", size = 14),
    axis.title.y = element_text(vjust = 0.5, margin = margin(r = 15, l = 15), face = "bold", size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 95), breaks = seq(0, 95, by = 15)) +
  scale_y_discrete(limits = unique(data_long$Condition))

p_kegg <- ggplot(data_long %>% filter(Type == "KEGG"), aes(x = value, y = Condition, fill = Antibiotic)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Total number of KEGG pathways",
    y = "",
    fill = "Antibiotic"
  ) +
  scale_fill_manual(values = mic_colors) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.x = element_line(color = "white"),
    panel.grid.major.y = element_line(color = "grey88", size = 0.2, linetype = 2),
    axis.title.x = element_text(vjust = 0.5, margin = margin(t = 15, b = 15), face = "bold", size = 14),
    axis.title.y = element_text(vjust = 0.5, margin = margin(r = 15, l = 15), face = "bold", size = 14),
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical",
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_x_continuous(limits = c(0, 13), breaks = seq(0, 13, by = 5))

# Combine plots
combined_plot <- grid.arrange(p_rxn, p_kegg, ncol = 2, widths = c(4, 3))

# Save the combined plot
ggsave("Histogram_RxnPathNumbers.png", plot = combined_plot, width = 14, height = 6)
