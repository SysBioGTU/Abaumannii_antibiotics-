# Load necessary libraries
library(logisticPCA)
library(ggplot2)
library(rARPACK)
library(readxl)
library(rgl)
library(shiny)
library(plotly)

# Read data from Excel file
data <- read_excel("logisticSVD.xlsx")

# Extract column names and create a binary matrix
column_names <- colnames(data)[-1]
binary_matrix <- as.matrix(data[, -1])

# Identify and remove columns with constant values
constant_columns <- apply(binary_matrix, 2, function(x) length(unique(x)) == 1)
updated_matrix <- binary_matrix[, !constant_columns]
updated_column_names <- column_names[!constant_columns]

# Create a new data frame with updated column names
new_data <- data.frame(Condition = data$Condition, updated_matrix)
colnames(new_data)[-1] <- updated_column_names
row.names(new_data) <- new_data$Condition
new_data$Condition <- NULL

# Modify the 'SubString' column to combine specific substrings
new_data$SubString <- gsub("^\\d+_([^_]+)_.*", "\\1", rownames(new_data))
new_data$SubString <- gsub("^(AMI\\d+).*", "Amikacin sulfate", new_data$SubString)
new_data$SubString <- gsub("^(CIP\\d+).*", "Ciprofloxacin", new_data$SubString)
new_data$SubString <- gsub("^(COL\\d+).*", "Colistin", new_data$SubString)
new_data$SubString <- gsub("^(MERO\\d+).*", "Meropenem", new_data$SubString)

# Define color mapping for each substring
color_mapping <- c("Amikacin sulfate" = "#008080", 
                   "Ciprofloxacin" = "#DAA520",
                   "Colistin" = "#FF6347",  # Tomato color
                   "Meropenem" = "#800080")  # Purple

# Perform logisticSVD with k=3
logsvd_model <- logisticSVD(new_data[, -ncol(new_data)], k = 3)

# Define UI elements for Shiny app
ui <- fluidPage(
  titlePanel("Antibiotics"),
  fluidRow(
    column(width = 10,
           plotly::plotlyOutput("plotlyPlot", height = "1500px", width = "2000px")
    )
  )
)

# Define server logic for Shiny app
server <- function(input, output, session) {
  output$plotlyPlot <- renderPlotly({
    # Create a 3D plot using plotly with custom colors and legend labels
    p <- plot_ly(
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 9),
      text = new_data$SubString
    )
    
    # Add traces for each level with desired colors
    unique_levels <- unique(new_data$SubString)
    for (level in unique_levels) {
      subset_data <- new_data[new_data$SubString == level, ]
      p <- add_trace(
        p,
        x = logsvd_model$A[new_data$SubString == level, 1],
        y = logsvd_model$A[new_data$SubString == level, 2],
        z = logsvd_model$A[new_data$SubString == level, 3],
        color = factor(subset_data$SubString),
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 12.5, color = color_mapping[level]),
        text = subset_data$SubString,
        name = level,
        showlegend = TRUE 
      )
    }
    
    # Update layout to show the legend
    p <- layout(
      p,
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      ),
      legend = list(
        x = 0.5,
        y = 1.0,
        xanchor = "center",
        yanchor = "bottom",
        orientation = "h",
        title = "Strains",
        traceorder = "normal",
        bgcolor = "rgba(255, 255, 255, 0)",
        bordercolor = "#FFFFFF",
        borderwidth = 2,
        itemsizing = "constant",
        itemwidth = 35
      ),
      margin = list(l = 100, r = 600, b = 500, t = 100)
    )
    
    return(p)
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
