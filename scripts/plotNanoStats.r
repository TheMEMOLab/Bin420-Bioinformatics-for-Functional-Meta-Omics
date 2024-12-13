# Load required libraries
library(tidyverse)

# Function to read and process data from a file
read_summary <- function(file, label) {
  read_lines(file) %>%
    # Split each line into name and value
    str_split_fixed(":\\s+", 2) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    # Name the columns
    rename(Metric = V1, Value = V2) %>%
    # Convert Value column to numeric
    mutate(Value = as.numeric(gsub(",", "", Value)),
           File = label) # Add file label
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure two arguments are provided
if (length(args) != 2) {
  stop("Usage: Rscript plot_metrics.R <file_A> <file_B>")
}

# File paths from arguments
file_a <- args[1]
file_b <- args[2]

# Extract file labels (names without .txt)
label_a <- tools::file_path_sans_ext(basename(file_a))
label_b <- tools::file_path_sans_ext(basename(file_b))

# Read data from files with labels
data_a <- read_summary(file_a, label_a)
data_b <- read_summary(file_b, label_b)

# Combine data from both files
combined_data <- bind_rows(data_a, data_b)

# Filter and scale the required metrics
filtered_data <- combined_data %>%
  filter(Metric %in% c("Mean read length", "Mean read quality",
                       "Median read length", "Number of reads", "Total bases")) %>%
  # Scale metrics
  mutate(Value = case_when(
           Metric == "Mean read length" ~ Value / 1e3,
           Metric == "Median read length" ~ Value / 1e3,
           Metric == "Total bases" ~ Value / 1e9,
           Metric == "Number of reads" ~ Value / 1e6,
           TRUE ~ Value
         ),
         Metric = case_when(
           Metric == "Mean read length" ~ "Mean read length (Thousands)",
           Metric == "Median read length" ~ "Median read length (Thousands)",
           Metric == "Total bases" ~ "Total bases (Billions)",
           Metric == "Number of reads" ~ "Number of reads (Millions)",
           TRUE ~ Metric
         ))

# Check unique File values to debug issues
print(unique(filtered_data$File))

# Plot the data
plot <- ggplot(filtered_data, aes(x = Metric, y = Value, fill = File)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of Selected Metrics",
       x = "Metric",
       y = "Value",
       fill = "File") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Use File names directly for colors
  scale_fill_manual(values = setNames(c("#1f78b4", "#33a02c"), c(label_a, label_b)))

# Save the plot to a file
output_file <- "comparison_NanoStats.pdf"
ggsave(output_file, plot, width = 10, height = 6)
cat("Plot saved to", output_file, "\n")
