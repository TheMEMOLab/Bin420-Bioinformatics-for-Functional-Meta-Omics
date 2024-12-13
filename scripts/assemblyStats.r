library(tidyverse)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript script_name.R <input_file>")
}
input_file <- args[1]

# Read data
data <- read_delim(
  input_file,
  delim = "\t",
  col_names = TRUE
)

# Extract filename base and scale values
processed_data <- data %>%
  mutate(Legend = basename(filename)) %>%
  mutate(across(
    where(is.numeric),
    ~ case_when(
      cur_column() == "number" ~ . / 1e6,        # Divide `number` by millions
      cur_column() == "total_length" ~ . / 1e6,  # Divide `total_length` by millions
      TRUE ~ . / 1e3                            # Others by thousands
    )
  )) %>%
  pivot_longer(
    cols = -c(filename, Legend),
    names_to = "Metric",
    values_to = "Value"
  )

# Find the max total_length (in millions) to define scaling
max_total_length <- max(processed_data$Value[processed_data$Metric == "total_length"], na.rm = TRUE)

# We want to bring total_length onto a 0–100 scale for the primary axis
scaling_factor <- 100 / max_total_length

# Apply scaling to total_length values only
processed_data <- processed_data %>%
  mutate(Value = case_when(
    Metric == "total_length" ~ Value * scaling_factor,
    TRUE ~ Value
  ))

# Now, when we plot:
# - The primary axis will show all metrics (including total_length) in a 0–100 range.
# - The secondary axis will be transformed back to the original scale by dividing by the scaling factor.
#   This means if total_length was originally in millions, the secondary axis will reflect that.

plot <- ggplot(processed_data, aes(x = Metric, y = Value, fill = Legend)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of Metrics",
       x = "Metric",
       y = "Kb",
       fill = "File") +
  scale_y_continuous(
    limits = c(0, 100),
    sec.axis = sec_axis(
      trans = ~ . / scaling_factor,        # Convert back to original millions
      name = "Mb"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(angle = 90)
  )

output_file <- "AssemblyStats.pdf"
ggsave(output_file, plot, width = 12, height = 8)
cat("Plot saved as", output_file, "\n")
