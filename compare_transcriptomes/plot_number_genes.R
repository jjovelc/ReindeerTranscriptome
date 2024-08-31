# Load necessary libraries
library(ggplot2)
library(reshape2)

# Data for the first series
hybrid_t <- c(17415, 17949, 17722, 17919)
names(hybrid_t) <- c("R_BM_1.genes", "R_BM_2.genes", "R_BM_A.genes", "R_BM_S.genes")

# Data for the second series
Btaurus_t <- c(13981, 14245, 14239, 14237)
names(Btaurus_t) <- c("R_BM_1.genes", "R_BM_2.genes", "R_BM_A.genes", "R_BM_S.genes")

# Combine data into a data frame
data <- data.frame(
  Library = names(first_series),
  Hybrid_t = hybrid_t,
  Btaurus_t = Btaurus_t
)

# Melt the data frame for ggplot2
data_melted <- melt(data, id.vars = "Library")

# Create the paired barplot
ggplot(data_melted, aes(x = Library, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Genes per Library",
       x = "Library",
       y = "Number of Genes",
       fill = "Experiment") +
  theme_minimal()

# Save the plot if needed
# ggsave("paired_barplot.png")