getwd()
setwd("~/Desktop/data science /DEseq2 project")

#after DESeq2 RNAseq analysis
# MA plot
plotMA(res0.01)

#custom data visualization of MA plot using ggplot2
#MA plot
library(ggplot2)
res0.01$significant <- ifelse(res0.01$padj < 0.01, "Significant", "Not Significant")
ggplot(res0.01, aes(x = baseMean, y = log2FoldChange, colour = significant)) +
  geom_point(alpha = 0.8, size = 1) +
  scale_x_log10() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(title = "MA Plot", x = "Mean expression (log scale)", y = "Log2 fold change") +
  scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "red")) +
  theme_minimal()

#Volcano plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = (abs(log2FoldChange) > 2 & padj < 0.05))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("red", "blue")) +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)",
       color = "Significant") +
  theme_minimal()

