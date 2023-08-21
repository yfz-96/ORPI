library(RColorBrewer)
#I/O---------------------------------------------------
table <- read.table("../data/Abundance_Stat.all.xls", sep = "\t", header = T, comment.char = "")
rownames(table) <- table$Specie
table <- table[, 8:ncol(table)]
table <- table[, order(colnames(table))]
meta <- read.table("../data/meta.txt", sep = "\t", header = T, row.names = 1)
meta <- meta[order(rownames(meta)), ]
pathogen_list <- read.table("../respiratory_pathogens.txt", sep = "\t", header = F)[,1]
identical(colnames(table), rownames(meta))

output_dir = '../results/orpi_of_sample_between_group_test_results/'
if(! dir.exists(output_dir)) {
  dir.create(output_dir)
}

#ORPI
cal_orpi <- function(pathogen_table) {
  return (colSums(pathogen_table))
}

pathogen_table <- subset(table, rownames(table) %in% pathogen_list)
orpi <- cal_orpi(pathogen_table)
orpi<- merge(orpi, meta, by = 0, all = T)
write.table(orpi, paste0(output_dir, "ORPI.txt"), sep = "\t", quote = F, row.name = F)

x <- subset(orpi, Group == "Clean_Prosthesis")$x
y <- subset(orpi, Group == "Unclean_Prosthesis")$x

shapiro.test(x)
shapiro.test(y)
quantile_orpi_clean_Prosthesis <- quantile(x)
round(quantile_orpi_clean_Prosthesis, 5)
quantile_orpi_unclean_Prosthesis <- quantile(y)
round(quantile_orpi_unclean_Prosthesis, 5)

wilcoxon_result <- wilcox.test(x, y)
wilcoxon_result

colnames(orpi) <- c("Pathogen", "value", "Group", "Group2")
plot <- ggplot(orpi, aes(x = Group, y = value)) +
  geom_boxplot() +
  scale_colour_manual(values = brewer.pal(7, "Pastel1")[c(3, 1)]) +
  geom_point(aes(color = Group), shape = 16, size = 4) +
  xlab("Group") + 
  ylab("ORPI") + 
  geom_text(x = 1.5, y = max(orpi$value) * 0.9, label = paste("p =", round(wilcoxon_result$p.value, 3))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
plot
ggsave(filename = paste0(output_dir, "orpi_boxplot.pdf"), plot = plot, width = 3.5, height = 4)

clean_prosthesis_orpi <- orpi %>% filter(Group == "Clean_Prosthesis")
clean_prosthesis_orpi <- clean_prosthesis_orpi[order(-clean_prosthesis_orpi$value), ]
clean_prosthesis_orpi[1:5, ]
unclean_prosthesis_orpi <- orpi %>% filter(Group == "Unclean_Prosthesis")
unclean_prosthesis_orpi <- unclean_prosthesis_orpi[order(-unclean_prosthesis_orpi$value), ]
unclean_prosthesis_orpi[1:5, ]
