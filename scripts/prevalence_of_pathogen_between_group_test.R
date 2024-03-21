p <- c("RColorBrewer", "reshape2", "ggplot2", "scales")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
## clean R environment
rm(list = ls())
setwd('./')

table <- read.table("../data/Abundance_Stat.all.xls", sep = "\t", header = T, comment.char = "")
rownames(table) <- table$Specie
table <- table[, 8:ncol(table)]
table <- table[, order(colnames(table))]
meta <- read.table("../data/meta.txt", sep = "\t", header = T, row.names = 1)
meta <- meta[order(rownames(meta)), ]
pathogen_list <- read.table("../respiratory_pathogens.txt", sep = "\t", header = F)[,1]
identical(colnames(table), rownames(meta))

output_dir = '../results/prevalence_of_pathogen_between_group_test_results/'
if(! dir.exists(output_dir)) {
  dir.create(output_dir)
}

count_positive <- function(row) {
  sum(row > 0)
}

pathogen_table <- subset(table, rownames(table) %in% pathogen_list)

Prevalence_of_ANY_type_of_orp  <- apply(pathogen_table, 1, count_positive) / ncol(pathogen_table)

identical(colnames(pathogen_table), rownames(meta))
clean_Prosthesis_pathogen_table <- pathogen_table[, which(meta[match(colnames(pathogen_table), rownames(meta)), "Group"] == "Clean_Prosthesis")]
Prevalence_of_clean_Prosthesis_of_orp  <- apply(clean_Prosthesis_pathogen_table, 1, count_positive) / ncol(clean_Prosthesis_pathogen_table)
print(Prevalence_of_clean_Prosthesis_of_orp)

unclean_Prosthesis_pathogen_table <- pathogen_table[, which(meta[match(colnames(pathogen_table), rownames(meta)), "Group"] == "Unclean_Prosthesis")]
Prevalence_of_unclean_Prosthesis_of_orp  <- apply(unclean_Prosthesis_pathogen_table, 1, count_positive) / ncol(unclean_Prosthesis_pathogen_table)
print(Prevalence_of_clean_Prosthesis_of_orp)

Prevalence_of_orp <- data.frame(all_type = Prevalence_of_ANY_type_of_orp, Clean_Prosthesis = Prevalence_of_clean_Prosthesis_of_orp, Unclean_Prosthesis = Prevalence_of_unclean_Prosthesis_of_orp)
print(Prevalence_of_orp)
write.table(Prevalence_of_orp, paste0(output_dir, "Prevalence_of_orp.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

x <- Prevalence_of_orp$Clean_Prosthesis
y <- Prevalence_of_orp$Unclean_Prosthesis
quantile_prevalence_clean_Prosthesis <- quantile(x)
round(quantile_prevalence_clean_Prosthesis, 5)
quantile_prevalence_unclean_Prosthesis <- quantile(y)
round(quantile_prevalence_unclean_Prosthesis, 5)
shapiro.test(x)
shapiro.test(y)

result <- wilcox.test(x, y, paired = TRUE)
print(result)

Prevalence_of_orp$Pathogen <- rownames(Prevalence_of_orp)
Prevalence_df <- melt(Prevalence_of_orp, id.vars = "Pathogen", variable.name = "type", value.name = "value")

palette1 <- brewer.pal(9, "Pastel1")
palette2 <- brewer.pal(8, "Set3")
palette3 <- brewer.pal(8, "Set2")
palette4 <- brewer.pal(8, "Set1")


# 将多个调色板合并成一个向量
color_vector <- c(palette1, palette2, palette3, palette4)

Prevalence_df <- subset(Prevalence_df, type != "all_type")

plot <- ggplot(Prevalence_df, aes(x = type, y = value)) +
  geom_boxplot() +
  scale_colour_manual(values = color_vector) +
  geom_point(aes(color = Pathogen), shape = 16, size = 4) +
  xlab("Group") + 
  ylab("Prevalence of opportunistic respiratory pathogen") + 
  geom_text(x = 1.5, y = 0.4, label = paste("p =", round(result$p.value, 3))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
plot
ggsave(filename = paste0(output_dir, "Prevalence_of_orp.pdf"), plot = plot, width = 7, height = 6)

Prevalence_of_orp <- Prevalence_of_orp[order(-Prevalence_of_orp$Clean_Prosthesis), ]
Prevalence_of_orp[1:5, c("Clean_Prosthesis", "Pathogen")]
Prevalence_of_orp <- Prevalence_of_orp[order(-Prevalence_of_orp$Unclean_Prosthesis), ]
Prevalence_of_orp[1:5, c("Unclean_Prosthesis", "Pathogen")]