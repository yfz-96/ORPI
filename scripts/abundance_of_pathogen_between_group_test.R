source("BetweenGroup.test.R")
#I/O---------------------------------------------------
table <- read.table("../data/Abundance_Stat.all.xls", sep = "\t", header = T, comment.char = "")
rownames(table) <- table$Specie
table <- table[, 8:ncol(table)]
table <- table[, order(colnames(table))]
meta <- read.table("../data/meta.txt", sep = "\t", header = T, row.names = 1)
meta <- meta[order(rownames(meta)), ]
pathogen_list <- read.table("../respiratory_pathogens.txt", sep = "\t", header = F)[,1]
identical(colnames(table), rownames(meta))

output_dir = '../results/abundance_of_pathogen_between_group_test_results/'
if(! dir.exists(output_dir)) {
  dir.create(output_dir)
}

#-------------------------------
# install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2", "dplyr", "ROCR", "doParallel", "parallel", "tidyverse", "foreach","viridis", "doMC", "plyr", "rlang")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
registerDoParallel(cores=4)

# identify differentially abundant features
t_table <- t(table)

x <- t_table
y <- factor(meta$Group)

sp_out <- BetweenGroup.test(x, y, clr_transform=TRUE, positive_class = "Clean_Prosthesis", q_cutoff = 0.05)
sp_out_sig <- sp_out %>% filter(IfSig=="Sig")
pathogen_out_sig <- sp_out_sig %>% filter(rownames(sp_out_sig) %in% pathogen_list)
pathogen_out_sig <- subset(sp_out_sig, rownames(sp_out_sig) %in% pathogen_list)
write.table(pathogen_out_sig, paste0(output_dir, "wilcox_out_sig_pathogens.tsv"), sep = "\t", quote = F, row.names = T, col.names = NA)
