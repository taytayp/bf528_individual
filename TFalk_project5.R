## Taylor Falk
## tfalk@bu.edu
## BF528 - Project 5

library(ggpubr)
library(knitr)
library(tidyverse)
options(scipen = 999)

### Step 5 - Programmer
# create histogram of cufflinks output
genes <- read.table(paste0("/projectnb2/bf528/users/hedgehog/",
                           ".tay/P0_1_cufflinks/genes.fpkm_tracking"), header = T)
sprintf("%i genes have non-zero FPKM values, of %i total genes (%f)", 
        length(which(genes$FPKM != 0)),
        length(genes$FPKM),
        length(which(genes$FPKM != 0))/length(genes$FPKM))

ggplot(genes[genes$FPKM > 0,], aes(FPKM)) + 
  geom_histogram(bins = 90, fill="#7D58F0", color="darkgray") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,10,100,1000,10000,100000), 
                     trans="log1p", expand=c(0,0)) +
  ggtitle("Histogram of non-zero FPKM counts in P0_1") +
  xlab("FPKM, log-transformed") +
  ylab("Counts") +
  theme_light() 

### Step 6 - Analyst
step_six <- function(cuffDiffFile) {
  df <- read.table(cuffDiffFile, header = T)
  df <- df[order(df$q_value),]
  print(kable(df[1:10,c("gene", "log2.fold_change.", "p_value", "q_value")], 
              format = "latex", row.names = F, 
              col.names = c("Gene", "Log2 fold change", "p-value", "q-value")))
  
  p1 <- ggplot(df, aes(log2.fold_change.)) + 
    geom_histogram(bins=90, fill = "#3CF086", color="darkgray") +
    scale_y_log10() +
    ggtitle("Histogram of log2 fold change values") +
    xlab("Log2 fold change") +
    ylab("Counts, log10-transformed") +
    theme_light() 
  
  p2 <- ggplot(df[df$significant == "yes",], aes(log2.fold_change.)) + 
    geom_histogram(bins=90, fill = "#3CF086", color="darkgray") +
    scale_y_log10() +
    xlim(c(-7.5, 7.5)) +
    ggtitle("Histogram of significant log2 fold change values") +
    xlab("Log2 fold change") +
    ylab("Counts, log10-transformed") +
    theme_light() 
  
  updf <- df[df$log2.fold_change. > 2 & df$significant == "yes",]
  downdf <- df[df$log2.fold_change. < -2 & df$significant == "yes",]
  print(sprintf("%i genes were found to be significant, with %i up and %i down reg.",
                dim(df[df$significant == "yes",])[1],
                dim(df[df$significant == "yes" & df$log2.fold_change. > 0,])[1],
                dim(df[df$significant == "yes" & df$log2.fold_change. < 0,])[1]))
  print(sprintf("There are %i up-regulated genes and %i down-regulated genes (%i total).", 
          dim(updf)[1], dim(downdf)[1], length(df$test_id)))
  write.csv(updf$gene, file="/projectnb2/bf528/users/hedgehog/.tay/upreg.csv", 
            row.names = F)
  write.csv(downdf$gene, file="/projectnb2/bf528/users/hedgehog/.tay/downreg.csv",
             row.names = F)
  
  ggarrange(p1, p2,
            labels = c("A", "B"),
            ncol = 1, nrow = 2)
}

step_six("/projectnb2/bf528/users/hedgehog/.tay/cuffdiff_out/gene_exp.diff")
