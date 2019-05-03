#!/usr/bin/env Rscript
library("edgeR")
library("dplyr")

####################### Data #######################
all.counts <- readRDS("/Users/tim.baker/Documents/loyola/FatFlies_scRNA/docker_run_20190425/zUMIs_output/expression/dnpf_scrna_run.dgecounts.rds")

# cell data supplied from SampleSheet.csv at beginning of pipeline
cell.data <- read.csv("/Users/tim.baker/Documents/loyola/FatFlies_scRNA/cell_data.csv", header=TRUE, sep=",")
cell.data$col_name <- paste(cell.data$gfp_state, cell.data$condition, cell.data$cell_id, sep=".")
rownames(cell.data) <- cell.data$col_name
cell.data$comb_condition <- paste(cell.data$condition, cell.data$gfp_state, sep=".")

# making count tables
umi.count <- as.matrix(all.counts$umicount$inex$all)
read.count <- as.matrix(all.counts$readcount$inex$all)

umi.count <- umi.count[,cell.data$barcode_sequence]
colnames(umi.count) <- cell.data[["col_name"]]

read.count <- read.count[,cell.data$barcode_sequence]
colnames(read.count) <- cell.data[["col_name"]]

# removing to save space
rm(all.counts)

# subsetting data all the way down testing group
cell.data <- cell.data %>% dplyr::filter(gfp_state=='pos')
cell.data <- as.data.frame(cell.data)
rownames(cell.data) <- cell.data$col_name

# vector of column names to subset the matrix
exp.names <- cell.data$col_name
exp.group <- cell.data$condition

# subsetting the matrix
umi.count <- umi.count[,exp.names]
read.count <- read.count[,exp.names]



umi.y <- DGEList(counts=umi.count, group=exp.group)
umi.keep <- rowSums(cpm(umi.y) > 1) >= 4
umi.y <- umi.y[umi.keep, , keep.lib.sizes=FALSE]
umi.y <- calcNormFactors(umi.y, method="TMMwzp")
design <- model.matrix(~exp.group)
umi.y <- estimateDisp(umi.y, design)
umi.fit <- glmFit(umi.y, design)
umi.lrt <- glmLRT(umi.fit, coef=2)
umi.glm_results <- topTags(umi.lrt, n = nrow(umi.count), sort.by = "none")
umi.sig_results <- umi.glm_results[umi.glm_results$table$FDR < 0.08,]
write.csv(as.data.frame(umi.sig_results), 
          file="umi_sig_deg_edgeR.csv")
read.y <- DGEList(counts=read.count, group=exp.group)
read.keep <- rowSums(cpm(read.y) > 1) >= 4
read.y <- read.y[read.keep, , keep.lib.sizes=FALSE]
read.y <- calcNormFactors(read.y, method="TMMwzp")
design <- model.matrix(~exp.group)
read.y <- estimateDisp(read.y, design)
read.fit <- glmFit(read.y, design)
read.lrt <- glmLRT(read.fit, coef=2)
read.glm_results <- topTags(read.lrt, n = nrow(read.count), sort.by = "none")
read.sig_results <- read.glm_results[read.glm_results$table$FDR < 0.08,]
write.csv(as.data.frame(read.sig_results), 
          file="read_sig_deg_edgeR.csv")