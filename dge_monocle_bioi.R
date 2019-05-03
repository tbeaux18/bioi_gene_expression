#!/usr/bin/env Rscript
library("monocle")
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

umi.keep <- rowSums(cpm(umi.count) > 1) >= 4
umi.count <- umi.count[umi.keep,]

read.keep <- rowSums(cpm(read.count) > 1) >= 4
read.count <- read.count[read.keep,]

pd <- data.frame(group = cell.data$condition)
rownames(pd) <- rownames(cell.data)
pd <- new("AnnotatedDataFrame", data = pd)

umi.obj <- newCellDataSet(as.matrix(umi.count), phenoData = pd, expressionFamily = negbinomial.size())
umi.obj <- estimateSizeFactors(umi.obj)
umi.obj <- estimateDispersions(umi.obj)
umi.res <- differentialGeneTest(umi.obj, fullModelFormulaStr = "~group")
umi.sig.genes <- subset(umi.res, qval < 0.08)
write.csv(as.data.frame(umi.sig.genes), 
          file="umi_sig_deg_monocle.csv")
read.obj <- newCellDataSet(as.matrix(read.count), phenoData = pd, expressionFamily = negbinomial.size())
read.obj <- estimateSizeFactors(read.obj)
read.obj <- estimateDispersions(read.obj)
read.res <- differentialGeneTest(read.obj, fullModelFormulaStr = "~group")
read.sig.genes <- subset(read.res, qval < 0.08)
write.csv(as.data.frame(read.sig.genes), 
          file="read_sig_deg_monocle.csv")
