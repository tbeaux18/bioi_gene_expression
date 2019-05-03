#!/usr/bin/env Rscript
library("MAST")
library("dplyr")
library("plyr")
library("data.table")
library("GenomicRanges")

FCTHRESHOLD <- log2(1.5)
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

umi.keep <- rowSums(edgeR::cpm(umi.count) > 1) >= 4
umi.count <- umi.count[umi.keep,]

read.keep <- rowSums(edgeR::cpm(read.count) > 1) >= 4
read.count <- read.count[read.keep,]

umi.count.cpm <- edgeR::cpm(umi.count)
read.count.cpm <- edgeR::cpm(read.count)

umi.log.count <- log(umi.count.cpm + 1) / log(2)
umi.fData <- data.frame(names = rownames(umi.log.count))
rownames(umi.fData) <- rownames(umi.log.count);
umi.cData <- data.frame(cond = exp.group)
rownames(umi.cData) <- colnames(umi.log.count)

umi.obj <- FromMatrix(as.matrix(umi.log.count), umi.cData, umi.fData)
colData(umi.obj)$cngeneson <- scale(colSums(assay(umi.obj) > 0))
cond <- factor(colData(umi.obj)$cond)
umi.zlmCond <- zlm(~ cond + cngeneson, umi.obj)
umi.summaryCond <- summary(umi.zlmCond, doLRT=TRUE)
umi.summaryDt <- umi.summaryCond$datatable
# umi.summaryDt <- as.data.frame(umi.summaryCond)
# print(umi.summaryCond, n=4)

umi.fcHurdle <- merge(umi.summaryDt[contrast=='condstarve' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  umi.summaryDt[contrast=='condstarve' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
umi.fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
umi.fcHurdleSig <- merge(umi.fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(umi.obj)), by='primerid')
setorder(umi.fcHurdleSig, fdr)


read.log.count <- log(read.count.cpm + 1) / log(2)
read.fData <- data.frame(names = rownames(read.log.count))
rownames(read.fData) <- rownames(read.log.count);
read.cData <- data.frame(cond = exp.group)
rownames(read.cData) <- colnames(read.log.count)

read.obj <- FromMatrix(as.matrix(read.log.count), read.cData, read.fData)
colData(read.obj)$cngeneson <- scale(colSums(assay(read.obj) > 0))
cond <- factor(colData(read.obj)$cond)
read.zlmCond <- zlm(~ cond + cngeneson, read.obj)
read.summaryCond <- summary(read.zlmCond, doLRT=TRUE)
read.summaryDt <- umi.summaryCond$datatable

read.fcHurdle <- merge(read.summaryDt[contrast=='condstarve' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  read.summaryDt[contrast=='condstarve' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
read.fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
read.fcHurdleSig <- merge(read.fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(read.obj)), by='primerid')
setorder(read.fcHurdleSig, fdr)