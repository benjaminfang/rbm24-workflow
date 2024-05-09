#! /usr/bin/env Rscript

args = commandArgs()
library(celda)
library(singleCellTK)
library(Matrix)

sce <- importCellRangerV3(sampleDirs=args[6], dataType="filter")
sce.raw <- importCellRangerV3(sampleDirs=args[6], dataType="raw")
sce <- decontX(sce, background=sce.raw)
mtx <- decontXcounts(sce)
mtx = Matrix(mtx, sparse=TRUE)
writeMM(mtx, args[7])
