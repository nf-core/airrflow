#!/usr/bin/env Rscript

library(rmarkdown)

args = commandArgs(trailingOnly=TRUE)

wd=getwd()

rmarkdown::render(args[1], output_file = "Bcellmagic_report.html", knit_root_dir = wd, output_dir = wd)
