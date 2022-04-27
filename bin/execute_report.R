#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)

option_list = list(
    make_option(c("-r", "--report_file"), type="character", default=NULL, help="report rmarkdown file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

wd=getwd()

rmarkdown::render(opt$report_file, output_file = "Bcellmagic_report.html", knit_root_dir = wd, output_dir = wd)
