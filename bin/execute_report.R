#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)

option_list = list(
    make_option(c("-r", "--report_file"), type="character", default=NULL, help="report rmarkdown file", metavar="character"),
    make_option(c("-rf", "--references"), type="character", default=NULL, help="references bibtex file", metavar="character"),
    make_option(c("-c","--css"), type ="character", default=NULL, help="css file for report formatting", metavar="character"),
    make_option(c("-l","--logo"), type="character", default=NULL, help="logo for report", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

wd=getwd()

rmarkdown::render(opt$report_file, output_file = "Bcellmagic_report.html", knit_root_dir = wd, output_dir = wd,
                    params = list(path_references = opt$references,
                                    path_css = opt$css,
                                    path_logo = opt$logo))
