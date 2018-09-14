#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.0.1
# Description: plot distribution of nuclear features per reorder(dataset, date) and cell type,
# 			   as a manual sanity check.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(ggplot2))

# INPUT ========================================================================

# Create arguent parser
scriptName = "nucleiDatasetStats.R"
parser = arg_parser('Generate density plots per dataset for nuclear features',
	name = scriptName)

# Define mandatory arguments
parser = add_argument(parser, arg = 'nucleiTable',
	help = 'Path to nuclei table.')
parser = add_argument(parser, arg = 'compTable',
	help = 'Path to compartment table.')
parser = add_argument(parser, arg = 'outputFolder',
	help = 'Path to output folder.')

# Define elective arguments
parser = add_argument(parser, arg = '--suffix', short = '-s',
	help = 'Suffix for output pdf.', default = "")

# Version argument
version_flag = "0.0.1"
parser = add_argument(parser, arg = '--version', short = '-V',
	help = 'Print version and quit.', flag = T)

args = commandArgs(trailingOnly=TRUE)
if ( "--version" %in% args ) {
	cat(sprintf("%s v%s\n", scriptName, version_flag))
	quit()
}

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)], warn.conflicts = F)

# RUN ==========================================================================

tnuc = read.delim(nucleiTable, as.is = T, header = T)
tcom = read.delim(compTable, as.is = T, header = T)

if ( !dir.exists(outputFolder) )
	stop(sprintf("output folder not found: %s", outputFolder))

if ( 0 != nchar(suffix) ) {
	if ( !grepl("^\\.", suffix) ) suffix = paste0(".", suffix)
}
pdf(sprintf("%s/nuclei_by_dataset%s.pdf", outputFolder, suffix),
	width = 20, height = 10)

ggplot(tnuc, aes(x = flat_size, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("Area in Z projection [px]"
	)

ggplot(tnuc, aes(x = size, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("Volume [vx]"
	)

ggplot(tnuc, aes(x = surf, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("Surface [px]"
	)

ggplot(tnuc, aes(x = sumI, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("Intensity integral over volume [a.u.]"
	)

ggplot(tnuc, aes(x = flat_sumI, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("Intensity integral in Z-projection [a.u.]"
	)

ggplot(tnuc, aes(x = meanI, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("Average intensity integral over volume [a.u.]"
	)

ggplot(tnuc, aes(x = shape, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("Sphericity [a.u.]"
	)

ggplot(tcom, aes(x = a, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("1st semi-axis [px]"
	)

ggplot(tcom, aes(x = b, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("2nd semi-axis [px]"
	)

ggplot(tcom, aes(x = c, color = reorder(dataset, date))
	) + geom_density() + theme_bw(
	) + facet_wrap(~cell_type
	) + scale_x_log10(
	) + guides(color = guide_legend(title = "dataset")
	) + xlab("3rd semi-axis (Z) [px]"
	)

graphics.off()

# END --------------------------------------------------------------------------

################################################################################
