#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.0.1
# Description:
# 	- Plot distribution of intensity interal over volume and flattened area.
# 	- Fit sum of Gaussians and select k*sigma range around the major Peak.
# 	- Plot scatterplot of the two and color the selected ones differently.
# 	- Subselect nuclei and filter dot table accordingly.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))

# INPUT ========================================================================

# Create arguent parser
scriptName = "nucleiSelectPopulation.R"
parser = arg_parser(paste0("Subselect nuclei based on area in Z-projection",
	" and intensity integral. Specifically, fits sum of Gaussians to the",
	" distribution, then select an interval of +-k*sigma around the mean",
	" of the first Gaussian. This selection is applied to the following",
	" tables: nuclei, nuclear compartments, and dots."), name = scriptName)

# Define mandatory arguments
parser = add_argument(parser, arg = 'dotsTable',
	help = 'Path to dots table.')
parser = add_argument(parser, arg = 'nucleiTable',
	help = 'Path to nuclei table.')
parser = add_argument(parser, arg = 'compTable',
	help = 'Path to compartment table.')
parser = add_argument(parser, arg = 'outputFolder',
	help = 'Path to output folder.')

# Define elective arguments
parser = add_argument(parser, arg = '--elective', short = '-e', type = class(0),
	help = 'Elective argument.', default = 1, nargs = 1)

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
attach(p['' != names(p)])

# RUN ==========================================================================

tdot = read.delim(dotTable, as.is = T, header = T)
tnuc = read.delim(nucleiTable, as.is = T, header = T)
tcom = read.delim(compTable, as.is = T, header = T)

# END --------------------------------------------------------------------------

################################################################################
