#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Description: plot number of nuclei per probe.
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
parser = add_argument(parser, arg = 'dotsTable',
	help = 'Path to dots table.')
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

if ( !file.exists(dotsTable) )
	stop(sprintf("dots table not found: %s", dotsTable))

if ( !dir.exists(outputFolder) )
	stop(sprintf("output folder not found: %s", outputFolder))

if ( 0 != nchar(suffix) ) {
	if ( !grepl("^\\.", suffix) ) suffix = paste0(".", suffix)
}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "",
	initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# FUNCTIONS ====================================================================

source(file.path(script.basename, "common.functions.R"))

# RUN ==========================================================================

tdot = read.delim(dotsTable, as.is = T, header = T)

data = do.call(rbind, by(tdot, tdot$probe_label, FUN = function(pt) {
	data.frame(
		probe = pt$probe_label[1],
		chrom = unlist(strsplit(pt$probe_label[1], '.', fixed = T))[1],
		nnuclei = length(unique(pt$nuID)),
		ndots = nrow(pt),
		stringsAsFactors = F
	)
}))
data = add_chrom_ID(data)

write.table(data, file.path(outputFolder,
	sprintf("nuclei_by_probe%s.tsv", suffix)),
	quote = F, row.names = F, sep = "\t")

nLowProbes = sum(data$nnuclei < 100)

pdf(sprintf("%s/nuclei_by_probe%s.pdf", outputFolder, suffix),
	width = 25, height = 10)

ggplot(data, aes(
		x = reorder(probe, chromID),
		y = nnuclei,
		fill = reorder(chrom, chromID)
	)
	) + geom_bar(stat = 'identity'
	) + xlab("Probe label") + ylab("Number of nuclei"
	) + geom_hline(yintercept = 100, linetype = "dashed", color = "black"
	) + ggtitle(sprintf("%d probes with less than 100 nuclei.", nLowProbes)
	) + theme_bw(
	) + theme(
		axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
	)

ggplot(data, aes(
		x = reorder(probe, chromID),
		y = ndots,
		fill = reorder(chrom, chromID)
	)
	) + geom_bar(stat = 'identity'
	) + xlab("Probe label") + ylab("Number of dots"
	) + theme_bw(
	) + theme(
		axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
	)

graphics.off()

# END --------------------------------------------------------------------------

################################################################################
