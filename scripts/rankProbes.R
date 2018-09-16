#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Description: rank probes per cell type.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(ggplot2))
suppressMessages(library(readr))

# INPUT ========================================================================

# Create arguent parser
scriptName = "rankProbes.R"
parser = arg_parser('Rank probes, per cell type.',
	name = scriptName)

# Define mandatory arguments
parser = add_argument(parser, arg = 'dotsTable', type = class(""),
	help = 'Path to dots table.')
parser = add_argument(parser, arg = 'probesBed', type = class(""),
	help = 'Path to probes bed file. No header, track line.')
parser = add_argument(parser, arg = 'outputFolder', type = class(""),
	help = 'Path to output folder.')
parser = add_argument(parser, arg = '--binSize', short = '-b', type = class(0),
	help = 'Space-separated sizes for window around probe center.',
	nargs = Inf)

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
	stop(sprintf("dots table not found: %s\n", dotsTable))
if ( !file.exists(probesBed) )
	stop(sprintf("probes bed file not found: %s\n", probesBed))
if ( !file.exists(outputFolder) )
	cat(sprintf("Output folder not found, created: %s\n", outputFolder))
	dir.create(outputFolder, showWarnings = F)

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

bpro = read_probe_bed(probesBed)
tdot = read_dot_table(dotsTable)

# Rank probes ------------------------------------------------------------------

trank_probes = do.call(rbind, by(tdot, tdot$cell_type, FUN = function(ct) {
	cell_type = ct$cell_type[1]

	ctProbes = bpro
	ctProbes$cell_type = cell_type

	available_probes = unique(ct$probe_label)
	ctProbes = ctProbes[ctProbes$name %in% available_probes,]

	ctProbes = do.call(rbind, by(ctProbes, ctProbes$name, FUN = function(pt) {
		probeDots = ct[ct$probe_label == pt$name[1],]
		pt$dln_min = min(probeDots$lamin_dist_norm, na.rm = T)
		pt$dln_mean = mean(probeDots$lamin_dist_norm, na.rm = T)
		pt$dln_median = median(probeDots$lamin_dist_norm, na.rm = T)
		pt$dln_max = max(probeDots$lamin_dist_norm, na.rm = T)
		return(pt)
	}))

	return(ctProbes)
}))

write.table(trank_probes, file.path(outputFolder,
	sprintf("ranks.probe%s.tsv", suffix)),
	quote = F, sep = "\t", row.names = F)

l = lapply(binSize, FUN = function(size) {
	hsize = ceiling(size/2)
	out = trank_probes
	mids = ceiling((trank_probes$end-trank_probes$start)/2+trank_probes$start)
	out$start = mids-hsize
	out$end = mids+hsize
	write.table(out, file.path(outputFolder,
		sprintf("ranks.%dprobe%s.tsv", size, suffix)),
		quote = F, sep = "\t", row.names = F)
})

# Rank chromosomes -------------------------------------------------------------

trank_chrom = do.call(rbind, by(tdot, tdot$cell_type, FUN = function(ct) {
	cell_type = ct$cell_type[1]

	ctChroms = unique(ct$chrom[order(ct$chromID)])

	ctChroms = do.call(rbind, lapply(ctChroms, FUN = function(chrom) {
		chromDots = ct[ct$chrom == chrom,]
		data.frame(
			chrom = chrom,
			chromID = chrom2chromID(chrom),
			dln_min = min(chromDots$lamin_dist_norm, na.rm = T),
			dln_mean = mean(chromDots$lamin_dist_norm, na.rm = T),
			dln_median = median(chromDots$lamin_dist_norm, na.rm = T),
			dln_max = max(chromDots$lamin_dist_norm, na.rm = T),
			stringsAsFactors = F
		)
	}))
	ctChroms$cell_type = cell_type

	return(ctChroms)
}))

write.table(trank_chrom, file.path(outputFolder,
	sprintf("ranks.chrom%s.tsv", suffix)),
	quote = F, sep = "\t", row.names = F)

# Plot distribution ------------------------------------------------------------

pdf(file.path(outputFolder, sprintf("probes.violin%s.pdf", suffix)),
	width = 30, height = 15)
l = by(tdot, tdot$cell_type, FUN = function(ct) {
	p = ggplot(tdot, aes(
				x = reorder(reorder(probe_label, probeID), chromID),
				y = lamin_dist_norm,
				color = reorder(chrom, chromID)
			)
		) + geom_jitter(size = .5, alpha = .3
		) + geom_violin(color = "black", fill = NA, draw_quantiles = c(.25, .5, .75)
		) + theme_bw(
		) + theme(
			axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
		) + facet_wrap(~reorder(chrom, chromID), scales = "free"
		) + ylim(0, 1+1e-12
		) + ylab("Median normalized lamina distance"
		) + xlab("Probe label"
		) + ggtitle(ct$cell_type[1]
		)
	suppressMessages(print(p))
})
graphics.off()

# END ==========================================================================

################################################################################