#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Description: select N brightest dots.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(parallel))
suppressMessages(library(readr))

# INPUT ========================================================================

# Create arguent parser
scriptName = "dotSelect.R"
parser = arg_parser('Select the top brightest dots, per cell-type.',
	name = scriptName)

# Define mandatory arguments
parser = add_argument(parser, arg = 'dotsTable', type = class(""),
	help = 'Path to dots table.')
parser = add_argument(parser, arg = '--ndots', short = '-N', type = class(""),
	help = 'Space-separated expected_Ndots:cell_type couples.', nargs = Inf)

# Define elective arguments
parser = add_argument(parser, arg = '--value-label', short = '-l',
	help = 'Dot intensity column label. Default: "Value"', default = "Value",
	nargs = 1, type = class(""),)
parser = add_argument(parser, arg = '--threads', short = '-t', type = class(0),
	help = 'Number of threads for parallelization.', default = 1, nargs = 1)

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

# FUNCTIONS ====================================================================

parseNdots = function(x) {
	x = unlist(strsplit(x, ':', fixed = T))
	outl = list(as.numeric(x[1]))
	names(outl) = x[2]
	return(outl)
}

selectBrightestDots = function(ct, expected_ndots) {
	cell_type = ct$cell_type[1]
	expN = expected_ndots[[cell_type]]

	if ( !cell_type %in% names(expected_ndots) ) {
		cat(paste0(sprintf("Skipped '%s' cell type,", cell_type),
			" no expected number of dots provided.\n"))
		return(NULL)
	}

	ct$ncuID = paste0(ct$nuID, '~', ct$Channel)
	ncuID = unique(ct$ncuID)

	# First identify good sets (channel-cells): those with ndots <= expN
	# And the bad sets, that need further action
	dots_count = table(ct$ncuID)
	good_sets = names(dots_count)[dots_count <= expN]
	bad_sets = ncuID[!ncuID %in% good_sets]
	
	# Select dots in bad sets
	ct2 = do.call(rbind, mclapply(bad_sets, FUN = function(i) {
		st = ct[ct$ncuID == i,]
		
		if ( expN >= nrow(st) ) stop("unexpected number of dots found.")
		
		st = st[order(st[, value_label], decreasing = T),]
		
		return(st[1:expN,])
	}, mc.cores = threads))

	# Add god sets
	ct2 = rbind(ct2, ct[ct$ncuID %in% good_sets,])

	cat(sprintf("%d/%d (%.2f%%) dots selected as brightest in %s.\n",
		nrow(ct2), nrow(ct), nrow(ct2)/nrow(ct)*100, cell_type))

	ct2 = ct2[, which(colnames(ct2) != "ncuID")]
	return(ct2)
}

# RUN ==========================================================================

expected_ndots = do.call(c, lapply(ndots, FUN = parseNdots))

cat(sprintf("Reading input dot table: '%s'\n", dotsTable))
tdot = suppressMessages(as.data.frame(
	read_delim(dotsTable, "\t", progress = F)))

cat(sprintf("Selecting dots...\n"))
tdot2 = do.call(rbind, by(tdot, tdot$cell_type,
	FUN = selectBrightestDots, expected_ndots))

outTable = file.path(dirname(dotsTable), sprintf("%s.%s.tsv",
	tools::file_path_sans_ext(basename(dotsTable)), "brightest"))
cat(sprintf("Writing output: '%s'\n", outTable))
write.table(tdot2, outTable, quote = F, row.names = F, sep = "\t")

# END --------------------------------------------------------------------------

################################################################################
