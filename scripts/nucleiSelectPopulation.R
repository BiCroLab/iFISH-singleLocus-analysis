#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Description:
# 	- Plot distribution of intensity interal over volume and flattened area.
# 	- Fit sum of Gaussians and select k*sigma range around the major Peak.
# 	- Plot scatterplot of the two and color the selected ones differently.
# 	- Subselect nuclei and filter dot table accordingly.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(ggplot2))
suppressMessages(library(GGally))
suppressMessages(library(parallel))
suppressMessages(library(readr))
suppressMessages(library(viridis))

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
parser = add_argument(parser, arg = '--ksigma', short = '-k', type = class(0),
	help = 'Sigma constant for interval definition.', default = 3, nargs = 1)
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
if ( !file.exists(nucleiTable) )
	stop(sprintf("nuclei table not found: %s", nucleiTable))
if ( !file.exists(compTable) )
	stop(sprintf("compartments table not found: %s", compTable))

if ( !dir.exists(outputFolder) )
	stop(sprintf("output folder not found: %s", outputFolder))

# FUNCTION =====================================================================

cross_x = function(xs, ys) {
	idx = which(unlist(lapply(1:(length(ys)-1), FUN = function(i) {
		ys[i] * ys[i+1] < 0
	})))

	P1 = data.frame(x = xs[idx], y = ys[idx])
	P2 = data.frame(x = xs[idx+1], y = ys[idx+1])

	m = (P2$y - P1$y) / (P2$x - P1$x)
	i = P2$y - m * P2$x

	types = unlist(lapply(P1$y, FUN = function(x) {
		if( x >= 0 ) {
			return("maxima")
		} else {
			return("minima")
		}
	}))

	return(data.frame(i = idx + 1, x = -i/m, type = types))
}

fwhm = function(m, xs, ys) {
	hm = ys[xs == m] / 2

	crossings = cross_x(xs, ys - hm)

	fir = crossings$x[crossings$x <= m]
	fir = fir[length(fir)]
	sec = crossings$x[crossings$x >= m][1]

	return(2 * min(m - fir, sec - m))
}

calc_feature_range = function(data, column, k_sigma = 3, plot = F) {

	# Calculate starting points for fitting

	ddens = density(data[, column])
	ddens = data.frame(x = ddens$x, density = ddens$y, stringsAsFactors = F)

	peaks = cross_x(ddens$x[-1], diff(ddens$density))
	maxima = peaks[peaks$type == "maxima",]
	maxima = maxima[order(ddens$density[maxima$i], decreasing = T),]
	peak1 = maxima[1,]
	peak2 = maxima[maxima$x > peak1$x,][1,]

	f1 = fwhm(ddens$x[peak1$i], ddens$x, ddens$density)
	sigma1 = f1 / sqrt(8 * log(2))
	f2 = fwhm(ddens$x[peak2$i], ddens$x, ddens$density)
	sigma2 = f2 / sqrt(8 * log(2))

	C1 = ddens$density[peak1$i]
	C2 = ddens$density[peak2$i]
	mean1 = peak1$x
	mean2 = peak2$x

	# Fit Sum of Gaussians
	
	model = tryCatch({
		nls(ddens$density ~ (C1 * exp(-(ddens$x-mean1)**2/(2 * sigma1**2)) +
			C2 * exp(-(ddens$x-mean2)**2/(2 * sigma2**2))), data=ddens,
			start=list(C1=C1, mean1=mean1, sigma1=sigma1,
			         C2=C2, mean2=mean2, sigma2=sigma2), algorithm="port") 
	}, error = function(e) {
		nls(ddens$density ~ (C1 * exp(-(ddens$x-mean1)**2/(2 * sigma1**2))),
			data=ddens, start=list(C1=C1, mean1=mean1, sigma1=sigma1),
			algorithm="port") 
	})
	cfs = as.data.frame(t(coef(model)), stringsAsFactors = T)
	ddens$G1 = cfs$C1*exp(-(ddens$x-cfs$mean1)**2/(2*cfs$sigma1**2))
	

	if ("C2" %in% names(cfs)) {
		# Enforce 2nd Gaussian to be AFTER and LOWER than the 1st
		# by reverting to single Gaussian fitting
		if (mean2 < (mean1 + 2 * sigma1) || C2 > C1) {
			model = nls(
				ddens$density ~ (C1 * exp(-(ddens$x-mean1)**2/(2 * sigma1**2))),
				data=ddens, start=list(C1=C1, mean1=mean1, sigma1=sigma1),
				algorithm="port") 
			cfs = as.data.frame(t(coef(model)), stringsAsFactors = T)
		}
	}
	if ("C2" %in% names(cfs))
		ddens$G2 = cfs$C2*exp(-(ddens$x-cfs$mean2)**2/(2*cfs$sigma2**2))

	return(list(
		range = c(mean1-k_sigma*sigma1, mean1+k_sigma*sigma1),
		density = ddens
	))
}

select_dataset_population = function(dataset, tnuc) {
	dt = tnuc[tnuc$dataset == dataset,]

	sumI = calc_feature_range(dt, "sumI", ksigma)
	area = calc_feature_range(dt, "flat_size", ksigma)

	dt$selected = 1
	dt$selected[dt$sumI < sumI$range[1] | dt$sumI > sumI$range[2]] = 0
	dt$selected[dt$flat_size < area$range[1] | dt$flat_size > area$range[2]] = 0
	n_selected = sum(dt$selected)
	dt$selected = factor(dt$selected)

	color = list(
		range = "#b2182b",
		G1 = "#01665e",
		G2 = "#8c510a",
		selected = "#2166ac",
		discarded = "#d6604d"
	)

	p1 = ggplot(dt, aes(x = flat_size, y = sumI, color = selected)
		) + geom_point(size = .5, alpha = .5
		) + theme_bw(
		) + scale_color_manual(values = c(color$selected, color$discarded)
		) + geom_vline(xintercept = area$range,
			color = color$range, linetype = "dashed"
		) + geom_hline(yintercept = sumI$range,
			color = color$range, linetype = "dashed"
		) + xlim(
				min(c(area$density$x, area$range[1])),
				max(c(area$density$x, area$range[2])
			)
		) + ylim(
				min(c(sumI$density$x, sumI$range[1])),
				max(c(sumI$density$x, sumI$range[2])
			)
		)

	p2 = ggplot(sumI$density, aes(x = x, y = density)) + geom_line(
		) + coord_flip() + scale_y_reverse(
		) + theme_bw(
		) + geom_vline(xintercept = sumI$range,
			color = color$range, linetype = "dashed"
		) + geom_line(aes(y = G1) , color = color$G1,
			alpha = 0.75, size = .5
		)
	if ( "G2" %in% colnames(sumI$density) )
		p2 = p2 + geom_line(aes(y = G2), color = color$G2,
			alpha = 0.75, size = .5)

	p3 = ggplot(area$density, aes(x = x, y = density)) + geom_line(
		) + scale_y_reverse(
		) + theme_bw(
		) + geom_vline(xintercept = area$range,
			color = color$range, linetype = "dashed"
		) + geom_line(aes(y = G1) , color = color$G1,
			alpha = 0.75, size = .5
		)
	if ( "G2" %in% colnames(area$density) )
		p3 = p3 + geom_line(aes(y = G2), color = color$G2,
			alpha = 0.75, size = .5)

	p = ggmatrix(list(p2, p1, NULL, p3), 2, 2,
		xProportions = c(.2, .8), yProportions = c(.8, .2),
		xlab = "Flattened area [px]", ylab = "Intensity integral [a.u.]",
		title = paste0(
			sprintf("%s - interval is %d·sigma around G1 peak\n",
				dt$dataset[1], ksigma),
			sprintf("Selected %d/%d (%.2f%%) nuclei (sumI & flat_area)\n",
				n_selected, nrow(dt), n_selected/nrow(dt)*100),
			sprintf("Flattened area interval: [%.2E, %.2E] [vx]\n",
				area$range[1], area$range[2]),
			sprintf("SumI interval: [%.2E, %.2E] [a.u.]",
				sumI$range[1], sumI$range[2])
		)
	)

	data = data.frame(
		dataset = dt$dataset[1],
		area_min = area$range[1],
		area_max = area$range[2],
		sumI_min = sumI$range[1],
		sumI_max = sumI$range[2],
		nselected = n_selected,
		nnuclei = nrow(dt),
		cell_type = dt$cell_type[1],
		stringsAsFactors = F
	)

	return(list(
		nuclei_table = dt,
		data = data,
		plot = p
	))
}

# RUN ==========================================================================

tdot = suppressMessages(as.data.frame(
	read_delim(dotsTable, "\t", progress = F)))
tnuc = suppressMessages(as.data.frame(
	read_delim(nucleiTable, "\t", progress = F)))
tcom = suppressMessages(as.data.frame(
	read_delim(compTable, "\t", progress = F)))

# Stats on dots
n_outer_dots = sum(is.na(tdot$cell_ID))
cat(sprintf("%d/%d (%.2f%%) dots outside of nuclei.\n",
	n_outer_dots, nrow(tdot), n_outer_dots/nrow(tdot)*100))
tdot = tdot[!is.na(tdot$cell_ID),]
cat("Removed dots outside of nuclei.\n")

# Add nuclear IDs
tdot$nuID = paste0(tdot$dataset, '~', tdot$date, '~', tdot$session, '~',
	tdot$File, '~', tdot$cell_ID)
tdot$nuID[is.na(tdot$cell_ID)] = NA
tcom$nuID = paste0(tcom$dataset, '~', tcom$date, '~', tcom$session, '~',
	tcom$File, '~', tcom$cell_ID)
tnuc$nuID = paste0(tnuc$dataset, '~', tnuc$date, '~', tnuc$session, '~',
	tnuc$s, '~', tnuc$n)

# Stats on nuclei counts
nuclei_wdots = unique(tdot$nuID[!is.na(tdot$nuID)])
n_nuclei_wdots = length(nuclei_wdots)
cat(sprintf("Found %d (%d) nuclei.\n", nrow(tnuc), nrow(tcom)))
cat(sprintf("Found %d (%.2f%%) nuclei with dots.\n",
	n_nuclei_wdots, n_nuclei_wdots/nrow(tnuc)*100))
cat("Summary of #nuclei per dataset:\n")
summary(unlist(by(tnuc, tnuc$dataset, nrow)))

# Discard nuclei without dots
tnuc = tnuc[tnuc$nuID %in% nuclei_wdots,]
tcom = tcom[tcom$nuID %in% nuclei_wdots,]
cat("Removed nuclei w/o dots.\n")
cat(sprintf("Reduced to %d (%d) nuclei.\n", nrow(tnuc), nrow(tcom)))

cat("Summary of #nuclei per dataset:\n")
summary(unlist(by(tnuc, tnuc$dataset, nrow)))

# Filter per dataset -----------------------------------------------------------

tmp = mclapply(unique(tnuc$dataset),
	FUN = select_dataset_population, tnuc, mc.cores = threads)

# Export dataset recap statistics
dataset_stats = do.call(rbind, lapply(tmp, FUN = function(x) x$data))
write.table(dataset_stats,
	file.path(outputFolder, "nuclei_select_population.datasets.tsv"),
	row.names = F, quote = F, sep = "\t")

# Export plot
pdf(file.path(outputFolder, "nuclei_select_population.selection.pdf"),
	width = 8, height = 9)
p = lapply(tmp, FUN = function(x) print(x$plot))
graphics.off()

pdf(file.path(outputFolder, "nuclei_select_population.datasets.pdf"),
	width = 20, height = 9)
ggplot(dataset_stats, aes(x = dataset, y = nselected, fill = dataset)
	) + geom_bar(stat = "identity"
	) + geom_hline(yintercept = 100, linetype = "dashed", color = "black"
	) + theme_bw(
	) + theme(
		axis.text.x = element_text(angle = 90, hjust = 1)
	) + facet_wrap(~cell_type
	) + ggtitle(sprintf("%d datasets below the threshold of 100 nuclei.",
		sum(dataset_stats$nselected <= 100))
	) + guides(fill = F
	)
graphics.off()

# Assemble selection
tnuc2 = do.call(rbind, lapply(tmp, FUN = function(x) x$nuclei_table))

n_selected = sum(dataset_stats$nselected)
cat(sprintf("%d/%d (%.2f%%) nuclei selected.\n",
	n_selected, nrow(tnuc2), n_selected/nrow(tnuc2)*100))

tnuc2 = tnuc2[as.numeric(as.character(tnuc2$selected)) == 1,]
tcom2 = tcom[tcom$nuID %in% tnuc2$nuID,]
tdot2 = tdot[tdot$nuID %in% tnuc2$nuID,]

cat(sprintf("%d (%d) nuclei after selection.\n", nrow(tnuc2), nrow(tcom2)))
cat("Summary of #nuclei per dataset:\n")
summary(unlist(by(tnuc2, tnuc2$dataset, nrow)))
cat(sprintf("%d/%d (%.2f%%) dots left after nuclear selection.\n",
	nrow(tdot2), nrow(tdot), nrow(tdot2)/nrow(tdot)*100))

write.table(tdot2, file.path(dirname(dotsTable), sprintf("%s%s.tsv",
	tools::file_path_sans_ext(basename(dotsTable)), ".G1selected")),
	quote = F, row.names = F, sep = "\t")
write.table(tnuc2, file.path(dirname(nucleiTable), sprintf("%s%s.tsv",
	tools::file_path_sans_ext(basename(nucleiTable)), ".G1selected")),
	quote = F, row.names = F, sep = "\t")
write.table(tcom2, file.path(dirname(compTable), sprintf("%s%s.tsv",
	tools::file_path_sans_ext(basename(compTable)), ".G1selected")),
	quote = F, row.names = F, sep = "\t")

# END --------------------------------------------------------------------------

################################################################################
