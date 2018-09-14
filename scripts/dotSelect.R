#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: ...
# Email: ...
# Description: ...
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))

# INPUT ========================================================================

# Create arguent parser
scriptName = "dotSelect.R"
parser = arg_parser('...', name = scriptName)

# Define mandatory arguments
parser = add_argument(parser, arg = 'mandatory',
	help = 'Mandatory positional argument.')

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
attach(p['' != names(p)], warn.conflicts = F)

# RUN ==========================================================================

# END --------------------------------------------------------------------------

################################################################################
