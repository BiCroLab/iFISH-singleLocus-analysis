#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: ...
# Email: ...
# Version: X.Y.Z
# Description: ...
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('...', name = 'script.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'mandatory',
	help = 'Mandatory positional argument.')

# Define elective arguments
parser = add_argument(parser, arg = '--elective', short = '-e', type = class(0),
	help = 'Elective argument.', default = 1, nargs = 1)

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

# END --------------------------------------------------------------------------

################################################################################
