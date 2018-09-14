#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.0.1
# Date: 20180914
# Project: iFISH single locus
# Description: remove discarded channels (after merge) and produce statistics.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd
import sys

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
	Remove discarded channels (after merge) and produce statistics before/after
	filter.
''', formatter_class = argparse.RawDescriptionHelpFormatter)

# Add mandatory arguments
parser.add_argument('inputFolder', type = str,
	help = 'Path to folder with gpseq_fromfish_merge output.')

# Version flag
version = "0.0.1"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

# FUNCTIONS ====================================================================

# RUN ==========================================================================

assert_msg = "input folder not found: %s" % args.inputFolder
assert os.path.isdir(args.inputFolder), assert_msg

tdot_path = "%s/dots.merged.tsv" % args.inputFolder
assert_msg = "merged dot table not found: %s" % tdot_path
assert os.path.isfile(tdot_path), assert_msg

tdot = pd.read_csv(tdot_path, sep = "\t")

ndots_before = tdot.shape[0]

nan = 'nan'.encode()
to_remove = np.where(nan == tdot['probe_label'].values.astype('S'))[0]

ndots_to_remove = len(to_remove)
tdot = tdot.iloc[np.where(nan != tdot['probe_label'].values.astype('S'))[0],:]

ndots_after = tdot.shape[0]

print("""
Dots from merged table: %d
Dots from discarded channels: %d (%.2f%%)
Dots after discarded channels removal: %d (%.2f%%)
""" % (ndots_before,
	ndots_to_remove, ndots_to_remove/ndots_before*100,
	ndots_after, ndots_after/ndots_before*100
))

tdot.to_csv("%s/dots.merged.noDiscardedChannels.tsv" % args.inputFolder,
	sep = "\t", index = False)

# END ==========================================================================

################################################################################
