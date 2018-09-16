#
# 20180912 - Gabriele Girelli
# Project: GPSeq manuscript figures
# 
# Aim:
#   - Collect common functions in one place.
#



# DEPENDENCIES =================================================================

suppressMessages(library(readr))

# GENERAL ======================================================================

chrom2chromID = function(chrom) {
    if ( grepl(":", chrom) ) {
        return(floor(as.numeric(gsub(":", ".", substr(chrom, 4, nchar(chrom))))))
    } else {
        chromID = substr(chrom, 4, nchar(chrom))
        if ( "X" == chromID ) chromID = 23
        if ( "Y" == chromID ) chromID = 24
        return(as.numeric(chromID))
    }
}

add_chrom_ID = function(data) {
    cid_table = data.frame(chrom = unique(data$chrom), stringsAsFactors = F)
    cid_table$chromID = unlist(lapply(cid_table$chrom, FUN = chrom2chromID))
    data$chromID = cid_table$chromID[match(data$chrom, cid_table$chrom)]
    return(data)
}

read_probe_bed = function(path) {
    # Requires BED4 file with no header and with track line
    bpro = read.delim(path, as.is = T, header = F, skip = 1)
    colnames(bpro) = c("chrom", "start", "end", "name")
    bpro = add_chrom_ID(bpro)
    return(bpro)
}

read_dot_table = function(path) {
    tdot = suppressMessages(as.data.frame(read_delim(path, "\t")))
    probe2chrom = do.call(rbind, lapply(unique(tdot$probe_label),
        FUN = function(label) {
            probeBits = unlist(strsplit(label, '.', fixed = T))
            data.frame(
                label = label,
                chrom = probeBits[1],
                probeID = as.numeric(probeBits[2]),
                stringsAsFactors = F
            )
        })
    )
    tdot$chrom = probe2chrom$chrom[
        match(tdot$probe_label, probe2chrom$label)]
    tdot$probeID = probe2chrom$probeID[
        match(tdot$probe_label, probe2chrom$label)]
    tdot = add_chrom_ID(tdot)
    return(tdot)
}

# END ==========================================================================

################################################################################