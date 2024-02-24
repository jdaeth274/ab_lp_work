###############################################################################
## Manipulating the AMR results from AcBa and lp ##############################
###############################################################################

require(dplyr)


################################################################################
## load up the acba results ####################################################
################################################################################

acba_res <- read.table("~/Dropbox/phd/acba_legion_paper/amrfinder_res/acba_amrfinder.tsv",
                     stringsAsFactors = FALSE, sep = "\t", 
                      header = TRUE)

acba_colnames <- colnames(acba_res)
acba_classes <- sub("^.*\\.\\.","",acba_colnames)
plyr::count(acba_classes)
