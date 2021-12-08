library(data.table)
##source("variance_component.R")
source("vc_test.R")
Sg <- as.matrix(fread("../data/sg.txt.gz"))
Ce <- as.matrix(fread("../data/re.txt.gz"))
                                        # pick only one row of data to test
df <- fread("../data/input.txt.gz", )[,-1]
n <- nrow(Sg)
sqrt_Sg_ginv <- sqrt_ginv(Sg)

system.time({res <- mclapply(1:nrow(df), function(i) {
    unlist(tau_optim(unlist(df[i,]), sqrt_Sg_ginv, Ce))
}, mc.cores = detectCores() - 2)})

res <- do.call("rbind", res)

