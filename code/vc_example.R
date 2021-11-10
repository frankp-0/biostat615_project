library(data.table)
source("variance_component.R")
Sg <- as.matrix(fread("../data/sg.txt.gz"))
Ce <- as.matrix(fread("../data/re.txt.gz"))
# pick only one row of data to test
df <- fread("../data/input.txt.gz", )[80000,-1]
n <- nrow(Sg)
ind <- seq(1, 2 * n, 2)
sqrt_Sg_ginv = sqrt_ginv(Sg)
run_vc_optimizer(df, sqrt_Sg_ginv, Ce, ind)


