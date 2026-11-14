# Run from any working directory: Rscript /path/to/demo.R
# Or in R/RStudio: source("/path/to/demo.R", chdir=TRUE)
args <- commandArgs(trailingOnly=FALSE)
file_arg <- grep("^--file=",args,value=TRUE)
base_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=","",file_arg[1]))) else getwd()
source(file.path(base_dir,"dgp.R"))
source(file.path(base_dir,"vq_test.R"))

# install.packages("cointReg")  # only if missing
# install.packages("fGarch")    # only for type="garch"

# Same errors and regressors for H0/H1, so the change is easy to compare.
for (delta in c(0, 0.9)) {
  set.seed(2026)
  dat <- dgp_main(N=250, delta=delta, kstar=0.5, type="normal")
  cat("\ndelta =",delta,"; coefficient change after index",dat$change_after,"\n")
  print(vq_test(dat$y,dat$x))
}

# Optional paper-definition Q (not identical to the original simulation code):
cat("\nSame H1 data, matrix-root Q:\n")
print(vq_test(dat$y,dat$x,q_method="matrix"))
# Other main error settings: type="ar", "garch", or "smooth".
