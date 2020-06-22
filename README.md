# irg
R Package for Granger-Causal Analysis of Irregularly Sampled Signals

For users who are interested in having the latest developments, the GitHub version is ideal although more dependencies are required to run a stable version of the package. Most importantly, users must have a (C++) compiler installed on their machine that is compatible with R (e.g. Clang).

# Install dependencies
install.packages(c("RcppArmadillo","devtools","knitr","rmarkdown"))

# Install the package from GitHub without Vignettes/User Guides
devtools::install_github("SMAC-Group/irg")

# Install the package with Vignettes/User Guides 
devtools::install_github("SMAC-Group/irg", build_vignettes = TRUE)

The setup to obtain the development version of simts is platform dependent.
