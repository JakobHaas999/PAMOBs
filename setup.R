##################################
## setup file for this project ###

# all required packages
packages <- c(
  "pammtools",
  "partykit"
)

# Function that checks if a package is installed and installs it if not
# Input:
#  - a character(1) which is the name of the package
#
check_installed <- function(pkg) {
  stopifnot(is.character(pkg) && length(pkg) == 1)

  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("'%s' was not installed", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

for (p in packages) {
  check_installed(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

rm(
  p,
  packages,
  check_installed
)
