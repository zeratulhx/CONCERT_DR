.onLoad <- function(libname, pkgname) {
  # Check if necessary packages are available
  required_packages <- c("BiocManager", "data.table")
  missing_packages <- setdiff(required_packages, installed.packages()[,"Package"])
  
  if (length(missing_packages) > 0) {
    warning("DRnew depends on the following packages that are not installed: ", 
            paste(missing_packages, collapse = ", "), 
            ". Please install them using install.packages() or BiocManager::install().")
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("DRnew: Drug Response Data Analysis Tools")
  packageStartupMessage("Version: ", utils::packageVersion("DRnew"))
}