% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmap_utils.R
\name{interactive_cmap_setup}
\alias{interactive_cmap_setup}
\title{Interactive CMap configuration and workflow}
\usage{
interactive_cmap_setup(
  siginfo_file = NULL,
  geneinfo_file = NULL,
  gctx_file = NULL,
  config_dir = "conf",
  output_dir = "output",
  verbose = TRUE
)
}
\arguments{
\item{siginfo_file}{Path to siginfo_beta.txt file (default: NULL to prompt user)}

\item{geneinfo_file}{Path to geneinfo_beta.txt file (default: NULL to prompt user)}

\item{gctx_file}{Path to GCTX file (default: NULL to prompt user)}

\item{config_dir}{Directory to save configuration files (default: "conf")}

\item{output_dir}{Directory to save output files (default: "output")}

\item{verbose}{Logical; whether to print progress messages (default: TRUE)}
}
\value{
A list with workflow results
}
\description{
Provides a user-friendly interface to configure and run the CMap workflow.
Tries to use the GUI interface if available, falls back to console version if not.
}
\examples{
\dontrun{
results <- interactive_cmap_setup()
}

}
