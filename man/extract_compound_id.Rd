% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/result_summary.R
\name{extract_compound_id}
\alias{extract_compound_id}
\title{Extract Compound ID from Signature Matching Results}
\usage{
extract_compound_id(
  compound_strings,
  method = "split_colon",
  regex_pattern = NULL,
  part_index = 2
)
}
\arguments{
\item{compound_strings}{Vector of compound identifier strings}

\item{method}{Method for extraction: "split_colon" (default), "split_underscore", or "regex"}

\item{regex_pattern}{Regular expression pattern for extraction (used when method="regex")}

\item{part_index}{Which part to extract when splitting (default: 2)}
}
\value{
Vector of extracted compound identifiers
}
\description{
Helper function to extract compound identifiers from complex compound strings
in signature matching results. This function can be customized based on the
specific format of compound identifiers in your data.
}
\examples{
\dontrun{
# Example compound strings
compounds <- c("CVD001_HEPG2_6H:BRD-K03652504-001-01-9:10.0497",
               "CVD001_HEPG2_6H:BRD-A37828317-001-03-0:10")

# Extract using colon splitting (default)
ids <- extract_compound_id(compounds)

# Extract using custom regex
ids <- extract_compound_id(compounds, method = "regex",
                          regex_pattern = "BRD-[A-Z0-9-]+")
}

}
