% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reactome_data_documentation.R
\docType{data}
\name{reactome_reactions}
\alias{reactome_reactions}
\title{Metadata about Reactome reactions}
\format{
A \code{tibble} with the following fields:
\describe{
  \item{reaction_id}{reaction ID from Reactome}
  \item{reaction_name}{human-readable reaction name}
  \item{compartment_id}{compartment ID from Reactome}
  \item{compartment_name}{human-readable compartment name}
  \item{query_id}{list of 1 or more query IDs corresponding to pathway diagrams containing the reaction}
}
}
\source{
\url{https://reactome.org}
}
\usage{
reactome_reactions
}
\description{
A tibble mapping from reaction IDs to reaction names, compartment info and IDs used
to query for pathway diagrams containing the reaction from the Reactome Diagram Exporter
(\url{https://reactome.org/dev/content-service/diagram-exporter}).
}
\seealso{
\code{\link{reactome_mapping_db}}, \code{\link{reactome_data_version}}
}
\keyword{datasets}
