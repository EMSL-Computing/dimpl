# Data object documentation for Reactome data 


#' Mapping from proteins to reactions and metabolites
#' 
#' A tibble that maps from proteins to reactions and the metabolites consumed and produced
#' from those reactions. The proteins are Uniprot IDs extracted from the reaction 
#' \code{modifier} field in the Reactome SBML source files. Metabolites are CHEBI IDs extracted
#' from the reaction \code{consumed} and \code{produced} fields of the SBML files. 
#' 
#' @format A \code{tibble} with the following fields:
#' \describe{
#'   \item{protein}{Uniprot ID (character string)}
#'   \item{reaction_id}{reaction ID from Reactome}
#'   \item{consumed}{list of double-valued CHEBI IDs of metabolites consumed in the corresponding reaction}
#'   \item{produced}{list of double-valued CHEBI IDs of metabolites produced in the corresponding reaction}
#' }
#' @source \url{https://reactome.org}
#' @seealso \code{\link{reactome_reactions}}, \code{\link{reactome_data_version}}
"reactome_mapping_db"



#' Metadata about Reactome reactions
#' 
#' A tibble mapping from reaction IDs to reaction names, compartment info and IDs used
#' to query for pathway diagrams containing the reaction from the Reactome Diagram Exporter
#' (\url{https://reactome.org/dev/content-service/diagram-exporter}).
#' 
#' @format A \code{tibble} with the following fields:
#' \describe{
#'   \item{reaction_id}{reaction ID from Reactome}
#'   \item{reaction_name}{human-readable reaction name}
#'   \item{compartment_id}{compartment ID from Reactome}
#'   \item{compartment_name}{human-readable compartment name}
#'   \item{query_id}{list of 1 or more query IDs corresponding to pathway diagrams containing the reaction}
#' }
#' @source \url{https://reactome.org}
#' @seealso \code{\link{reactome_mapping_db}}, \code{\link{reactome_data_version}}
"reactome_reactions"



#' Reactome version used to produce data objects for mapping
#' 
#' A string containing the version of the Reactome data that was parsed to form the
#' \code{reactome_mapping_db} and \code{reactome_reactions} objects.
#' 
#' @source \url{https://reactome.org}
#' @seealso \code{\link{reactome_mapping_db}}, \code{\link{reactome_reactions}}
"reactome_data_version"