#' Find matches in regulation of proteins and metabolites against a database
#' 
#' [TODO: Insert more information about the relevance of this]
#'
#' @param db Database against which to compare protein and metabolite regulation data
#' @param prot_reg Protein regulation vector with names. Values should be < 0 for downregulated elements and > 0 for upregulated elements.
#' @param metab_reg Metabolite regulation vector. Values should be < 0 for downregulated elements and > 0 for upregulated elements.
#'
#' @return A list of matching results, each list element corresponds to a single protein. Each list item contains elements:
#'   \item{consumed}{consumed metabolites in \code{db} that are downregulated in \code{metab_reg}}
#'   \item{produced}{produced metabolites in \code{db} that are upregulated in \code{metab_reg}}
#'   \item{nondir_down}{non-directional metabolites in \code{db} that are downregulated in \code{metab_reg}}
#'   \item{nondir_up}{non-directional metabolites in \code{db} that are upregulated in \code{metab_reg}}
#' Note that non-directional metabolites are those that are listed as being both consumed and produced in the database.
#' @importFrom dplyr %>%
#' @export
find_matches <- function(db, prot_reg, metab_reg) {

  # get the portions of the database corresponding to proteins positively or negatively regulated
  prot_reg <- prot_reg[prot_reg != 0]
  db <- db %>% 
    dplyr::filter(protein %in% names(prot_reg))
  db <- db %>% 
    dplyr::left_join(data.frame(protein=names(prot_reg), protein_regulation=prot_reg))

  metab_reg <- metab_reg[metab_reg != 0]
  
  intersect.null <- function(x,y) {
    result <- intersect(x,y)
    if (length(result) == 0) result <- NULL
    return(result)
  }
  fn <- function(consumed, produced, prot_reg, metab_reg) {
    nondir  <- as.character(intersect.null(consumed, produced))
    consumed <- setdiff(as.character(consumed), nondir)
    produced <- setdiff(as.character(produced), nondir)
    
    metab_reg <- metab_reg * prot_reg
    cons <- intersect.null(consumed, names(metab_reg[metab_reg < 0]))
    prod <- intersect.null(produced, names(metab_reg[metab_reg > 0]))
    ndd <- intersect.null(nondir, names(metab_reg[metab_reg < 0]))
    ndu <- intersect.null(nondir, names(metab_reg[metab_reg > 0]))
    
    if (length(c(cons, prod, ndd, ndu)) > 0) {
      return(list(consumed=cons, produced=prod, nondir_down=ndd, nondir_up=ndu))
    } else {
      return(NULL)
    }
    
  }
  
  matches <- lapply(1:nrow(db), function(i) fn(db[i, "consumed"][[1]], db[i, "produced"][[1]], 
                                                db[i, "protein_regulation"][[1]], metab_reg))
  result <- db %>% dplyr::ungroup()
  result$matches <- matches
  result <- result %>% dplyr::filter(unlist(purrr::map(matches, length)) > 0)

  return(result)
}
