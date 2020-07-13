# This script parses the Reactome sbml files into R data objects which are 
# used in mapping between proteins and metabolites. Resulting data objects
# will be stored in the dimpl/data directory using the usethis::use_data
# function. 
# 
# To use this file to update package data: 
# 1) Download archive of updated Reactome SBML files (https://reactome.org/download-data)
#    See Specialized Data Formats section --> All Species reactions in SBML
# 2) Unzip archive
# 3) Specify parameters in the top section of the script


###############################################################################
## PARAMETERS

# Root of dimpl package directory
dimpl_pkg_dir <- "D:/Files/iPMART/R/dimpl"

# Directory of Reactome SBML files
data_dir <- "D:/Files/iPMART/Data/all_species.3.1.sbml"

# Reactome data version (from downloaded data filename)
reactome_data_version <- "3.1"

# If TRUE intermediate data steps will be saved as RDS files in data_dir
save_intermediate_data <- TRUE

# Number of cores to use for parallel processing
num_cores <- 10

## END PARAMETERS

wd <- getwd()
setwd(dimpl_pkg_dir)

library(xml2)
library(dplyr)
library(doParallel)
library(usethis)

# Cluster for parallel processing SBML files:
cl <- makeCluster(num_cores)
registerDoParallel(cl)


###############################################################################
## FUNCTIONS FOR EXTRACTING DATA FROM XML 

# extract CHEBI ids
get_chebi_ids<-function(x){
  s<-grep("CHEBI",x)
  if(length(s>0)){
    x<-x[s]
    x<-as.numeric(gsub("^.*CHEBI:","",x))
    return(x)
  } else {
    return(NULL)
  }
}

# extract Uniprot ids
get_uniprot_ids<-function(x){
  s<-grep("uniprot",x)
  if(length(s>0)){
    x<-x[s]
    x<-gsub("^.*uniprot\\/","",x)
    return(x)
  } else {
    return(NULL)
  }
}

###############################################################################
## REACTIONS

get_reaction_data <- function(root_node) {
  rxns <- xml_children(xml_find_all(root_node, ".//listOfReactions"))
  rxn.data <- xml_attrs(rxns)
  rxn.data <- do.call(rbind, rxn.data)
  rxn.data <- matrix(rxn.data[, c("id", "name", "compartment")], ncol=3)  #in case there's only 1 row
  
  colnames(rxn.data) <- c("reaction_id", "reaction_name", "compartment_id")
  rxn.data <- data.frame(rxn.data, stringsAsFactors = FALSE)
  
  query.ids <- lapply(rxns, get_reaction_query_id)
  names(query.ids) <- rxn.data$reaction_id
  
  species.data <- get_species_data(root_node)
  reactants <- lapply(rxns, get_reactant_data, species_df=species.data)
  names(reactants) <- rxn.data$reaction_id
  products <- lapply(rxns, get_product_data, species_df=species.data)
  names(products) <- rxn.data$reaction_id
  modifiers <- lapply(rxns, get_modifier_data, species_df=species.data)
  names(modifiers) <- rxn.data$reaction_id
  
  #return(list(reactions=rxn.data, reactants=reactants, products=products, modifiers=modifiers))
  
  compartments <- get_compartment_data(root_node) 
  
  reactions_tbl <- rxn.data %>%
    left_join(compartments) %>%
    left_join(tibble(reaction_id=names(query.ids), query_id=query.ids)) %>%
    left_join(tibble(reaction_id=names(reactants), reactants=reactants)) %>%
    left_join(tibble(reaction_id=names(products), products=products)) %>%
    left_join(tibble(reaction_id=names(modifiers), modifiers=modifiers)) #%>%
  #     select(-c(compartment_id, reaction_id)) 
  
  reactions_tbl <- reactions_tbl %>%
    filter(unlist(purrr::map(modifiers, function(x) !is.null(x))))
  
  return(reactions_tbl)
}

# Get the ID used for querying for reaction diagram (root_node is reaction node)
get_reaction_query_id <- function(root_node) {
  nodes <- xml_find_all(root_node, ".//annotation/rdf:RDF/rdf:Description/bqbiol:is/rdf:Bag/rdf:li")
  ids <- xml_attr(nodes, "resource")
  ind <- grepl("https://reactome.org/content/detail", ids)
  ids <- unique(gsub("https://reactome.org/content/detail/", "", ids[ind]))
  return(ids)
}

get_reactant_data <- function(root_node, species_df=NULL) {
  # root_node is the root of a single reaction
  rxnts <- xml_children(xml_find_all(root_node, ".//listOfReactants"))
  if (length(rxnts) == 0) return(NULL)
  rxnt.data <- xml_attrs(rxnts)
  rxnt.data <- do.call(rbind, rxnt.data)
  rxnt.data <- matrix(rxnt.data[, c("id", "species")], ncol=2) #in case there's only 1 row
  colnames(rxnt.data) <- c("reactant_id", "species_id")
  rxnt.data <- data.frame(rxnt.data, stringsAsFactors = FALSE)
  
  if (!(all(is.null(species_df)))) {
    rxnt.data <- left_join(rxnt.data, species_df) %>%
      select(-c(species_id, reactant_id))
  }
  
  # extract CHEBI IDs from the is_a column and drop the has_part column
  chebi <- unique(unlist(lapply(rxnt.data$is_a, get_chebi_ids)))
  chebi2 <- unique(unlist(lapply(rxnt.data$part_of, get_chebi_ids)))

  return(unique(c(chebi, chebi2)))
}

get_product_data <- function(root_node, species_df=NULL) {
  # root_node is the root of a single reaction
  prods <- xml_children(xml_find_all(root_node, ".//listOfProducts"))
  if (length(prods) == 0) return(NULL)
  prod.data <- xml_attrs(prods)
  prod.data <- do.call(rbind, prod.data)
  prod.data <- matrix(prod.data[, c("id", "species")], ncol=2) #in case there's only 1 row
  colnames(prod.data) <- c("product_id", "species_id")
  prod.data <- data.frame(prod.data, stringsAsFactors = FALSE)
  
  if (!(all(is.null(species_df)))) {
    prod.data <- left_join(prod.data, species_df) %>%
      select(-c(species_id, product_id))
  }
  
  # extract CHEBI IDs from the is_a column and drop the has_part column
  chebi <- unique(unlist(lapply(prod.data$is_a, get_chebi_ids)))
  chebi2 <- unique(unlist(lapply(prod.data$part_of, get_chebi_ids)))

  return(unique(c(chebi, chebi2)))
}

get_modifier_data <- function(root_node, species_df=NULL) {
  # root_node is the root of a single reaction
  mods <- xml_children(xml_find_all(root_node, ".//listOfModifiers"))
  if (length(mods) == 0) return(NULL)
  mod.data <- xml_attrs(mods)
  mod.data <- do.call(rbind, mod.data)
  mod.data <- matrix(mod.data[, c("id", "species")], ncol=2) #in case there's only 1 row
  colnames(mod.data) <- c("modifier_id", "species_id")
  mod.data <- data.frame(mod.data, stringsAsFactors = FALSE)
  
  if (!(all(is.null(species_df)))) {
    mod.data <- left_join(mod.data, species_df) %>%
      select(-c(species_id, modifier_id))
  }
  
  # extract Uniprot IDs from the is_a and has_part columns 
  is_a_uniprot <- unlist(lapply(mod.data$is_a, get_uniprot_ids))
  has_part_uniprot <- unlist(lapply(mod.data$has_part, get_uniprot_ids))
  uniprot_ids <- unique(c(is_a_uniprot, has_part_uniprot))

  return(uniprot_ids)
}

get_compartment_data <- function(root_node) {
  cmps <- xml_children(xml_find_all(root_node, ".//listOfCompartments"))
  if (length(cmps) == 0) return(NULL)
  comp.data <- xml_attrs(cmps)
  comp.data <- do.call(rbind, comp.data)
  comp.data <- matrix(comp.data[, c("id", "name")], ncol=2) #in case there's only 1 row
  colnames(comp.data) <- c("compartment_id", "compartment_name")
  return(data.frame(comp.data, stringsAsFactors = FALSE))
}


###############################################################################
## SPECIES

get_species_data <- function(root_node) {
  spec <- xml_children(xml_find_all(root_node, ".//listOfSpecies"))
  if (length(spec) == 0) return(NULL)
  species.data <- xml_attrs(spec)
  species.data <- lapply(species.data, function(x) c(x["id"], x["name"]))
  species.data <- do.call(rbind, species.data)
  colnames(species.data) <- c("species_id", "species_name")
  species.data <- data.frame(species.data, stringsAsFactors = FALSE)
  
  parts <- lapply(spec, get_parts_data)
  names(parts) <- species.data$species_id
  is <- lapply(spec, get_is_data)
  names(is) <- species.data$species_id
  partOf <- lapply(spec, get_part_of_data)
  names(partOf) <- species.data$species_id
  
  #return(list(species=species.data, parts=parts, is=is))
  
  species_tbl <- species.data %>%
    left_join(tibble(species_id=names(is), is_a=is)) %>%
    left_join(tibble(species_id=names(parts), has_part=parts)) %>%
    left_join(tibble(species_id=names(partOf), part_of=partOf))
  return(species_tbl)
}

get_parts_data <- function(root_node) {
  parts <- xml_find_all(root_node, ".//bqbiol:hasPart//rdf:li")
  if (length(parts) == 0) return(NULL)
  parts_ids <- xml_attr(parts, "resource")
  return(parts_ids)
}

get_part_of_data <- function(root_node) {
  parts <- xml_find_all(root_node, ".//bqbiol:isPartOf//rdf:li")
  if (length(parts) == 0) return(NULL)
  parts_ids <- xml_attr(parts, "resource")
  if (!is.null(parts_ids)) message("Found an isPartOf element!")
  return(parts_ids)
}

get_is_data <- function(root_node) {
  is <- xml_find_all(root_node, ".//bqbiol:is//rdf:li")
  if (length(is) == 0) return(NULL)
  is_ids <- xml_attr(is, "resource")
  return(is_ids)
}


###############################################################################
## PARSE ALL FILES

fnames <- list.files(data_dir, full.names = TRUE, pattern="*.sbml")
print(sprintf("%d files to process", length(fnames)))

parse_file <- function(fname) {
  require(xml2)
  require(dplyr)
  suppressMessages({
    file.data <- read_xml(fname)
    file.data <- xml_ns_strip(file.data)
    reactions <- get_reaction_data(file.data)
    #species <- get_species_data(file.data)
  })  
  return(reactions)
}


# Parse files in parallel:
print(paste("Starting file parsing:", Sys.time()))
all_reactions <- foreach(i = 1:length(fnames)) %dopar% parse_file(fnames[i])
print(paste("Finished file parsing:", Sys.time()))

if (save_intermediate_data) 
  saveRDS(all_reactions, file=file.path(data_dir, "reactome_data_intermediate1.RDS"))

# Join all results together
all_reactions2 <- do.call(rbind, all_reactions)
if (save_intermediate_data) saveRDS(all_reactions2, file=file.path(data_dir, "reactome_data_intermediate.RDS"))


# Reorganize data to collect proteins ("modifiers") across all reactions
.data <- tidyr::unnest(all_reactions2, modifiers) 
.data <- .data %>%  group_by(modifiers)
.data <- .data %>% 
  filter(unlist(purrr::map(reactants, function(x) !is.null(x))) | unlist(purrr::map(products, function(x) !is.null(x))))

reactome_reactions <- .data %>% 
  ungroup() %>% 
  select(reaction_id, reaction_name, compartment_id, compartment_name, query_id) %>% 
  unique()
if (save_intermediate_data) saveRDS(reactome_reactions, file=file.path(data_dir, "reactome_reactions.RDS"))

reactome_mapping_db <- .data %>%
  group_by(modifiers, reaction_id) %>%
  summarize(consumed=list(unique(unlist(reactants))), produced=list(unique(unlist(products)))) %>%
  rename(protein=modifiers)
if (save_intermediate_data) saveRDS(reactome_mapping_db, file=file.path(data_dir, "reactome_mapping_db"))

# save data in dimpl/data directory for distribution with package
usethis::use_data(reactome_reactions, reactome_mapping_db, reactome_data_version, overwrite=TRUE)

# remove large data objects
rm(.data, all_reactions, all_reactions2, reactome_mapping_db, reactome_reactions)

# reset working directory to original
setwd(wd)

