#' Load data as OTU table and convert to a phyloseq object.
#'
#' Load data as OTU table and convert to a phyloseq object.
#'
#' @usage amp_load(otutable, metadata)
#'
#' @param otutable (required) A OTU table generated from workflow scripts v.4+.
#' @param metadata (required) A metadata file with sample names in first column.
#' @param refseq Reference sequences for all OTUs.
#' @param rarefy Rarefy all samples to the same sequencing depth.
#' 
#' @return A phyloseq object.
#' 
#' @export
#' @import phyloseq
#' @import dplyr
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_load <- function(otutable, metadata, rarefy = NULL, percent = FALSE){
  # Remove whitespace from the otutable as this will break the structure of the taxonomy
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  otutable$Kingdom<-trim(as.character(otutable$Kingdom))
  otutable$Phylum<-trim(as.character(otutable$Phylum))
  otutable$Class<-trim(as.character(otutable$Class))
  otutable$Order<-trim(as.character(otutable$Order))
  otutable$Family<-trim(as.character(otutable$Family))
  otutable$Genus<-trim(as.character(otutable$Genus))
  otutable$Species<-trim(as.character(otutable$Species))
  
  #metadata: order rows by rownames
  rownames(metadata) <- metadata[,1]
  metadata = suppressWarnings(as.data.frame(as.matrix(metadata)))
  metadata <- metadata[order(rownames(metadata)), ]
  
  #abund: all columns from otutable except the last 7 to numeric and order rows by rownames:
  abund <- as.data.frame(otutable[,1:(ncol(otutable) - 7)])/1
  abund <- abund[order(rownames(abund)),order(colnames(abund))]
  #abundances to percent
  if(percent == TRUE) {
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  #tax: the last 7 columns from otutable to factor, order rows by rownames and order columns by taxonomic rank(not alphabetically)
  tax <- data.frame(otutable[, (ncol(otutable) - 6):ncol(otutable)] %>% transform(as.factor) 
                    ,OTU = rownames(otutable))
  tax <- tax[order(rownames(tax)), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
  
  #data: return the data in a combined list
  data <- list(abund = abund, tax = tax, metadata = metadata)
  
  #rarefy function
  #  if(!is.null(rarefy)){data <- rarefy_even_depth(data, sample.size = rarefy, rngseed = 712)}metadata <- metadata[order(rownames(metadata)), ]
  
  #check if metadata and otutable match
  if(!all(rownames(data$metadata) %in% colnames(data$abund))) {
    stop("The sample names in metadata do not match those in otutable")
  }
  return(data)
}
