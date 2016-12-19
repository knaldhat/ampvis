#' Load data as OTU table and convert to a phyloseq object.
#'
#' Load data as OTU table and convert to a phyloseq object.
#'
#' @usage amp_load(otutable, metadata)
#'
#' @param otutable (required) A OTU table generated from workflow scripts v.4+.
#' @param metadata (required) A metadata file with sample names in first column.
#' @param refseq Reference sequences for all OTUs. Must be loaded with readDNAStringSet() from the biostrings package.
#' @param rarefy Rarefy all samples to the same sequencing depth.
#' @param percent Transform abundances from raw counts to percentages.
#' 
#' @return A phyloseq object.
#' 
#' @export
#' @import dplyr
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_load <- function(otutable, metadata, refseq = NULL, rarefy = NULL, percent = FALSE){
  #Must be data frames
  otutable <- as.data.frame(otutable)
  metadata <- as.data.frame(metadata)
  
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
  metadata = suppressWarnings(as.data.frame(as.matrix(metadata)))
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[order(rownames(metadata)), ]
  
  #Only alphanumeric characters in metadata column names, replace others with _
  colnames(metadata) <- str_replace_all(colnames(metadata), "[^[:alnum:]]", "_")
  
  #abund: all columns from otutable except the last 7 to numeric and order rows by rownames:
  abund <- as.data.frame(otutable[,1:(ncol(otutable) - 7)])/1
  abund <- abund[order(rownames(abund)),order(colnames(abund))]
  
  #rarefy function
  if(!is.null(rarefy)){
    abund <- rarefy(abund, sample = rarefy)
  }
  
  #abundances to percent, must be done AFTER rarefy
  if(percent == TRUE) {
    abund <- as.data.frame(sapply(abund, function(x) x/sum(x)*100))
  }
  
  #tax: the last 7 columns from otutable to factor, order rows by rownames and order columns by taxonomic rank(not alphabetically)
  tax <- data.frame(otutable[, (ncol(otutable) - 6):ncol(otutable)] 
                    ,OTU = rownames(otutable))
  tax <- tax[order(rownames(tax)), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")]
  
  #data: return the data in a combined list w or w/o refseq. Load
  if(!is.null(refseq) & class(refseq) == "DNAStringSet") {
    data <- list(abund = abund, tax = tax, metadata = metadata, refseq = refseq)
  } else if(!is.null(refseq) & !class(refseq) == "DNAStringSet") {
    stop("The reference sequences must be loaded with readDNAStringSet() from the biostrings package.")
  } else if(is.null(refseq)) {
    data <- list(abund = abund, tax = tax, metadata = metadata)
  }
  
  #check if metadata and otutable match
  if(!all(rownames(data$metadata) %in% colnames(data$abund))) {
    stop("The sample names in metadata do not match those in otutable")
  }
  return(data)
}
