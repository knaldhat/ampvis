amp_subset_taxa <- function(list, ...) {
  #Check the data first
  if(!is.list(list) | 
     !any(names(list) == "abund") | 
     !any(names(list) == "tax") | 
     !any(names(list) == "metadata") | 
     !is.data.frame(list[["abund"]]) |
     !is.data.frame(list[["tax"]]) |
     !is.data.frame(list[["metadata"]])
  ) {
    stop("The data must be a list with three dataframes named abund, tax and metadata")
  }
  
  #extract data from the list
  metadata <- list$metadata
  abund <- list$abund
  tax <- list$tax
  
  #subset tax table based on ... and only keep rows in abund and metadata matching the rows in the subsetted tax table
  newtax <- subset(tax, ...)
  newabund <- abund[rownames(newtax), , drop=FALSE]
  newmetadata <- metadata[colnames(newabund), , drop=FALSE]
  
  #return a new list
  newlist <- list(abund = newabund, tax = newtax, metadata = newmetadata)
  return(newlist)
}