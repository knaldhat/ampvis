amp_subset <- function(list, ...) {
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
  
  #subset metadata based on ... and only keep columns in otutable matching the rows in the subsetted metadata
  newmetadata <- subset(metadata, ...)
  newabund <- abund[, rownames(newmetadata), drop=FALSE]
  
  #return a new list
  newlist <- list(abund = newabund, tax = tax, metadata = newmetadata)
  return(newlist)
}
