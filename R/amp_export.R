#' Export sequences from a phyloseq object
#'
#' Export sequences from a phyloseq object
#'
#' @usage amp_export(data)
#'
#' @param data (required) A data list loaded with amp_load.
#' @param file Name of the file containing the exported sequences.
#' @param tax Add taxonomic strings to the output (default: T).
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_export <- function(data, file = "exported_sequences.fa", tax = T){
  
  t <- data[["refseq"]]
  
  if (tax == T){
    tax <- as.data.frame(data[["tax"]][,1:7])
    tax_s <- data.frame(lapply(tax, as.character), stringsAsFactors=FALSE)
    df_args <- c(tax_s, sep="; ")
    tax_sf <- do.call(paste, df_args)
    names(t) <- paste(names(t), tax_sf, sep = "; ")
  }
  
  writeXStringSet(t, file = file)
}
