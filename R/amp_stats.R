#' Calculate basic statistics for each sample
#'
#' Calculate basic statistics for each sample and combine it with metadata.
#'
#' @usage amp_stats(data)
#'
#' @param data (required) A phyloseq object.
#' @param measure Alpha-diversity measures to be included (default:observed).
#' 
#' @export
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

function (data, split = TRUE, measures = "Observed") 
{
  Reads <- rowSums(data[["abund"]])
  dmeta <- data[["metadata"]]
  if (!any(data[["abund"]] == 1)) {
    warning("The data you have provided does not have\n", 
            "any singletons. This is highly suspicious. Results of richness\n", 
            "estimates (for example) are probably unreliable, or wrong, if you have already\n", 
            "trimmed low-abundance taxa from the data.\n", "\n", 
            "We recommended that you find the un-trimmed data and retry.")
  }
  if (!split) {
    OTU <- rowSums(data[["abund"]])
  }
  else if (split) {
    OTU <- as(data[["abund"]], "matrix")
  }
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
                "InvSimpson", "Fisher")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", 
                        "simpson", "invsimpson", "fisher")
  if (is.null(measures)) {
    measures = as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% 
                                                            measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  outlist = vector("list")
  estimRmeas = c("Chao1", "Observed", "ACE")
  if (any(estimRmeas %in% measures)) {
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = diversity(OTU, 
                                                   index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = diversity(OTU, 
                                                   index = "simpson")))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = diversity(OTU, 
                                                      index = "invsimpson")))
  }
  if ("Fisher" %in% measures) {
    fisher = tryCatch(fisher.alpha(OTU, se = TRUE), warning = function(w) {
      warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
      suppressWarnings(fisher.alpha(OTU, se = TRUE)[, 
                                                    c("alpha", "se")])
    })
    if (!is.null(dim(fisher))) {
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    }
    else {
      outlist <- c(outlist, Fisher = list(fisher))
    }
  }
  rich = do.call("cbind", outlist)
  namechange = intersect(colnames(rich), names(renamevec))
  colnames(rich)[colnames(rich) %in% namechange] <- renamevec[namechange]
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, 
                   colnames(rich), ignore.case = TRUE)
  rich = rich[, sort(unique(unlist(colkeep))), drop = FALSE]
  rich <- as.data.frame(rich)
  combined <- cbind.data.frame(dmeta, Reads, rich) %>% arrange(Reads)
  return(combined)
}

