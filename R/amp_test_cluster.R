#' Tests if there is a significant difference between samples.
#'
#' A nice long description
#'
#' @usage amp_test_cluster(data, group)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param group (required) The group variables to be tested.
#' @param method The distance measure (default: "bray").
#' @param plot.label A vector specifying the sample data to be used as labels.
#' @param plot.color A vector specifying the sample data to be used as colors.
#' @param plot.theme Chose different standard layouts choose from "normal" or "clean" (default: "normal").
#' 
#' @return A number of things.
#' 
#' @export
#' @import ggplot2
#' @import ape
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_test_cluster <- function(data, group, method = "bray", plot.label = NULL, plot.color = NULL, plot.theme = "none"){
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  sample <- data[["metadata"]]
  
  ## Calculate beta-diversity
  betad <- vegdist(x = abund, method = method)
  
  ## Use Adonis to test for overall differences
  test_formula <- as.formula(paste("betad", paste(group, collapse="+"), sep="~ "))  
  res_adonis <- adonis(test_formula, sample)  
  
  ## Cluster the samples 
  hc <- hclust(betad)
    
  ## Add custom colors and labels to the cluster data
  if (!is.null(plot.label)){
    if (length(plot.label) > 1){
      hc$labels <- apply(sample[,plot.label], 1, paste, collapse= "; ")
    } else{
      hc$labels <- as.character(sample[,plot.label])
      names(hc$labels) <- rownames(sample)
    }
  } else {names(hc$labels) <- rownames(sample)}  
  
  hc_d <- dendro_data(as.dendrogram(hc))
    
    
  if (!is.null(plot.color)){    
    if (length(plot.color) > 1){
      tcol <- data.frame(sample = rownames(sample), group = as.factor(apply(sample[,plot.color], 1, paste, collapse= "; ")))
    } else{
      tcol <- data.frame(sample = rownames(sample), group = sample[,plot.color])
    } 
    ncols <- length(levels(tcol$group))
    
    hues = seq(15, 375, length=ncols+1)
    pal <- hcl(h=hues, l=65, c=100)[1:ncols]    
    tcol$color <- pal[as.numeric(tcol$group)]
    tcol <- tcol[rownames(label(hc_d)),]
    hc_d$labels$color <- tcol$color
    hc_d$labels$group <- factor(tcol$group, levels = unique(tcol$group))
  } else {
    cols <- "black"
  }
  
  ## Plot clusters
  
  p1 <- ggplot(data = segment(hc_d)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    coord_flip() +
    scale_x_discrete(labels=label(hc_d)$label) +
    ylab(paste("Distance (beta diversity = ", method, ")", sep ="")) +
    theme(axis.text.y = element_text(color = hc_d$labels$color),
          axis.title.y = element_blank())
  
  if (!is.null(plot.color)){ 
  p1 <- p1 + geom_point(data=hc_d$label, aes(x = x, y = y, color = group), inherit.aes =F, alpha = 0) + 
    scale_color_manual(labels = rev(levels(hc_d$label$group)), values = rev(unique(hc_d$label$color)), name = "") +
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))
  }
  
  if(plot.theme == "clean"){
    p1 <- p1 + theme(text = element_text(size = 8, color = "black"),
                     axis.text = element_text(size = 8, color = "black"),
                     axis.text.y = element_text(hjust = 1),
                     plot.margin = unit(c(0,0,0,0), "mm"),
                     axis.line = element_line(color = "black"),
                     axis.line.y = element_blank(),
                     panel.grid = element_blank(),
                     legend.key = element_blank(),
                     axis.ticks = element_line(color = "black"),
                     axis.ticks.y = element_blank(),
                     panel.background = element_blank()
    )
  }
  
  return(list(betad = betad, adonis = res_adonis, hc = hc, plot_cluster = p1, temp  = hc_d))  
}
