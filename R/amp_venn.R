#' Generates a ggplot2 style venn diagram
#'
#' Calculates the number of "core" OTUs shared by groups given tresholds for how frequent the OTUs should be above a certain abundance. Also returns the average abundance of the OTUs in a particular group.
#'
#' @usage amp_venn(data)
#'
#' @param data (required) A phyloseq object.
#' @param group Group the data based on a sample variable.
#' @param cut_a Abundance cutoff in percent (default: 0.1).
#' @param cut_f Frequency cutoff in percent (default: 80).
#' @param text.size Size of the plotted text.
#' @param output Either "plot" or "complete" (default: "plot").
#' 
#' @return A ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_venn <- function(data, group = NULL,cut_a = 0.1, cut_f = 80, text.size = 5, output = "plot"){
  
  ## Extract the data into separate objects for readability
  abund <- data[["abund"]]
  tax <- data[["tax"]]
  OTU <- data[["tax"]]["OTU"]
  sample <- data[["metadata"]]
  
  ## Test for number of groups
  if (length(levels(sample[,group])) > 3){
    stop(paste("A maximum of 3 levels are supported. The grouping variable contains:",
               paste(levels(sample[,group]), collapse = ", ")))
  }
  
  ## Select grouping variable
  
  colnames(sample)[1] <- "SeqID"
  if (!is.null(group)){
    sample <- sample[,c("SeqID", group)]
    colnames(sample)[2] <- "GRP"  
  } else {
    sample <- data.frame(SeqID = sample[,1], GRP = "Core")
  }
  
  ## Add OTU names to the abundance information
  abund1 <- cbind.data.frame(abund, OTU)
  
  ## Melt the dataframe for subsequent processing
  abund2 <- melt(abund1, id.var = "OTU", value.name= "Abundance", variable.name = "SeqID")
  
  ## Merge sample information with the abundance data
  abund3 <- merge(abund2, sample, by = "SeqID")
  
  ## Add frequent abundant column
  abund4 <- mutate(abund3, 
                   Freq = ifelse(Abundance > cut_a, 1, 0))
  
  ## Evaluate if OTUs are part of the core
  abund5 <- group_by(abund4, OTU, GRP) %>%
    summarise(Abundance = mean(Abundance),
              cFreq = sum(Freq)/n()*100,
              Core = ifelse(sum(Freq)/n()*100 >= cut_f, 1, 0))
  
  ## Convert back into matrix format
  a <- dcast(abund5, OTU~GRP, value.var = "Core")
  
  ################### PLOT #############################
  
  ## 1 group
  if(ncol(a) == 2){
    c_A <- subset(a, a[,2] == 1)
    c_n <- subset(a, a[,2] == 0)
    a_A <- subset(abund5, OTU %in% c_A$OTU) %>% group_by(OTU) %>% 
           summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_n <- subset(abund5, OTU %in% c_n$OTU) %>% group_by(OTU) %>% 
           summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    
    AD <- data.frame(counts = c(A = nrow(c_A), N = nrow(c_n)), 
                     abund = round(c(A = a_A$Sum, N = a_n$Sum), 1))
    
    ## Plot  
    p <- ggplot(data.frame(), aes(x=0, y=0)) +
           annotate("text", x=c(0, 0), y = c(0,-0.5), 
                    label = c(paste(AD[1,1], "\n(", AD[1,2],")", sep = ""), 
                            paste("Non-core: ", AD[2,1]," (", AD[2,2],")", sep = "")), 
                    size = text.size) +
           annotate("text", x=c(0), y = 0.45, label = colnames(a)[2], size = text.size) +
           xlim(-0.65,0.65) +
           ylim(-0.65,0.65) +
           annotate("path", 
                    x=00.4*cos(seq(0,2*pi,length.out=100)),
                    y=0+0.4*sin(seq(0,2*pi,length.out=100))) +
           theme(panel.grid.major=element_blank(), 
                 panel.grid.minor=element_blank(), 
                 axis.text=element_blank(),
                 axis.title=element_blank(),
                 axis.ticks=element_blank(),
                 panel.border=element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 plot.margin = unit(c(0,0,0,0), "mm"))
    
    ## Generate lists of species in each group
    
    ot <- cbind.data.frame(OTU = rownames(tax),tax, abund) %>%
           mutate(Shared = "Non-core") %>%
           mutate(Shared = ifelse(OTU %in% as.character(unique(c_A$OTU)), colnames(a)[2], Shared))
    
  res <- list(plot = p, 
                A = as.character(unique(c_A$OTU)),
                Noncore = as.character(unique(c_n$OTU)),
                Otutable = ot)
    
    names(res)[2] <- colnames(a)[2]
    
  }
  
  ## 2 groups
  if(ncol(a) == 3){
    c_AB <- subset(a, a[,2] == 1 & a[,3] == 1)
    c_A <- subset(a, a[,2] == 1 & a[,3] == 0)
    c_B <- subset(a, a[,2] == 0 & a[,3] == 1)
    c_n <- subset(a, a[,2] == 0 & a[,3] == 0)
    
    a_AB <- subset(abund5, OTU %in% c_AB$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_A <- subset(abund5, OTU %in% c_A$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_B <- subset(abund5, OTU %in% c_B$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_n <- subset(abund5, OTU %in% c_n$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    
    AD <- data.frame(counts = c(A = nrow(c_A), AB = nrow(c_AB), B = nrow(c_B), N = nrow(c_n)), 
                     abund = round(c(A = a_A$Sum, AB = a_AB$Sum, B = a_B$Sum, N = a_n$Sum), 1))
    
    p <- ggplot(data.frame(), aes(x=c(-0.2,0.2), y=0)) +
           annotate("text", x=c(-0.4, 0, 0.4, 0), y = c(0,0,0,-0.5), 
                    label = c(paste(AD[1:3,1], "\n(", AD[1:3,2],")", sep = ""), 
                              paste("Non-core: ", AD[4,1]," (", AD[4,2],")", sep = "")), 
                    size = text.size) +
           annotate("text", x=c(-0.2, 0.2), y = 0.45, label = colnames(a)[2:3], size = text.size) +
           xlim(-0.65,0.65) +
           ylim(-0.65,0.65) +
           annotate("path",
                     x=0.2+0.4*cos(seq(0,2*pi,length.out=100)),
                     y=0+0.4*sin(seq(0,2*pi,length.out=100))) +
           annotate("path",
                     x=-0.2+0.4*cos(seq(0,2*pi,length.out=100)),
                     y=0+0.4*sin(seq(0,2*pi,length.out=100))) +
           theme(panel.grid.major=element_blank(), 
                 panel.grid.minor=element_blank(), 
                 axis.text=element_blank(),
                 axis.title=element_blank(),
                 axis.ticks=element_blank(),
                 panel.border=element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 plot.margin = unit(c(0,0,0,0), "mm"))
    
    ot <- cbind.data.frame(OTU = rownames(tax),tax, abund) %>%
      mutate(Shared = "Non-core") %>%
      mutate(Shared = ifelse(OTU %in% as.character(unique(c_A$OTU)), colnames(a)[2], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_B$OTU)), colnames(a)[3], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_AB$OTU)), "Core", Shared))
    
    res <- list(plot = p, 
                A = as.character(unique(c_A$OTU)),
                B = as.character(unique(c_B$OTU)),
                Core = as.character(unique(c_AB$OTU)),
                Noncore = as.character(unique(c_n$OTU)),
                Otutable = ot)
    
    names(res)[2:4] <- c(colnames(a)[2],
                         colnames(a)[3], 
                         "Core")
    
  }
  
  ## 3 groups
  if(ncol(a) == 4){
    c_ABC <- subset(a, a[,2] == 1 & a[,3] == 1 & a[,4] == 1)
    c_AB  <- subset(a, a[,2] == 1 & a[,3] == 1 & a[,4] == 0)
    c_AC  <- subset(a, a[,2] == 1 & a[,3] == 0 & a[,4] == 1)
    c_BC  <- subset(a, a[,2] == 0 & a[,3] == 1 & a[,4] == 1)
    c_A   <- subset(a, a[,2] == 1 & a[,3] == 0 & a[,4] == 0)
    c_B   <- subset(a, a[,2] == 0 & a[,3] == 1 & a[,4] == 0)
    c_C   <- subset(a, a[,2] == 0 & a[,3] == 0 & a[,4] == 1)
    c_n   <- subset(a, a[,2] == 0 & a[,3] == 0 & a[,4] == 0)
    
    
    a_ABC <- subset(abund5, OTU %in% c_ABC$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_AB  <- subset(abund5, OTU %in% c_AB$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_AC  <- subset(abund5, OTU %in% c_AC$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_BC  <- subset(abund5, OTU %in% c_BC$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_A   <- subset(abund5, OTU %in% c_A$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_B   <- subset(abund5, OTU %in% c_B$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_C   <- subset(abund5, OTU %in% c_C$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    a_n   <- subset(abund5, OTU %in% c_n$OTU) %>% group_by(OTU) %>% 
      summarise(Mean = mean(Abundance)) %>% as.data.frame() %>% summarise(Sum = sum(Mean))
    
    AD <- data.frame(counts = c(ABC = nrow(c_ABC), AB = nrow(c_AB),AC = nrow(c_AC), BC = nrow(c_BC), 
                                A = nrow(c_A), B = nrow(c_B), C = nrow(c_C), N = nrow(c_n)), 
                     abund = round(c(ABC = a_ABC$Sum, AB = a_AB$Sum, AC = a_AC$Sum, BC = a_BC$Sum, 
                                     A = a_A$Sum, B = a_B$Sum, C = a_C$Sum, N = a_n$Sum), 1))
    
    p <- ggplot(data.frame(), aes(x=c(-0.2,0.2), y=0)) +
      annotate("text", x=c(0, 0, -0.25, 0.25, -0.4, 0.4, 0, 0.5), y = c(0.05, 0.4, -0.05, -0.05, 0.3, 0.3, -0.3, -0.6), 
               label = c(paste(AD[1:7,1], "\n(", AD[1:7,2],")", sep = ""), 
                         paste("Non-core:\n", AD[8,1]," (", AD[8,2],")", sep = "")), 
               size = text.size) +
      annotate("text", x=c(-0.2, 0.2, 0), y = c(0.65, 0.65, -0.65), label = colnames(a)[2:4], size = text.size) +
      xlim(-0.65,0.65) +
      ylim(-0.65,0.65) +
      annotate("path",
               x=0.2+0.4*cos(seq(0,2*pi,length.out=100)),
               y=0.2+0.4*sin(seq(0,2*pi,length.out=100))) +
      annotate("path",
               x=-0.2+0.4*cos(seq(0,2*pi,length.out=100)),
               y=0.2+0.4*sin(seq(0,2*pi,length.out=100))) +
      annotate("path",
               x=0+0.4*cos(seq(0,2*pi,length.out=100)),
               y=-0.2+0.4*sin(seq(0,2*pi,length.out=100))) +  
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            axis.text=element_blank(),
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            panel.border=element_blank(),
            panel.background = element_blank(),
            legend.key = element_blank(),
            plot.margin = unit(c(0,0,0,0), "mm"))
    
    ot <- cbind.data.frame(OTU = rownames(tax),tax, abund) %>%
      mutate(Shared = "Non-core") %>%
      mutate(Shared = ifelse(OTU %in% as.character(unique(c_ABC$OTU)), "Core", Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_AB$OTU)), paste(colnames(a)[2], colnames(a)[3], sep = "_"), Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_AC$OTU)), paste(colnames(a)[2], colnames(a)[4], sep = "_"), Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_BC$OTU)), paste(colnames(a)[3], colnames(a)[4], sep = "_"), Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_A$OTU)), colnames(a)[2], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_B$OTU)), colnames(a)[3], Shared),
             Shared = ifelse(OTU %in% as.character(unique(c_C$OTU)), colnames(a)[4], Shared))
    
    res <- list(plot = p, 
                Core = as.character(unique(c_ABC$OTU)),
                AB = as.character(unique(c_AB$OTU)),
                AC = as.character(unique(c_AC$OTU)),
                BC = as.character(unique(c_BC$OTU)),
                A = as.character(unique(c_A$OTU)),
                B = as.character(unique(c_B$OTU)),
                C = as.character(unique(c_C$OTU)),
                Noncore = as.character(unique(c_n$OTU)),
                Otutable = ot)
    
    names(res)[2:8] <- c("Core",
                         paste(colnames(a)[2], colnames(a)[3], sep = "_"),
                         paste(colnames(a)[2], colnames(a)[4], sep = "_"),
                         paste(colnames(a)[3], colnames(a)[4], sep = "_"),
                         colnames(a)[2],
                         colnames(a)[3],
                         colnames(a)[4]
    )
    
  }  
  
  ## Export data
  if(output == "complete"){ return(res) }
  if(output == "plot"){ return(p) }
}
