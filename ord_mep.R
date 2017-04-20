#A wrapper around the vegan package to make beautiful ggplot2 ordination plots 
#to analyse and compare microbial community compositions. The input data must be a list loaded with
#amp_load() from the ampvis package, but any OTU table-like matrix (OTU's in rows,
#sampleID's in columns and abundances in the corresponding cells, aka a contingency table) alongside metadata 
#for the same samples can be used. Simply choose an ordination type and a plot is returned.
#Required packages:
#dplyr, ampvis, vegan, ape
ord_mep <- function(
  ##### data, ordination type, transformations etc... #####
  datalist, #Data list loaded with amp_load() containing the elements: abund, metadata and tax.
  filter_lowabund = 0.1, #Remove low abundant OTU's across all samples below this threshold in percent. Recommended minimum: 0.1%.
  type = "PCA", #Ordination type; Principal Components Analysis(PCA), Redundancy Analysis(RDA), non-metric Multidimensional Scaling(NMDS), metric Multidimensional Scaling(MMDS, aka PCoA), Correspondence Analysis(CA), Canonical Correspondence Analysis(CCA) or Detrended Correspondence Analysis(DCA).
  metric = NULL, #Distance metric used for the distance-based ordination methods (nMDS/PCoA), any of the following:
  #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "jsd" (Jensen-Shannon Divergence),
  #"altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"
  #or simply "none" or "sqrt". See details in ?vegdist, and for JSD http://enterotype.embl.de/enterotypes.html.
  transform = NULL, #Transform the abundance table with the decostand() function, fx "normalize", "chi.square", "hellinger" or "sqrt", see details in ?decostand. Using the hellinger transformation is always a good choice and is recommended for PCA/RDA/nMDS/PCoA to obtain a more ecologically meaningful result (learn about the double-zero problem).
  constrain = NULL, #Variable(s) in the metadata for constrained analyses (RDA and CCA). Multiple variables can be provided by a vector, fx c("Year", "Temperature"), but keep in mind the more variables the more the result will be similar to unconstrained analysis.
  
  ##### Geometric options, colors, labels etc... #####
  x_axis = 1, #Which axis from the ordination results to plot as the first axis. Have a look at the $screeplot with output = "detailed" to validate axes.
  y_axis = 2, #Which axis from the ordination results to plot as the second axis. Have a look at the $screeplot with output = "detailed" to validate axes.
  color_by = NULL, #Color points by a variable in the metadata.
  color_order = NULL, #Order the colors in color_by by a vector.
  label_by = NULL, #Label by a variable in the metadata.
  shape_by = NULL, #Shape by a variable in the metadata.
  colorframe = FALSE, #Frame the points with a polygon colored by the color_by argument.
  opacity = 0.8, #Sample points and colorframe opacity, 0:invisible, 1:opaque.
  trajectory = NULL, #Make a trajectory between sample points by a variable in the metadata
  trajectory_group = trajectory, #Make a trajectory between sample points by the trajectory argument, but within individual groups.
  
  label_samples_by = NULL, #Label sample points by a variable in the metadata
  sample_label_size = 4, #Sample labels text size
  sample_label_segment_color = "black", #Sample labels repel-segment color

  #### Species options ####
  plot_species_points = FALSE, #Plot species points.
  nspecies_labels = 0, #Number of most extreme species labels to plot.
  nspecies_labels_size = 3, #Size of the species text labels.
  label_species_by = "Genus", #Taxonomic level by which to label the species points.
  scale_species = FALSE, #rescale species
  species_point_size = 2, #Size of the species points
  species_point_shape = 1, #The shape of the points, fx 1 for hollow circles or 20 for dots
  
  ##### fit environmental variables #####
  envfit_factor = NULL, 
  envfit_numeric = NULL,
  envfit_signif_level = 0.001, 
  envfit_textsize = 3,
  envfit_color = "darkred",
  envfit_numeric_arrows_scale = 1,
  envfit_show = TRUE,
  
  ##### miscellaneous #####
  tax.empty = "best", # "OTU" or "best"
  
  ##### output #####
  output = "plot", #"plot" or "detailed"; output as list with additional information(model, scores, inputmatrix etc) or just the plot
  ... #Pass additional arguments to the vegan ordination functions, fx rda(...), cca(...), metaMDS(...), see vegan help
) {
  #Check the data
  datalist <- amp_rename(data = datalist, tax.empty = tax.empty)
  
  #First transform to percentages
  abund_pct <- as.data.frame(sapply(datalist$abund, function(x) x/sum(x) * 100))
  rownames(abund_pct) <- rownames(datalist$abund) #keep rownames
  datalist$abund <- abund_pct
  
  #Then filter low abundant OTU's where ALL samples have below the threshold set with filter_lowabund in percent
  datalist$abund <- datalist$abund[!apply(datalist$abund, 1, function(row) all(row <= filter_lowabund)),] #remove low abundant OTU's 
  datalist$tax <- datalist$tax[rownames(datalist$abund),] #same with taxonomy
  datalist$metadata <- datalist$metadata[colnames(datalist$abund),] #same with metadata
  
  #to fix user argument characters, so fx PCoA/PCOA/pcoa are all valid
  type <- tolower(type)
  output <- tolower(output)
  
  if(!is.null(metric)) {
    metric <- tolower(metric)
  } else if(is.null(metric)) {
    if(type == "nmds" | type == "mmds" | type == "pcoa" | type == "dca") {
      warning("No distance metric selected, using raw data. If this is not deliberate, please provide one with the argument: metric")
    }
    metric <- "none"
  }
  
  #data transformation with decostand()
  if(!is.null(transform)) {
    transform <- tolower(transform)
    if(transform == "sqrt") {
      datalist$abund <- t(sqrt(t(datalist$abund)))
    } else {
      datalist$abund <- t(decostand(t(datalist$abund), method = transform))
    }
  } 
  
  #Calculate distance matrix with vegdist()
  if (metric == "none") {
    inputmatrix <- t(datalist$abund)
  } else if(metric == "jsd") {
    #This is based on http://enterotype.embl.de/enterotypes.html
    #Abundances of 0 will be set to the pseudocount value to avoid 0-value denominators
    dist.JSD <- function(inMatrix, pseudocount=0.000001) {
      KLD <- function(x,y) sum(x *log(x/y))
      JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
      matrixColSize <- length(colnames(inMatrix))
      matrixRowSize <- length(rownames(inMatrix))
      colnames <- colnames(inMatrix)
      resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
      
      inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
      
      for(i in 1:matrixColSize) {
        for(j in 1:matrixColSize) { 
          resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                                 as.vector(inMatrix[,j]))
        }
      }
      colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
      as.dist(resultsMatrix)->resultsMatrix
      attr(resultsMatrix, "method") <- "dist"
      return(resultsMatrix) 
    }
    inputmatrix <- dist.JSD(datalist$abund)
  } else if(any(metric == c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"))) {
    inputmatrix <- vegdist(t(datalist$abund), method = metric)
  }
  #################################### end of block ####################################
  
  #Generate data depending on the chosen ordination type
  if(type == "pca") {
    #make the model
    model <- rda(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("PC", x_axis)
    y_axis_name <- paste0("PC", y_axis)
    
    #Calculate the amount of inertia explained by each axis
    totalvar <- round(model$CA$eig/model$tot.chi * 100, 1)
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "rda") {
    if(is.null(constrain)) 
      stop("Argument constrain must be provided when performing constrained/canonical analysis.")
    #make the model
    codestring <- paste0("rda(inputmatrix~", paste(constrain, collapse = "+"), ", datalist$metadata, ...)") #function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <-  eval(parse(text = codestring))
    
    #axes depend on the results
    x_axis_name <- paste0("RDA", x_axis)
    if (model$CCA$rank <= 1){
      y_axis_name <- "PC1"
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig/model$CA$tot.chi * 100, 1)) #UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("RDA", y_axis)
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig/model$CCA$tot.chi * 100, 1)) #constrained of total constrained space
    }
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "nmds") {
    #make the model
    model <- metaMDS(inputmatrix, trace = FALSE, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("NMDS", x_axis)
    y_axis_name <- paste0("NMDS", y_axis)
    
    #Calculate species- and site scores
    #Speciesscores may not be available with MDS
    sitescores <- scores(model, display = "sites")
    if(!length(model$species) > 1) {
      speciesscores <- warning("Speciesscores are not available.")
    } else {
      speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
    }
  } else if(type == "mmds" | type == "pcoa") {
    #make the model
    model <- pcoa(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("PCo", x_axis)
    y_axis_name <- paste0("PCo", y_axis)
    
    #Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$values$Relative_eig * 100, 1)
    names(totalvar) <- c(paste0("PCo", seq(1:length(totalvar))))
    
    #Calculate species- and site scores
    #Speciesscores are not available with pcoa
    sitescores <- as.data.frame(model$vectors)
    colnames(sitescores) <- c(paste0("PCo", seq(1:length(sitescores))))
    speciesscores <- warning("Speciesscores are not available.")
  } else if(type == "ca") {
    #make the model
    model <- cca(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("CA", x_axis)
    y_axis_name <- paste0("CA", y_axis)
    
    #Calculate the percentage of eigenvalues explained by the axes
    totalvar <- round(model$CA$eig/model$tot.chi * 100, 1)
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "cca") {
    if(is.null(constrain)) 
      stop("Argument constrain must be provided when performing constrained/canonical analysis.")
    #make the model
    codestring <- paste0("cca(inputmatrix~", paste(constrain, collapse = "+"), ", datalist$metadata, ...)") #function arguments written in the format "rda(x ~ y + z)" cannot be directly passed to rda(), now user just provides a vector
    model <-  eval(parse(text = codestring))
    
    #axes depend on the results
    x_axis_name <- paste0("CCA", x_axis)
    if (model$CCA$rank <= 1){
      y_axis_name <- "CA1"
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CA$eig/model$CA$tot.chi * 100, 1)) #UNconstrained of total UNconstrained space
    } else if (model$CCA$rank > 1) {
      y_axis_name <- paste0("CCA", y_axis)
      
      #Calculate the amount of inertia explained by each axis
      totalvar <- c(round(model$CCA$eig/model$tot.chi * 100, 1), #constrained of total space
                    round(model$CA$eig/model$tot.chi * 100, 1)  #UNconstrained of total space
      )
      constrainedvar <- c(round(model$CCA$eig/model$CCA$tot.chi * 100, 1)) #constrained of total constrained space
    }
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  } else if(type == "dca") {
    #make the model
    model <- decorana(inputmatrix, ...)
    
    #axis (and data column) names
    x_axis_name <- paste0("DCA", x_axis)
    y_axis_name <- paste0("DCA", y_axis)
    
    #Calculate the percentage of eigenvalues explained by the axes
    #totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)
    
    #Calculate species- and site scores
    sitescores <- scores(model, display = "sites", choices = c(x_axis, y_axis))
    speciesscores <- scores(model, display = "species", choices = c(x_axis, y_axis))
  }
  #################################### end of block ####################################
  
  #Make data frames for ggplot
  dsites <- cbind.data.frame(datalist$metadata, sitescores)
  
  if (!is.null(color_order)) {
    dsites[, color_by] <- factor(dsites[, color_by], levels = color_order)
  }
  
  if(length(speciesscores) > 1) {
    dspecies <- merge(data.frame(speciesscores, OTU = rownames(speciesscores)), datalist$tax, by.x = "OTU")
    dspecies$dist <- dspecies[, x_axis_name]^2 + dspecies[, y_axis_name]^2
    dspecies <- arrange(dspecies, desc(dist))
    rownames(dspecies) <- dspecies$OTU
    if (scale_species == TRUE) {
      maxx <- max(abs(dsites[, x_axis_name]))/max(abs(dspecies[,x_axis_name]))
      dspecies[, x_axis_name] <- dspecies[, x_axis_name] * maxx * 0.8
      maxy <- max(abs(dsites[, y_axis_name]))/max(abs(dspecies[,y_axis_name]))
      dspecies[, y_axis_name] <- dspecies[, y_axis_name] * maxy * 0.8
    }
  } else {
    dspecies = NULL
  }
  
  #Generate a nice ggplot with the coordinates from scores
  plot <- ggplot(dsites,
                 aes_string(x = x_axis_name,
                            y = y_axis_name,
                            color = color_by,
                            shape = shape_by
                 )
  ) +
    geom_point(size = 2, alpha = opacity) + 
    theme_minimal() +
    theme(axis.line = element_line(colour = "black", size = 0.5))
  
  #Only eigenvalue-based ordination methods can be displayed with % on axes
  if(type == "pca" | type == "ca" | type == "pcoa" | type == "mmds") {
    plot <- plot +
      xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) + 
      ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "%]", sep = ""))
  } else if(type == "rda" | type == "cca") {
    if(model$CCA$rank > 1) {
      plot <- plot + 
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "% / ", constrainedvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    } else if(model$CCA$rank <= 1) {
      plot <- plot + 
        xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) +
        ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "% / ", constrainedvar[y_axis_name], "%]", sep = ""))
    }
  } else if (type == "nmds") {
    plot <- plot +
      annotate("text", size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste0("Stress value = ", round(model$stress, 3)))
  }
  
  #Plot species points
  if (plot_species_points == TRUE) {
    plot <- plot + 
      geom_point(data = dspecies,
                 color = "darkgrey",
                 shape = species_point_shape,
                 size = species_point_size,
                 alpha = 0.8,
                 aes(text = paste0(label_species_by, ": ", dspecies$Genus)) #To allow hovering points with ggplotly()
      )
  }
  
  #Generate a color frame around the chosen color group
  if(colorframe == TRUE) {
    if(is.null(color_by)) stop("Please provide the argument color_by")
    splitData <- split(dsites, dsites[, color_by]) %>% 
      lapply(function(df) {
        df[chull(df[, x_axis_name], df[, y_axis_name]), ]
      })
    hulls <- do.call(rbind, splitData)
    plot <- plot + geom_polygon(data = hulls, aes_string(fill = color_by, group = color_by), alpha = 0.2*opacity)
  }
  
  #Plot text labels
  if (!is.null(label_by)) {
    temp <- data.frame(group = dsites[, label_by], 
                       x = dsites[, x_axis_name],
                       y = dsites[, y_axis_name]) %>% 
      group_by(group) %>%
      summarise(cx = mean(x), cy = mean(y)) %>% 
      as.data.frame()
      temp2 <- merge(dsites, temp,
                     by.x = label_by, 
                     by.y = "group")
      temp3 <- temp2[!duplicated(temp2[, label_by]), ]
      plot <- plot +
        geom_text_repel(data = temp3,
                        aes_string(x = "cx", 
                                   y = "cy",
                                   label = label_by),
                        size = 3,
                        color = "black",
                        fontface = 2
      )
  }
  
  #Trajectory
  if (!is.null(trajectory)) {
    traj <- dsites[order(dsites[, trajectory]), ]
    plot <- plot + geom_path(data = traj, aes_string(group = trajectory_group))
  }
  
  #Sample point labels
  if(!is.null(label_samples_by)) {
    plot <- plot + geom_text_repel(aes_string(label = label_samples_by), 
                             size = sample_label_size, color = "grey40", segment.color = sample_label_segment_color)
  }
  
  #Plot species labels
  if (nspecies_labels > 0) {
    plot <- plot +
      geom_text_repel(data = dspecies[1:nspecies_labels,],
                      aes_string(x = x_axis_name,
                                 y = y_axis_name,
                                 label = label_species_by
                      ),
                      colour = "red",
                      size = nspecies_labels_size,
                      fontface = 4,
                      inherit.aes = FALSE
      )
  }
  
  ######## Fit environmental variables ########
  # Categorial fitting
  if(!is.null(envfit_factor)) {
    evf_factor_model <- envfit(model,
                               datalist$metadata[,envfit_factor, drop = FALSE],
                               permutations = 999,
                               choices = c(x_axis_name, y_axis_name)
    )
    evf_factor_data <- data.frame(Name = rownames(evf_factor_model$factors$centroids),
                                  Variable = evf_factor_model$factors$var.id,
                                  evf_factor_model$factors$centroids,
                                  pval = evf_factor_model$factors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_factor_data) > 0 & envfit_show == TRUE) {
      plot <- plot + geom_text_repel(data = evf_factor_data,
                                     aes_string(x = x_axis_name,
                                                y = y_axis_name,
                                                label = "Name"),
                                     colour = envfit_color, 
                                     inherit.aes = FALSE,
                                     size = envfit_textsize,
                                     fontface = "bold"
      )
    }
    if (nrow(evf_factor_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.")
    }
  } else {
    evf_factor_model <- NULL
  }
  
  # Numerical fitting
  if (!is.null(envfit_numeric)) {
    evf_numeric_model <- envfit(model,
                                datalist$metadata[,envfit_numeric, drop = FALSE],
                                permutations = 999,
                                choices = c(x_axis_name, y_axis_name)
    )
    evf_numeric_data <- data.frame(Name = rownames(evf_numeric_model$vectors$arrows),
                                   evf_numeric_model$vectors$arrows * sqrt(evf_numeric_model$vectors$r) * envfit_numeric_arrows_scale,
                                   pval = evf_numeric_model$vectors$pvals
    ) %>% subset(pval <= envfit_signif_level)
    if (nrow(evf_numeric_data) > 0 & envfit_show == TRUE) {
      plot <- plot + geom_segment(data = evf_numeric_data,
                                  aes_string(x = 0,
                                             xend = x_axis_name,
                                             y = 0,
                                             yend = y_axis_name
                                  ),
                                  arrow = arrow(length = unit(3, "mm")),
                                  colour = "darkred",
                                  size = 1,
                                  inherit.aes = FALSE) + 
        geom_text(data = evf_numeric_data,
                  aes_string(x = x_axis_name,
                             y = y_axis_name,
                             label = "Name"),
                  colour = envfit_color,
                  inherit.aes = FALSE,
                  size = envfit_textsize,
                  hjust = 1.2,
                  vjust = 1.2,
                  fontface = "bold"
        )
    } 
    if (nrow(evf_numeric_data) == 0) {
      warning("No environmental variables fit below the chosen significant level.")
    }
  } else {
    evf_numeric_model <- NULL
  }
  
  #################################### end of block ####################################
  
  #return plot or additional details
  if(output == "plot")
    return(plot)
  else if(output == "detailed")
    if (type == "nmds") {
      screeplot <- NULL
    } else {
      ### screeplot ###
      constrained_eig <- model$CCA$eig/model$tot.chi*100
      unconstrained_eig <- model$CA$eig/model$tot.chi*100
      if (length(constrained_eig) > 10) {
        constrained_eig <- constrained_eig[1:10]
      }
      if (length(unconstrained_eig) > 10) {
        unconstrained_eig <- unconstrained_eig[1:10]
      }
      eigenvalues <- c(constrained_eig, unconstrained_eig) #constrained combined with unconstrained
      screeplot <- ggplot(data.frame(axis = factor(names(eigenvalues), levels = names(eigenvalues)), eigenvalues), aes(x = axis, y = eigenvalues)) +
        geom_col() +
        geom_text(label = round(eigenvalues, 2), vjust = -1, size = 3)  +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab("Axis (max. 10 axes will be shown)") +
        ylab("Eigenvalue in percent of total inertia")
    }
    return(list(plot = plot,
                screeplot = screeplot,
                model = model,
                dsites = dsites,
                dspecies = dspecies,
                evf_factor_model = evf_factor_model,
                evf_numeric_model = evf_numeric_model
    )
    )
}
