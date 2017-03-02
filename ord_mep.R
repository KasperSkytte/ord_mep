#A wrapper around the vegan package to make beautiful ggplot2 ordination plots 
#for use in microbial ecology profiling. The input data list must be loaded with
#amp_load() from the ampvis package, but any OTU table-like matrix (OTU's in rows,
#sampleID's in columns and abundances in the corresponding cells) alongside metadata 
#for the same samples can be used. Simply choose the desired distance metric and 
#ordination method and a plot is returned.
#Required packages:
#dplyr, ampvis, vegan, 
ord_mep <- function(
    datalist, #with 3 elements loaded with amp_load: abund, metadata and tax, first two required
    color_by = NULL, #Color points by a variable in the metadata
    type = "PCA", #Ordination type; PCA, RDA, NMDS, MMDS(or PCOA), CA, CCA, DCA
    metric = NULL, #Distance metric to use to calculate distance matrix for input to the ordination functions, any of the following:
    #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "jsd" (Jensen-Shannon Divergence),
    #"altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"
    #or simply "none" or "sqrt"
    constrain = NULL, #Variable(s) in the metadata for constrained analyses; RDA and CCA
    colorframe = TRUE, #Frame the points with a polygon colored by the color_by argument
    label_by = NULL, #Display text labels frin a variable in the metadata, can also be used to plot environmental variables
    plot_species_points = FALSE,
    output = "plot", #"plot" or "detailed"; output as list with additional information(model, scores, inputmatrix etc) or just the plot
    ... #Pass additional arguments to the vegan ordination functions, fx rda(...), cca(...), metaMDS(...), see vegan help
) {
    if(is.null(metric)) {
        warning("No metric selected. If this is not deliberate, please provide one with the argument: metric")
        metric == "none"
    }
    #Check the data
    datalist <- amp_rename(data = datalist, tax.empty = "best")
    
    #to fix user argument characters, so fx PCoA/PCOA/pcoa are all valid
    type <- tolower(type)
    metric <- tolower(metric)
    output <- tolower(output)
    
    #Calculate distance/variance/covariance/dissimilarity/similarity matrix from vegdist()
    if(metric == "sqrt") {
        inputmatrix <- t(sqrt(datalist$abund))
    } else if(metric == "none") {
        inputmatrix <- t(datalist$abund)
    } else if(metric == "jsd") {
        #This is based on http://enterotype.embl.de/enterotypes.html
        #Abundances of 0 will be set to the pseudocount value to avoid 0-value denominators
        dist.JSD <- function(inMatrix, pseudocount=0.000001) {
            KLD <- function(x,y) sum(x *log(x/y))
            JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
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
        x_axis_name <- "PC1"
        y_axis_name <- "PC2"
        
        #Calculate the percentage of eigenvalues explained by the axes
        totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)
        
        #Calculate species- and site scores
        sitescores <- scores(model, display = "sites")
        speciesscores <- scores(model, display = "species")
    } else if(type == "rda") {
        if(is.null(constrain)) 
            stop("Argument constrain must be provided when performing constrained/canonical analysis.")
        #make the model
        codestring <- paste0("rda(inputmatrix~", paste(constrain, collapse = "+"), ", datalist$metadata, ...)") #function arguments written as in rda(x ~ y + z) hard to pass to rda(), now user just provides a vector
        model <-  eval(parse(text = codestring))
        
        #axes depends on the results
        x_axis_name <- "RDA1"
        y_axis_name <- "PC1"
        if (model$CCA$rank > 1) {
            y_axis_name <- "RDA2"
        }
        
        #Calculate the percentage of eigenvalues explained by the axes
        totalvar <- c(round(model$CA$eig/model$tot.chi * 100, 1),
                      round(model$CCA$eig/model$CA$tot.chi * 100, 1)
        )
        
        #Calculate species- and site scores
        sitescores <- scores(model, display = "sites")
        speciesscores <- scores(model, display = "species")
    } else if(type == "nmds") {
        #make the model
        model <- metaMDS(inputmatrix, trace = FALSE, distance = metric, ...)
        
        #axis (and data column) names
        x_axis_name <- "NMDS1"
        y_axis_name <- "NMDS2"
        
        #Calculate species- and site scores
        #Speciesscores may not be available with MDS
        sitescores <- scores(model, display = "sites")
        if(ncol(scores(model)) == 2) {
            speciesscores <- warning("Speciesscores are unfortunately not available")
        } else {
            speciesscores <- scores(model, display = "species")
        }
    } else if(type == "mmds" | type == "pcoa") {
        #make the model
        model <- cmdscale(inputmatrix, eig = TRUE, ...)
        
        #axis (and data column) names
        x_axis_name <- "PCoA1"
        y_axis_name <- "PCoA2"
        
        #Calculate species- and site scores
        #Speciesscores may not be available with MDS
        sitescores <- scores(model, display = "sites")
        colnames(sitescores) <- c(x_axis_name, y_axis_name)
        if(ncol(scores(model)) == 2) {
            speciesscores <- warning("Speciesscores are unfortunately not available")
        } else {
            speciesscores <- scores(model, display = "species")
        }
    } else if(type == "ca") {
        #make the model
        model <- cca(inputmatrix, ...)
        
        #axis (and data column) names
        x_axis_name <- "CA1"
        y_axis_name <- "CA2"
        
        #Calculate the percentage of eigenvalues explained by the axes
        totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)
        
        #Calculate species- and site scores
        sitescores <- scores(model, display = "sites")
        speciesscores <- scores(model, display = "species")
    } else if(type == "cca") {
        if(is.null(constrain)) 
            stop("Argument constrain must be provided when performing constrained/canonical analysis.")
        #make the model
        codestring <- paste0("cca(inputmatrix~", paste(constrain, collapse = "+"), ", datalist$metadata, ...)") #function arguments written as in rda(x ~ y + z) hard to pass to rda(), now user just provides a vector
        model <-  eval(parse(text = codestring))
        
        #axes depends on the results
        x_axis_name <- "CCA1"
        y_axis_name <- "CA1"
        if (model$CCA$rank > 1) {
            y_axis_name <- "CCA2"
        }
        
        #Calculate the percentage of eigenvalues explained by the axes
        totalvar <- c(round(model$CA$eig/model$tot.chi * 100, 1),
                      round(model$CCA$eig/model$CA$tot.chi * 100, 1)
        )
        
        #Calculate species- and site scores
        sitescores <- scores(model, display = "sites")
        speciesscores <- scores(model, display = "species")
    } else if(type == "dca") {
        #make the model
        model <- decorana(inputmatrix, ...)
        
        #axis (and data column) names
        x_axis_name <- "DCA1"
        y_axis_name <- "DCA2"
        
        #Calculate the percentage of eigenvalues explained by the axes
        #totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)
        
        #Calculate species- and site scores
        sitescores <- scores(model, display = "sites")
        speciesscores <- scores(model, display = "species")
    }
    #################################### end of block ####################################
    
    #Make data frames for ggplot
    if(length(speciesscores) > 1) {
        dspecies <- cbind.data.frame(speciesscores)
    }
    dsites <- cbind.data.frame(datalist$metadata, sitescores)
    
    #Generate a nice ggplot with the coordinates from scores
    plot <- ggplot(dsites,
                   aes_string(x = x_axis_name, y = y_axis_name, color = color_by)
    ) +
        geom_point(size = 2) + 
        theme_minimal() +
        theme(axis.line = element_line(colour = "black", size = 0.5))
    
    #only some ordination methods can be displayed with % on axes
    if(type == "pca" | type == "rda" | type == "ca" | type == "cca") {
        plot <- plot +
            xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) + 
            ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "%]", sep = ""))
    }
    
    #Generate a color frame around the chosen color group
    if(colorframe == TRUE) {
        splitData <- split(dsites, dsites[, color_by]) %>% 
            lapply(function(df) {
                df[chull(df[, x_axis_name], df[, y_axis_name]), ]
            })
        hulls <- do.call(rbind, splitData)
        plot <- plot + geom_polygon(data = hulls, aes_string(fill = color_by, group = color_by), alpha = 0.2)
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
        plot <- plot + geom_text(data = temp3,
                                 aes_string(x = "cx", 
                                            y = "cy",
                                            label = label_by),
                                 size = 3,
                                 color = "black",
                                 fontface = 2
                                 )
    }
    
    #Plot species points
    if (plot_species_points == TRUE) {
        plot <- plot + geom_point(data = dspecies, color = "grey", 
                            shape = 20, size = 2)
    }
    #################################### end of block ####################################
    
    #return plot or additional details
    if(output == "plot")
        return(plot)
    else if(output == "detailed")
        return(list(plot = plot, model = model, sitescores = sitescores, speciesscores = speciesscores))
}