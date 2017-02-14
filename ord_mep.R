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
    ncomp = 3, 
    color = NULL, #Color the dots by some variable in the metadata
    type = "PCA", #ordination type; PCA, RDA, NMDS, MMDS(or PCOA), CA, CCA
    metric = "sqrt", #distance metric to use to calculate distance matrix for input to the ordination functions
    constrain = NULL, #variable(s) for constrained analyses; RDA and CCA
    output = "plot", #"plot" or "detailed"; output as list with additional information(model, scores, inputmatrix etc) or just the plot
    ... #Pass additional arguments to the vegan ordination functions, fx rda(...), cca(...), metaMDS(...), see vegan help
) {
    #Check the data
    datalist <- amp_rename(data = datalist, tax.empty = "best")
    
    #Calculate distance/variance/covariance/dissimilarity/similarity matrix from vegdist()
    if(metric == "sqrt") {
        inputmatrix <- t(sqrt(datalist$abund))
    } else if(metric == "none") {
        inputmatrix <- t(datalist$abund)
    } else if(any(metric == c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"))) {
        inputmatrix <- vegdist(t(datalist$abund), method = metric)
    }
    
    #Generate data depending on the chosen ordination type
    if(type == "PCA") {
        #make the model
        model <- rda(inputmatrix, ...)
        
        #axis (and data column) names
        x_axis_name <- "PC1"
        y_axis_name <- "PC2"
        
        #Calculate the percentage of eigenvalues explained by the axes
        totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)
        
        #calculate site scores and combine with metadata
        scores <- scores(model, choices = 1:ncomp)$sites
        d <- cbind.data.frame(datalist$metadata, scores)
    } else if(type == "RDA") {
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
        
        #calculate site scores and combine with metadata
        scores <- scores(model, choices = 1:ncomp)$sites
        d <- cbind.data.frame(datalist$metadata, scores)
    } else if(type == "NMDS") {
        #make the model
        model <- metaMDS(inputmatrix, trace = FALSE, distance = metric, ...)
        
        #axis (and data column) names
        x_axis_name <- "NMDS1"
        y_axis_name <- "NMDS2"
        
        #calculate site scores and combine with metadata
        scores <- scores(model, choices = 1:ncomp)
        d <- cbind.data.frame(datalist$metadata, scores)
    } else if(type == "MMDS" | type == "PCOA") {
        #make the model
        model <- cmdscale(inputmatrix, eig = TRUE, ...)
        
        #axis (and data column) names
        x_axis_name <- "MDS1"
        y_axis_name <- "MDS2"
        
        #calculate site scores and combine with metadata
        scores <- scores(model, choices = 1:ncomp)
        colnames(scores) <- c(x_axis_name, y_axis_name)
        d <- cbind.data.frame(datalist$metadata, scores)
    } else if(type == "CA") {
        #make the model
        model <- cca(inputmatrix, ...)
        
        #axis (and data column) names
        x_axis_name <- "CA1"
        y_axis_name <- "CA2"
        
        #Calculate the percentage of eigenvalues explained by the axes
        totalvar <- round(model$CA$eig/model$CA$tot.chi * 100, 1)
        
        #calculate site scores and combine with metadata
        scores <- scores(model, choices = 1:ncomp)$sites
        d <- cbind.data.frame(datalist$metadata, scores)
    } else if(type == "CCA") {
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
        
        #calculate site scores and combine with metadata
        scores <- scores(model, choices = 1:ncomp)$sites
        d <- cbind.data.frame(datalist$metadata, scores)
    }
    # Generate a nice ggplot with the coordinates from scores
    plot <- ggplot(d,
                   aes_string(x = x_axis_name, y = y_axis_name, color = color)
    ) +
        geom_point(size = 2) + 
        theme_minimal() +
        theme(axis.line = element_line(colour = "black", size = 0.5))
    if(type == "PCA" | type == "RDA" | type == "CA" | type == "CCA") {
        plot <- plot +
            xlab(paste(x_axis_name, " [", totalvar[x_axis_name], "%]", sep = "")) + 
            ylab(paste(y_axis_name, " [", totalvar[y_axis_name], "%]", sep = ""))
    }
    if(output == "plot")
        return(plot)
    else if(output == "detailed")
        return(list(plot = plot, model = model, scores = scores))
}