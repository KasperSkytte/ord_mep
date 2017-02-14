# ord_mep
A wrapper around the vegan package to make beautiful ggplot2 ordination plots for use in microbial ecology profiling.
The input data list must be loaded with amp_load() from the ampvis package, but any OTU table-like matrix
(OTU's in rows, sampleID's in columns and abundances in the corresponding cells) alongside metadata
for the same samples can be used. Simply choose the desired distance metric and ordination
method and a plot is returned.
