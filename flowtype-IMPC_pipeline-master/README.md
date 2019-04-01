Note: code marked with larger numbers are dependant on code marked with smaller numbers

* indicates no other code is dependant on this code



00_collateFT.R :: collates flowtype files into a matrix; imports and standardizes sampleMeta matrix

01a_normalize.R :: normalizes flowtype count matrix into normalized count matrix

01_childparent_ind.R :: organizes phenoMeta (organizes cell hierarchy meta)

02_count_stats.R* :: plots cell count

02_plot.R :: plots cell count over time (includes changepoint analysis), plots dimensionality reduction

02_pvalue_single.R :: (dependant on 02_plot) calculates pvalues & log ratios for single files

03_childparent_matrix.R :: creates feature matrices
04_checkmatrix.R :: overviews matrices made so far

04_dist.R :: calculates standard distances for different features

04_dist_w.R :: calculates weighted distances for different features

04_dist_lin.R :: linear combination of 04_dist or 04_dist_w distances

05_dist_plot.R :: plot dimensionality reduction of distances

05_dist_score.R :: calculates NCA score for distances

