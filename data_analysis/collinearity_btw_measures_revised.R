#' This calculates the collinearity between the different QC measures
#' 
#' We first load the abide and corr dataset
#' then we compute the correlation between the measures within each site
#' we average across sites after doing an r->t->z (this is the difference with the other script)
#' and finally we plot that difference
#' 
#+ setup
library(corrplot)
library(RColorBrewer)
library(plyr)
library(fdrtool)

#' # ABIDE
#' 
#' ## Setup
#+ abide-setup
# Read
qc_anat <- read.csv("abide_anat_rate.csv", row.names=1)
qc_func <- read.csv("abide_func.csv", row.names=1)

#names(qc_anat)[1:4] = c('subject','session','scan','site')
#names(qc_func)[1:4] = c('subject','session','scan','site')

# Fetch relevant columns
anat_cols <- names(qc_anat)[5:10]#c("cnr",  "efc",  "fber", "fwhm", "qi1",  "snr")
func_cols <- c('EFC','FBER','FWHM','SNR','Fr.Outliers','GCOR','Quality','RMSD','DVARS' )#names(qc_func)[c(5:11,13:22)] #c("efc","fber","fwhm","snr","dvars","mean_fd","quality",'zgcor')#names(qc_func)[5:13]

#' ## Correlate
#+ abide-correlate
# Anatomical
all_zmats <- daply(qc_anat, .(site), function(x) {
  cmat <- cor(x[,anat_cols], use="pairwise.complete.obs")
  cmat[cmat>0.9999999] <- 0.9999999 # clamp
  zmat <- atanh(cmat) * sqrt(nrow(x) - 3)
  zmat
})
ns <- daply(qc_anat, .(site), nrow)
mean_df <- mean(sqrt(ns-3))
abide_anat <- apply(all_zmats, c(2,3), function(x) tanh(mean(x)/mean_df))
abide_anat_p <- pcor0(abs(abide_anat), mean(ns), lower.tail=F)
write.csv(abide_anat, file='abide_anat_cor.csv')
write.csv(abide_anat_p, file='abide_anat_cor_pval.csv')
# Functional
all_zmats <- daply(qc_func, .(site), function(x) {
  cmat <- cor(x[,func_cols], use="pairwise.complete.obs")
  cmat[cmat>0.9999999] <- 0.9999999 # clamp
  zmat <- atanh(cmat) * sqrt(nrow(x) - 3)
  zmat
})
ns <- daply(qc_func, .(site), nrow)
mean_df <- mean(sqrt(ns-3))
abide_func <- apply(all_zmats, c(2,3), function(x) tanh(mean(x)/mean_df))
abide_func_p <- pcor0(abs(abide_func), mean(ns), lower.tail=F)
write.csv(abide_func, file='abide_func_cor.csv')
write.csv(abide_func_p, file='abide_func_cor_pval.csv')

# in case you want the original way
#abide_anat_orig <- cor(qc_anat[,anat_cols], use="pairwise.complete.obs")
#abide_func_orig <- cor(qc_func[,func_cols], use="pairwise.complete.obs")

#' # CoRR
#' 
#' ## Setup
#+ corr-setup
# Read
qc_anat <- read.csv("corr_anat.csv", row.names=1)
qc_func <- read.csv("corr_func.csv", row.names=1)

#names(qc_anat)[1:4] = c('subject','session','scan','site')
#names(qc_func)[1:4] = c('subject','session','scan','site') 
# Fetch relevant columns
anat_cols <- names(qc_anat)[5:11]#c("cnr",  "efc",  "fber", "fwhm", "qi1",  "snr")# names(qc_anat)[6:11]
func_cols <- c('EFC','FBER','FWHM','SNR','Fr.Outliers','GCOR','Quality','RMSD','DVARS' )#names(qc_func)[c(5:13,15:30 )]#c("efc","fber","fwhm","snr","dvars","mean_fd","quality",'gcor')#names(qc_func)[6:13] #c(4:9,12)]
# Only include session 1 and scan 1
qc_anat <- subset(qc_anat, session == 1)
qc_func <- subset(qc_func, session == 1 & series == 1)
# Fix these weird labels

#' ## Correlate
#+ corr-correlate
# anatomical
all_zmats <- daply(qc_anat, .(site), function(x) {
  cmat <- cor(x[,anat_cols], use="pairwise.complete.obs")
  cmat[cmat>0.9999999] <- 0.9999999 # clamp
  zmat <- atanh(cmat) * sqrt(nrow(x) - 3)
  zmat
})
ns <- daply(qc_anat, .(site), nrow)
inds <- ns>10
mean_df <- mean(sqrt(ns[inds]-3))
corr_anat <- apply(all_zmats[inds,,], c(2,3), function(x) tanh(mean(x)/mean_df))
corr_anat_p <- pcor0(abs(corr_anat), mean(ns), lower.tail=F)
write.csv(corr_anat, file='corr_anat_cor.csv')
write.csv(corr_anat_p, file='corr_anat_cor_pval.csv')
# functional
all_zmats <- daply(qc_func, .(site), function(x) {
  cmat <- cor(x[,func_cols], use="pairwise.complete.obs")
  cmat[cmat>0.9999999] <- 0.9999999 # clamp
  zmat <- atanh(cmat) * sqrt(nrow(x) - 3)
  zmat #zmat
})
ns <- daply(qc_func, .(site), nrow)
inds <- ns>10
mean_df <- mean(sqrt(ns[inds]-3))
corr_func <- apply(all_zmats[inds,,], c(2,3), function(x) tanh(mean(x)/mean_df))
corr_func_p <- pcor0(abs(corr_func), mean(ns), lower.tail=F)

write.csv(corr_func, file='corr_func_cor.csv')
write.csv(corr_func_p, file='corr_func_cor_pval.csv')

# old way
#corr_anat_orig <- cor(qc_anat[,anat_cols], use="pairwise.complete.obs")
#corr_func_orig <- cor(qc_func[,func_cols], use="pairwise.complete.obs")


#' # Visualize
#' 
#' We look at the abide and corr results together.
#' 
#' ## Fixes
#' Since CoRR is missing outliers, we remove that from the abide results
#+ viz-fix-rmcol
#rmcol <- which(colnames(abide_func) == "outlier")
#abide_func <- abide_func[-rmcol,-rmcol]
#abide_func_p <- abide_func_p[-rmcol,-rmcol]

#' We also want to ensure the column names are in the same order.
#' The anatomical is actually fine so only do func cols.
#' We match in reference to the abide functional colnames
#+ viz-fix-ordercol
# anat
all(colnames(abide_anat) == colnames(corr_anat)[-2])
# func
#ref <- colnames(abide_func)
#trg <- colnames(corr_func)
#inds <- sapply(ref, function(x) which(trg == x))
#corr_func0 <- corr_func # just in case
#corr_func <- corr_func0[inds,inds]

#' Now we can combine the two together for one big happy family
#'+ viz-fix-combine
# anat
anat <- abide_anat
anat[upper.tri(anat)] <- corr_anat[upper.tri(anat)]
anat_p <- abide_anat_p
anat_p[upper.tri(anat_p)] <- corr_anat_p[upper.tri(anat_p)]
# func (remove the num_fd column since redundant)
func <- abide_func
func[upper.tri(func)] <- corr_func[upper.tri(func)]
#rmind <- which(colnames(func) %in% c("num_fd", "perc_fd"))#== "num_fd")
#func <- func[-rmind,-rmind]
func_p <- abide_func_p
func_p[upper.tri(func_p)] <- corr_func_p[upper.tri(func_p)]
#rmind <- which(colnames(func_p) == "num_fd")
#func_p <- func_p[-rmind,-rmind]

#+ relabel
## anat
#old_labels <- c("cnr", "efc", "fber", "fwhm", "qi1", "snr")
#new_labels <- c("CNR", "EFC", "FBER", "FWHM", "Qi1", "SNR")
#inds <- sapply(old_labels, function(x) which(colnames(anat)==x))
#colnames(anat) <- new_labels[inds]
#rownames(anat) <- new_labels[inds]
## func
#old_labels <- c( "efc","fber","fwhm","snr","dvars","mean_fd","quality",'zgcor') #"perc_fd",
#new_labels <- c("EFC", "FBER", "FWHM", "SNR", "DVARS", "Mean FD",  "Quality",'Global corr') #"Percent FD",
#inds <- sapply(old_labels, function(x) which(colnames(func)==x))
#colnames(func) <- new_labels[inds]
#rownames(func) <- new_labels[inds]

#' ## Plot
#' Generate the corrplot
#+ viz-plot
cols <- rev(brewer.pal(10, "RdBu"))
# the anatomical
png("bysite_anat_corrplot.png");
corrplot(anat, method="color", diag=F, outline=F, col=cols, cl.length=length(cols)+1)
dev.off()

#savePlot("bysite_anat_corrplot.png")
# the functional
png("bysite_func_corrplot.png");
corrplot(func, method="color", diag=F, outline=F, col=cols, cl.length=length(cols)+1)
dev.off()
#savePlot("bysite_func_corrplot.png")

#+ viz-plot-sig
cols <- rev(brewer.pal(10, "RdBu"))
# the anatomical
png("figure/bysite_sig_anat_corrplot.png")
corrplot(anat, method="color", diag=F, outline=F,  col=cols, cl.length=length(cols)+1, 
         tl.cex=2, tl.col="black", p.mat=anat_p)
dev.off()
#quartz.save("bysite_sig_anat_corrplot.png", width=5, height=5)
# the functional
png("figure/bysite_sig_func_corrplot.png")
corrplot(func, method="color", diag=F, outline=F, col=cols, cl.length=length(cols)+1, 
          tl.col="black", p.mat=func_p)
#quartz.save("bysite_sig_func_corrplot.png", width=4, height=4)
dev.off()

#+ viz-plot-thr
cols <- rev(brewer.pal(10, "RdBu"))
# the anatomical
png("figure/bysite_thr_anat_corrplot.png")
corrplot(anat, method="color", diag=F, outline=F,  col=cols, cl.length=length(cols)+1, p.mat=anat_p, insig='blank')
dev.off()
#quartz.save("bysite_thr_anat_corrplot.png", width=5, height=5)
# the functional
png("figure/bysite_thr_func_corrplot.png")
corrplot(func, method="color", diag=F, outline=F, col=cols, cl.length=length(cols)+1, p.mat=func_p, insig='blank')
dev.off()
#quartz.save("bysite_thr_func_corrplot.png", width=4, height=4)

