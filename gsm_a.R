rm(list=ls())
# Ancestry matching algorithm first developed by AP Levine.
# updated 2024-6 by GT Doctor

###### CHOICES #######


# PATHS #
DATADIR=""  # directory holding PCA.EIGENVAL, PCA.EIGENVEC and SAMPLEFILE
OUTDIR=""   # directory the distance matrix is written to (created if absent)


# INPUTS #
PCA.EIGENVAL="10pca.svec.tsv" # eigenvalues (lambda, proportional to squared singular values), 1 per line.
PCA.EIGENVEC="10pca.tsv"  ##NB EXPECTS COLUMNS: FID, IID (FID CAN BE PLACEHOLDER) AND 10 PCS.
SAMPLEFILE="samplesincluded.tsv" ## EXPECTS COLUMNS FID, IID and then any number must include a phenotpye column.

SAMPLEFILE_IDCOL=2 #2 # equivalent to plink IID, can be changed
SAMPLEFILE_PHENOCOL=6 #6 # equivalent to plink PHENO/COHORT col. Expects controls to be enocded 1, cases to be encoded 2
# Accepts a headerless plink .fam (FID, IID in cols 1-2) or a headed plink2 .psam.
# Extra trailing columns are allowed: only the indexed columns are read.

# OUTPUTS #

FILEOUTSTEM= "distancematrix" # save the R workspace image. euc/man + wtd/unwtd + scaled will be appended
SCREENAME="Scree plot of principal components"

# Parameters #
PCmax=10 ## distances can be calcuated on up to 10 PCs.
euclidean <- TRUE # default. apply "FALSE" to calculate manhattan distances
eigvalweight <- FALSE # should distances be weighted by the eigvalue of the PC.
scale <- FALSE
SMALLSAMPLE <- FALSE # subsets data to a small proportion of cases and controls for testing
nchunks = 40 # number of case chunks for parallel processing. Reduce RAM usage if large dataset. 
MCCORES = 8 # parallel workers. Each holds a (chunk cases x all controls) block.


######### CODE ############

GSM_A <- function(
  DATADIR,
  OUTDIR,
  PCA.eigenval,
  PCA.eigenvec,
  SAMPLEFILE,
  SAMPLEFILE_IDCOL= 2 ,
  SAMPLEFILE_PHENOCOL= 6,
  FILEOUTSTEM,
  smallsample = FALSE,
  euclidean = TRUE,
  eigvalweight = FALSE,
  scale=FALSE,
  nchunks = 10,          # Default chunk size
  mccores = 4,
  PCmax = PCmax,
  screename){



  library(ggplot2)
  library(parallel)
  library(data.table)

setwd(DATADIR)
dir.create(OUTDIR,recursive = TRUE, showWarnings = FALSE)
FILEOUTSTEM <- paste0(FILEOUTSTEM,"_",PCmax,"pcs")


#Load PC data
d=fread(PCA.eigenvec) # header to avoid problems if plink1.9 vs 2 is used.
names(d)=c("FID","IID",paste0("PC",1:10)) # take care that these match and that first line of the file isn't removed.

#Load weightings
ev = fread(PCA.eigenval)[,1] # fread to read in.
ev <- as.numeric(ev[[1]])   # force vector

ev.weight = ev/sum(ev) # weighting as proportion of top 10 eigenvals


# read case control. check here that the correct columns (e.g. XIID, PC1-10 are being selected)
cohort=fread(SAMPLEFILE)

## data structure checks
print("Head of PCA file:")
print(head(d,2))

print("Head of SAMPLEFILE:")
print(head(cohort,2))


print(paste0("Sense check: nrow of PC scores loaded: ", nrow(d)))
print(paste0("Sense check: nrow of sample IDs loaded: ", nrow(cohort)))


id_pca<- d$IID
id <- cohort[[SAMPLEFILE_IDCOL]]

## no NA / missing
if (any(is.na(id_pca)) || any(id_pca == "")) {
  stop("NA or missing values in id_pca", call. = FALSE)
}
if (any(is.na(id)) || any(id == "")) {
  stop("NA or missing values in cohort IDs", call. = FALSE)
}

## no duplicates (within each)
if (any(duplicated(id_pca))) {
  stop("Duplicated IDs in id_pca", call. = FALSE)
}
if (any(duplicated(id))) {
  stop("Duplicated IDs in cohort IDs", call. = FALSE)
}

## identical order and values
if (!identical(id_pca, id)) {
  stop("ID vectors differ in values and/or order between PCA and cohort files",
       call. = FALSE)
}


## Per-PC coordinate multiplier. Single place where whitening and eigenvalue
## weighting are applied, so that Euclidean and Manhattan rescale the
## coordinates identically.
## SDs are taken over cases and controls pooled.
pc_cols <- paste0("PC", 1:PCmax)
mult <- rep(1, PCmax)

if (scale==TRUE){
mult <- mult / vapply(pc_cols, function(cc) sd(d[[cc]]), numeric(1))
FILEOUTSTEM=paste0(FILEOUTSTEM,"_scaled")
}

if (eigvalweight==TRUE){
mult <- mult * sqrt(ev.weight[1:PCmax])
}


case_ids <- cohort[[SAMPLEFILE_IDCOL]][cohort[[SAMPLEFILE_PHENOCOL]] == 2]
idx <- match(case_ids, d$IID)
cases <- d[idx, 2:12]

control_ids <- cohort[[SAMPLEFILE_IDCOL]][cohort[[SAMPLEFILE_PHENOCOL]] == 1]
idx <- match(control_ids, d$IID)
controls <- d[idx, 2:12]



# small sample
if (smallsample == TRUE) {
  cases= cases[1:500,]
  controlnames=c(paste0("C",1:nrow(controls)))
  controls$IID=controlnames
}


## Scree Plot - to decide on PCs needed
scree <- data.frame(PC = 1:length(ev.weight), VarianceExplained = ev.weight * 100)
#
# # Plotting the scree plot
# ggplot(scree, aes(x = factor(PC), y = VarianceExplained)) +
#   geom_bar(stat = "identity") +
#   xlab("Principal Component") +
#   ylab("Variance Explained (%)") +
#   ggtitle("Scree Plot") +
#   ylim(0, 100)

# Create case and control matrices, rescaled by the coordinate multiplier
cases_matrix <- as.matrix(cases[, 2:(PCmax+1)])  # Exclude ID columns
controls_matrix <- as.matrix(controls[, 2:(PCmax+1)])  # Exclude ID columns
cases_matrix <- sweep(cases_matrix, 2, mult, "*")
controls_matrix <- sweep(controls_matrix, 2, mult, "*")
ncases=nrow(cases_matrix)
nctrls=nrow(controls_matrix)


# Create a sequence to split the data

chunk_size = ceiling(ncases / nchunks)
split_seq <- rep(1:nchunks, each = chunk_size, length.out = ncases) # vector of assignments for each line of case data -- 1repeated chunksize times, 2 repeated chunksize times etc.
# indices for each case row
idx_chunks <- split(seq_len(ncases), split_seq)
casechunks <- lapply(idx_chunks, function(idx) cases_matrix[idx, , drop = FALSE])



## Distance kernels. Weighting/whitening is already in the coordinates, so
## each metric has exactly one implementation.

dist_euc <- function(X, Y) {
  #d^2= sum_k(x_k-y_k)^2)
  # =  sum_k(x^2) + sum_k(y^2) - 2 x.y
  x2 <- rowSums(X^2)
  y2 <- rowSums(Y^2)
  d2 <- outer(x2, y2, "+") - 2 * (X %*% t(Y))
  sqrt(pmax(d2, 0))   # guard tiny negative from FP error
}

dist_man <- function(X, Y) {
  #  |x_k - y_k| one PC at a time; outer gives an m x n block
  out <- matrix(0, nrow(X), nrow(Y))
  for (k in seq_len(ncol(X))) {
    out <- out + abs(outer(X[, k], Y[, k], "-"))
  }
  out
}

filesave <- paste0(FILEOUTSTEM, "_",
                   if (euclidean==TRUE) "euc" else "man",
                   if (eigvalweight==TRUE) "wtd" else "unwtd",
                   ".Rdata")
kern <- if (euclidean==TRUE) dist_euc else dist_man

print(sprintf("Computing %s distances: %s cases x %s controls in %s chunks on %s cores",
              if (euclidean==TRUE) "Euclidean" else "Manhattan",
              ncases, nctrls, length(casechunks), mccores))

results <- parallel::mclapply(
  casechunks,
  kern,
  Y = controls_matrix,
  mc.cores = min(mccores, parallel::detectCores())
)
ccdistances = do.call(rbind, results)


setwd(OUTDIR)
#ggsave(plot=last_plot(), filename = paste0(screename, ".jpg"), device = "jpeg")
print(paste0("Saving ",filesave))
save(ccdistances, file=filesave)
}

GSM_A(
  DATADIR=DATADIR,
               OUTDIR=OUTDIR,
               PCA.eigenval=PCA.EIGENVAL,
               PCA.eigenvec= PCA.EIGENVEC,
               SAMPLEFILE=SAMPLEFILE,
  SAMPLEFILE_IDCOL = SAMPLEFILE_IDCOL,
  SAMPLEFILE_PHENOCOL = SAMPLEFILE_PHENOCOL,
               FILEOUTSTEM = FILEOUTSTEM,
               smallsample = SMALLSAMPLE ,
               euclidean = euclidean ,
               eigvalweight = eigvalweight ,
               scale = scale,
               nchunks = nchunks ,
               mccores = MCCORES,
              screename =SCREENAME,
               PCmax = PCmax)
