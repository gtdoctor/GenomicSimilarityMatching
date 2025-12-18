rm(list=ls())
# Ancestry matching algorithm first developed by AP Levine.
# updated 2024-5 by GT Doctor

###### CHOICES #######


# PATHS # 
DATADIR=""
OUTDIR=""


# INPUTS #
PCA.EIGENVAL="10pca.svec.tsv" # list of eigenvals 1 per line.
PCA.EIGENVEC="10pca.tsv"  ##NB EXPECTS COLUMNS: FID, IID (FID CAN BE PLACEHOLDER) AND 10 PCS. NO HEADERS (OR HEADERLINE WITH #).
SAMPLEFILE="samplesincluded.tsv" ## EXPECTS COLUMNS FID, IID and then any number, must include a phenotpye column. (i.e. this determines which PC rows belong to cases and controls.  

SAMPLEFILE_IDCOL=2 #2 # equivalent to plink IID, can be changed
SAMPLEFILE_PHENOCOL=6 #6 # equivalent to plink PHENO/COHORT col. Expects controls to be enocded 1, cases to be encoded 2 
#"/mnt/nas1/projects/gabriel/INS/Oct24_gsaonly/ukbgsa270ins.forpca.fam""/mnt/nas1/projects/gabriel/INS/Dec2023/Merged.30may/mergedmegaUkbforpca.30may.fam" # plink.psam/orfam/" # plinkfile: headerless table with FID, IID columns (col1 is ignored entirely); or headed plink2 psam

# OUTPUTS #

FILEOUTSTEM= "distancematrix" # save the R workspace image. euc/eig/wtd will be appended
SCREENAME="Scree plot of principal components"    

# Parameters #
SMALLSAMPLE <- FALSE # subsets data to a small proportion of cases and controls for testing
euclidean <- TRUE # default. apply "FALSE" to calculate manhattan distances  
eigvalweight <- FALSE # should distances be weighted by the eigvalue of the PC.
nchunks = 30 # number of case chunks for parallel processing.

PCmax=10 ## distances can be calcuated on up to 10 PCs. 

######### CODE ############

GSM_A <- function(
  DATADIR, 
  OUTDIR, 
  PCA.eigenval, 
  PCA.eigenvec, 
  SAMPLEFILE,
  SAMPLEFILE_IDCOL,
  SAMPLEFILE_PHENOCOL,
  FILEOUTSTEM, 
  smallsample = FALSE,    # Optional argument with a default value
  euclidean = TRUE,     
  eigvalweight = FALSE,      
  nchunks,          # Default chunk size
  PCmax = 10,
  screename){
  
  
  
  library(ggplot2)
  library(parallel)
  
an = function(x){as.numeric(as.character(x))}
setwd(DATADIR)
FILEOUTSTEM <- paste0(FILEOUTSTEM,"_",PCmax,"pcs")

#Load PC data
d=read.table(PCA.eigenvec,header=F, comment.char = "#") # header to avoid problems if plink1.9 vs 2 is used. 
names(d)=c("FID","IID",paste("PC",1:10,sep="")) # take care that these match and that first line of the file isn't removed. 


#Load weightings
ev = an(read.table(PCA.eigenval, header=F)[,1])
ev.weight = ev/sum(ev) # weighting as proportion of top 10 eigenvals


# read case control. check here that the correct columns (e.g. XIID, PC1-10 are being selected)
cohort=read.table(SAMPLEFILE, header=T,comment.char = "#", )

# cases_i<-which(cohort[,SAMPLEFILE_PHENOCOL]==2)
# ctrls_i<-which(cohort[[SAMPLEFILE_PHENOCOL]]==1)
# cases <- d[cases_i,2:12]  # check
# controls <- d[ctrls_i,2:12]  # check
cases <- d[d$IID %in% cohort[[SAMPLEFILE_IDCOL]][cohort[[SAMPLEFILE_PHENOCOL]] == 2],2:12]
controls  <- d[d$IID %in% cohort[[SAMPLEFILE_IDCOL]][cohort[[SAMPLEFILE_PHENOCOL]] == 1],2:12]  # check that case control assignment is in col6 and IID in col 2



# small sample
if (smallsample == TRUE) {
  cases= cases[1:200,]
  controls=controls[1:2000,]
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

# Create case and control matrices
cases_matrix <- as.matrix(cases[, 2:(PCmax+1)])  # Exclude ID columns
controls_matrix <- as.matrix(controls[, 2:(PCmax+1)])  # Exclude ID columns
ncases=nrow(cases_matrix)
nctrls=nrow(controls_matrix)

# Initialize a matrix to store the distances
ccdistances <- matrix(nrow = ncases, ncol = nctrls)


# Create a sequence to split the data
chunk_size = ceiling(ncases / nchunks)
split_seq <- rep(1:nchunks, each = chunk_size, length.out = ncases) # vector of assignments for each line of case data -- 1repeated chunksize times, 2 repeated chunksize times etc.
casechunks <- split(cases_matrix, split_seq)

# Restructure each vector back into a matrix
casechunks <- lapply(casechunks, function(chunk) {
  matrix(chunk, ncol = ncol(cases_matrix), byrow = FALSE)
})

## Distance calculations

# Compute weighted Manhattan distances
if (euclidean == FALSE) {
# Create  (weighted) distances matrices
if (eigvalweight == TRUE) {
  filesave=paste0(FILEOUTSTEM,"_manwtd.Rdata")
  ev.weight = ev.weight[1:PCmax]

  print("Computing weighted Manhattan distances")

  distancefunction = function(casechunk, nctrls, controls_matrix, ev.weight) {
    # Initialize matrix to store results for this chunk with appropriate dimensions
    chunk_results = matrix(0, nrow(casechunk), nctrls)
    
    for (i in 1:nrow(casechunk)) {
      for (j in 1:nctrls) {
        # Calculate abs weighted difference for all PCs between case i and control j
        abs_diffs <- abs(casechunk[i, ] - controls_matrix[j, ]) * ev.weight
        chunk_results[i, j] <- sum(abs_diffs)
      }
    }
    chunk_results
  }
  results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix, ev.weight = ev.weight, mc.cores = detectCores())
  ccdistances = do.call(rbind, results)
  
  } # end Man weighted 

if (eigvalweight==FALSE){
  filesave=paste0(FILEOUTSTEM,"_manunwtd.Rdata")  
# Compute uneighted Manhattan distances
print("Computing unweighted Manhattan")


distancefunction = function(casechunk, nctrls, controls_matrix) {
  # Initialize matrix to store results for this chunk with appropriate dimensions
  chunk_results = matrix(0, nrow(casechunk), nctrls)
  
  for (i in 1:nrow(casechunk)) {
    for (j in 1:nctrls) {
      # Calculate abs  difference for all PCs between case i and control j
      abs_diffs <- abs(casechunk[i, ] - controls_matrix[j, ]) 
      chunk_results[i, j] <- sum(abs_diffs)
    }
  }
  chunk_results
}
results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix,  mc.cores = detectCores())
ccdistances = do.call(rbind, results)


} # End man weighted
} # end euclidean == "FALSE"

if (euclidean == TRUE) {  

if (eigvalweight == FALSE) {
  print("Computing unweighted Euclidean distances")
    filesave=paste0(FILEOUTSTEM,"_eucunwtd.Rdata")
    

  distancefunction = function(casechunk, nctrls, controls_matrix) {
      # Initialize matrix to store results for this chunk with appropriate dimensions
      chunk_results = matrix(0, nrow(casechunk), nctrls)
      for (i in 1:nrow(casechunk)) {
        for (j in 1:nctrls) {
          # Calculate squared difference for all PCs between case i and control j
          squared_distances <- (casechunk[i, ] - controls_matrix[j, ])^2 # squared_distances is a vector as is ev weight; multiples matching elements
          # Sum and take the square root to get weighted Euclidean distance
          chunk_results[i, j] <- sqrt(sum(squared_distances))
        }
      }
      chunk_results
    }
    results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix, mc.cores = detectCores())
    ccdistances = do.call(rbind, results)
    
    } # end euc unweighted

if (eigvalweight ==TRUE) {
  print("Computing weighted Euclidean distances")
  filesave=paste0(FILEOUTSTEM,"_eucwtd.Rdata")
  ev.weight = ev.weight[1:PCmax]
  
    distancefunction = function(casechunk, nctrls, controls_matrix,ev.weight) {
    # Initialize matrix to store results for this chunk with appropriate dimensions
    chunk_results = matrix(0, nrow(casechunk), nctrls)
    for (i in 1:nrow(casechunk)) {
      for (j in 1:nctrls) {
        # Calculate squared weighted difference for all PCs between case i and control j
        squared_distances <- (casechunk[i, ] - controls_matrix[j, ])^2 * ev.weight # squared_distances is a vector as is ev weight; multiples matching elements
        # Sum and take the square root to get weighted Euclidean distance
        chunk_results[i, j] <- sqrt(sum(squared_distances))
      }
    }
    chunk_results
  }
  results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix, ev.weight = ev.weight, mc.cores = detectCores())
  ccdistances = do.call(rbind, results)
}  
}  # end euclidean == "TRUE"

# Clean up the cluster
if (exists("cl")) {
  parallel::stopCluster(cl)
}

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
               nchunks = nchunks ,
      screename =SCREENAME,
               PCmax = PCmax)
