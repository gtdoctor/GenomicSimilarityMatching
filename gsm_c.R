rm(list=ls())
# INPUTS #
args<-commandArgs(trailingOnly=TRUE)
PCA.eigenvec=args[1] # PCA output, headed or headerless, with FID, IID, PCs
PSAM <- args[2]      # sample file, as used in parts A and B
KEPT <- args[3]      # full path to the kept.*.tsv written by part B

# CHOICES # 
PCMAX = 10 # number of PCs in the dataset. 
NORMPLOT <- "y" # choose y or n, to plot normalised data (which puts all PCs onto same scale)

##### function #####
GSM_C_plotting <- function(PCA.eigenvec = PCA.eigenvec, PSAM = PSAM, KEPT = KEPT, PCMAX = PCMAX, NORMPLOT = NORMPLOT){
  load_or_install <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

load_or_install(c("ggplot2","grid","gridExtra","data.table"))

# arrangeGrob measures grobs against the current device; with none open R
# starts the default png device and leaves Rplot00N.png in the output directory.
pdf(NULL)
on.exit(dev.off(), add = TRUE)


if (!exists("cohort")) {
# PSAM may carry extra trailing columns (e.g. genotyping chip) — name only the first 6.
cohort = fread(PSAM, header = FALSE)
setnames(cohort, 1:6, c("FID", "IID", "F", "M", "Sex", "PHENO1"))
}
if (!exists("d")){
  d <- fread(PCA.eigenvec, col.names = c("FID", "IID", paste0("PC",1:PCMAX)))
  if(NORMPLOT=="y"){ 
    d[,3:(PCMAX+2)] <- lapply(d[,3:(PCMAX+2)], function(col) col / sqrt(sum(col^2)))
  }
}

text <- readLines(KEPT, n = 2)
kept <- readLines(KEPT)[-c(1, 2)]

case_ids <- cohort$IID[cohort$PHENO1==2]  # Assuming the case IDs are stored in 'IID' column of 'cases_use'
control_ids <- cohort$IID[cohort$PHENO1==1]  # Assuming the case IDs are stored in 'IID' column of 'cases_use'

cases <- d[d$IID %in% case_ids,]
controls <- d[d$IID %in% control_ids,]
cases$Group <- ifelse(cases$IID %in% kept, "Case kept", "Case removed")
controls$Group <- ifelse(controls$IID %in% kept, "Control kept", "Control removed")

#### RESULTS #####
graphstem <-sub(".*kept\\.(.*)\\.tsv$", "\\1", KEPT)
GRAPHPRE=paste0("plot.allsamples.", graphstem, ".jpg")
GRAPHSPLIT=paste0( "plot.splitsamples.",graphstem, ".jpg")
GRAPHKEPT=paste0( "plot.onlykept.", graphstem, ".jpg")


print(text)

# writing out files
outdir <- sub("(.*)/kept\\..*", "\\1", KEPT)
cat("Saving to ", outdir)

setwd(outdir)



  ## plotting
  mycols_before <- c('Controls' = "blue",
                     'Cases' = "red")
  mycolors <- c('Control kept' = "blue",
                'Control removed' ="black",
                'Case kept' = "red",
                'Case removed' = "pink")
  
  
  theme_set(theme_minimal(base_size = 12) +
              theme(
                plot.title = element_text(size = 16),
                axis.title = element_text(size = 16),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 16),
                legend.key.size = unit(0.2, "lines"),  # Adjusts the size of the legend keys
                legend.spacing.x = unit(0.2, "cm"),  # Adjusts spacing between legend entries (horizontally)
                legend.spacing.y = unit(0.5, "cm"),  # Adjusts spacing between legend entries (vertically)
                legend.box.margin = margin(1, 1, 1, 1),  # Adjusts the margin around the entire legend box
                axis.text = element_blank(), # Removes axis text
                axis.ticks = element_blank() # Removes axis ticks
              ))
  
  
  
  g1a <- ggplot() +
    geom_point(data = controls, aes(x = PC1, y = PC2, color = "Controls"), size=0.4) +
    geom_point(data = cases, aes(x = PC1, y = PC2, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC2", color = "Group") +
    coord_fixed()
  
  
  g1b <-ggplot() +
    geom_point(data=subset(controls, Group == "Control kept"), aes(x = PC1, y = PC2, color = Group), size=0.4) +
    geom_point(data=cases, aes(x = PC1, y = PC2, color = Group), size=0.4) +
    scale_color_manual(values = mycolors) +
    labs( x = "PC1", y = "PC2") +
    coord_fixed()
  
  g1c <-ggplot() +
    geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC2, color = "Controls"), size=0.4) +
    geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC2, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC2", color = "Group") +
    coord_fixed()
  
  
  g2a <- ggplot() +
    geom_point(data = controls, aes(x = PC1, y = PC3, color = "Controls"), size=0.4) +
    geom_point(data = cases, aes(x = PC1, y = PC3, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC3", color = "Group") +
    coord_fixed()
  
  g2b <- ggplot() +
    geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC3, color = Group), size=0.4) +
    geom_point(data = cases, aes(x = PC1, y = PC3, color = Group), size=0.4) +
    scale_color_manual(values = mycolors) +
    labs( x = "PC1", y = "PC3") +
    coord_fixed()
  
  g2c <-ggplot() +
    geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC3, color = "Controls"), size=0.4) +
    geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC3, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC3", color = "Group")  +
    coord_fixed()
  
  g3a <- ggplot() +
    geom_point(data = controls, aes(x = PC1, y = PC4, color = "Controls"), size=0.4) +
    geom_point(data = cases, aes(x = PC1, y = PC4, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC4", color = "Group") +
    coord_fixed()
  
  g3b <- ggplot() +
    geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC4, color = Group), size=0.4) +
    geom_point(data = cases, aes(x = PC1, y = PC4, color = Group), size=0.4) +
    scale_color_manual(values = mycolors) +
    labs( x = "PC1", y = "PC4") +
    coord_fixed()
  
  g3c <-ggplot() +
    geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC4, color = "Controls"), size=0.4) +
    geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC4, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC4", color = "Group") +
    coord_fixed()
  
  
  g4a <- ggplot() +
    geom_point(data = controls, aes(x = PC1, y = PC5, color = "Controls"), size=0.4) +
    geom_point(data = cases, aes(x = PC1, y = PC5, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC5", color = "Group") +
    coord_fixed()
  
  g4b <- ggplot() +
    geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC5, color = Group), size=0.4) +
    geom_point(data = cases, aes(x = PC1, y = PC5, color = Group), size=0.4) +
    scale_color_manual(values = mycolors) +
    labs( x = "PC1", y = "PC5") +
    coord_fixed()
  
  g4c <-ggplot() +
    geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC5, color = "Controls"), size=0.4) +
    geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC5, color = "Cases"), size=0.4) +
    scale_color_manual(values = mycols_before) +
    labs(x = "PC1", y = "PC5", color = "Group") +
    coord_fixed()

  before= arrangeGrob(g1a, g2a, g3a, g4a, nrow = 2)
  after = arrangeGrob(g1b, g2b, g3b, g4b, nrow = 2 )
  after2 = arrangeGrob(g1c, g2c, g3c, g4c, nrow = 2 )
  
  ggsave(GRAPHPRE, plot = before,device = "jpeg", path = outdir,  width = 20, height = 16, dpi = 300)
  # arrangeGrob, not grid.arrange: the latter also draws to the current device,
  # littering the output directory with Rplot00N.png
  ggsave(GRAPHSPLIT,
         plot = arrangeGrob(after, top = textGrob(paste0(text),
                                                 gp = gpar(fontsize = 20, font = 2))),
         device = "jpeg", path = outdir, width = 20, height = 16, dpi = 300)
  ggsave(GRAPHKEPT, plot = after2 ,device = "jpeg", path = outdir,  width = 20, height = 16, dpi = 300)
}


##### run_function #####
GSM_C_plotting(PCA.eigenvec = PCA.eigenvec, PSAM = PSAM, KEPT = KEPT, PCMAX = PCMAX, NORMPLOT = NORMPLOT)



