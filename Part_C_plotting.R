rm(list=ls())
# INPUTS #
PCA.eigenvec="" # plink2 --PCA output, headed, with FID, IID, PCAs
PSAM <- "samplesincluded.tsv"
KEPT <- "kept.XXX"  #output of previous step (full path)

# CHOICES # 

PCMAX = 10 # number of PCs in the dataset. 
NORMPLOT <- "y" # choose y or n, to plot normalised data (which puts all PCs onto same scale)

GSM_C_plotting(PCA.eigenvec = PCA.eigenvec, PSAM = PSAM, KEPT = KEPT, PCMAX = PCMAX, NORMPLOT = NORMPLOT)



##### function #####
GSM_C_plotting <- function(PCA.eigenvec = PCA,eigenvec, PSAM = PSAM, KEPT = KEPT, PCMAX = PCMAX, NORMPLOT = NORMPLOT){
library(ggplot2)
library(grid)
library(gridExtra)


if (!exists("cohort")) {
cohort=read.table(PSAM, header=F, comment.char = "#", col.names = c("FID", "IID", "F", "M", "Sex", "PHENO1"))
}
if (!exists("d")){
  d <- read.table(PCA.eigenvec, header=F, comment.char = "#", col.names = c("FID", "IID", paste0("PC",1:PCMAX)))
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
current_date <- Sys.Date()
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
  
  
  theme_set(theme_minimal(base_size = 6) +
              theme(
                plot.title = element_text(size = 8),
                axis.title = element_text(size = 8),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 8),
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
  ggsave(GRAPHSPLIT,
         grid.arrange(after,top= textGrob(paste0(text),
                                          gp = gpar(fontsize = 10, font =2))),
         device = "jpeg", path = outdir, width = 20, height = 16, dpi = 300)
  ggsave(GRAPHKEPT, plot = after2 ,device = "jpeg", path = outdir,  width = 20, height = 16, dpi = 300)
}


##### run_function #####
GSM_C_plotting(PCA.eigenvec = PCA.eigenvec, PSAM = PSAM, KEPT = KEPT, PCMAX = PCMAX, NORMPLOT = NORMPLOT)
