## R SCRIPT TO ANALYSE SCG DATASET Y. GE ET AL., 2021 ##
##            Ruben Methorst, LUMC, 2021              ##
##            Szymon M. Kielbasa, LUMC, 2021          ##

VERSION="v1.0.1"
LASTEDITDATE="2022-01-20"
SCRIPTNAME="Analysis of cervical ganglia pilot dataset LUMC"
AUTHOR="Ruben Methorst | r.methorst@lumc.nl"
THISYEAR = format(as.Date(as.POSIXlt(Sys.time())), "%Y")

# Set the working path and paths to the input files/directories
path <- "/Users/rubenmethorst/Desktop/Protocol paper/Results/"

SCG_pilot_path <- "/Users/rubenmethorst/Desktop/Protocol paper/Data/"

# Data that should be in this directory
# 1. CellRanger folder with the barcode, matrix and feature files (e.g. ../SCG_pilot/your_three_files)
# 2. An excel/csv/text file with your marker genes (Copy of gene marker list-1.xlsx)
#
# Data generated throughout this script will be saved in the same directory

## USAGE
# You can just run it in Rstudio, but it is easier to just run it in the command line
# Go to your working directory with the above specified files and run
# 'Rscript Rscript_SCG_dataset.R'

# Log file
sink(file = paste0(path, "SCG_analysis.Rlog"), split=TRUE)

# Opening statement
cat(paste0(
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
",SCRIPTNAME,"
",VERSION," - ",LASTEDITDATE,"
",AUTHOR," | 1998-",THISYEAR,".
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))

# Lets time who is faster in analyzing the data
timeStart <- Sys.time()

## Lets first install the required packages
# FUNCTION TO INSTALL PACKAGES: Courtsey of S.W. vanderlaan (swvanderlaan.github.io)
install.packages.auto <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented. 
    #update.install.packages.auto(ask = FALSE) 
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE, repos = \"https://cloud.r-project.org/\")", x)))
  }
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
    BiocManager::install() # this would entail updating installed packages, which in turned may not be warrented
    
    eval(parse(text = sprintf("BiocManager::install(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

# FUNCTION TO DO THINGS QUIETLY
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Required packages
quiet(
  suppressPackageStartupMessages(c(
    install.packages.auto("tidyverse"),
    install.packages.auto("data.table"),
    install.packages.auto("ggpubr"),
    install.packages.auto("Seurat"),
    install.packages.auto("patchwork"),
    install.packages.auto("data.table"),
    install.packages.auto("glmGamPoi"),
    install.packages.auto("crayon"),
    install.packages.auto("readxl"),
    install.packages.auto("topGO"),
    install.packages.auto("org.Mm.eg.db")
  )))

cat(bold(green("Required packages are installed and loaded!\n\nLoading data...\n")))

# Manually annotated file with marker genes
annotated_genes <- read_xlsx(paste0(SCG_pilot_path, "selected_genes.xlsx"))

## Load data and visualize QC Metrics
# Most of this analysis is according to the Seurat Guided Clustering Tutorial (https://satijalab.org/seurat/articles/SCG3k_tutorial.html#assigning-cell-type-identity-to-clusters-1)

# Read the CellRanger 10X output data (for more info: https://support.10xgenomics.com/single-cell-gene-expression/SCGftware/pipelines/latest/using/count#cr-count)
SCG.data = Read10X(data.dir = SCG_pilot_path)

# Replace the names of from sample name to logical name (Mm1-L-SCG -> HTO1f)
rownames(SCG.data$`Antibody Capture`) <- c("HTO1f", "HTO2f", "HTO3f", "HTO4f", "HTO5m", "HTO6m", "HTO7m", "HTO8m")

SCG = CreateSeuratObject(counts = SCG.data$`Gene Expression`, project = "ANS", min.cells = 3, min.features = 200)

# Calculate the amount of mitochondrial genes
SCG[["percent.mt"]] = PercentageFeatureSet(SCG, pattern = "^mt-")

cat(paste0(bold(green("\nData loaded! Now performing demultiplexing for you...\n"))))

## Demultiplexing of HTO
# Demultiplexing
SCG[["Antibody"]] <- CreateAssayObject(counts = SCG.data$`Antibody Capture`)
SCG <- NormalizeData(SCG, assay = "Antibody", normalization.method = "CLR")
SCG <- HTODemux(SCG, assay = "Antibody",verbose = FALSE)


# Clean the data (i.e. replace "Singlet" with the sample ID into a new meta.data column)
d <- bind_cols(SCG[["Antibody_classification"]], SCG[["Antibody_classification.global"]] )

d[["Antibody_classification"]] <- ifelse(
  d[["Antibody_classification.global"]] == "Doublet", "Doublet", as.character( d[["Antibody_classification"]] )
)

SCG[["Antibody_classification.cleaned"]] <- d[["Antibody_classification"]]


# Distribution of hashtag counts and heatmap of normalized counts
RidgePlot <- RidgePlot(SCG, assay = "Antibody", features = rownames(SCG[["Antibody"]]), ncol = 4)
HTO <- HTOHeatmap(SCG, assay = "Antibody",ncells = ncol(SCG[["Antibody"]])) +
  scale_fill_continuous( low = "black", high = "violet")

ggsave(HTO, file = paste0(path, "htoHeatmap.pdf"), width = 4.5*3/4*2*3, height = 8, units = "cm", dpi = 600, scale = 1 )


cat(paste0(bold(green("\nData loaded! Now performing QC for you...\n"))))

# Visualize QC metrics
fs <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
violin_plot <- VlnPlot(SCG, features = fs, ncol = length(fs), pt.size = .05, group.by = "Antibody_classification.cleaned")

ggsave(violin_plot, file = paste0(path, "QC_metrics.pdf"), width = 4.5*3/4*2*3, height = 6, units = "cm", dpi = 600, scale = 1.5 )

scatter_features <- FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
scatter_mito <-  FeatureScatter(SCG, feature1 = "nCount_RNA", feature2 = "percent.mt")


##  Perform QC on dataset
# Filter and normalize dataset
SCG <- subset(SCG, subset = nCount_RNA < 20000 & nFeature_RNA > 500 & percent.mt < 5 & Antibody_classification.global == "Singlet")
SCG <- NormalizeData(SCG)


# Counts of singlets in each sample
cat(paste0(bold(green("\nWow! We found hundreds of cells in each sample\n"))))
sort(table(SCG$Antibody_classification.cleaned),decreasing = TRUE)

cat(paste0(bold(green("\nPerforming analysis...\n"))))

# Find variable features
SCG <- FindVariableFeatures(SCG, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SCG), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SCG)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + guides(colour = F) 

ggsave(plot2, file = paste0(path, "Variabel_features.pdf"), width = 4.5*3/4*2, height = 6, units = "cm", dpi = 600, scale = 2)

 # Scale the data
all.genes <- rownames(SCG)
SCG <- ScaleData(SCG, features = all.genes)

#Linear dimension reduction
SCG <- RunPCA(SCG, features = VariableFeatures(object = SCG))

# Examine and visualize PCA results a few different ways
cat(paste0(bold(green("\nThese are the PC defining genes\n"))))
print(SCG[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(SCG, dims = 1:2, reduction = "pca")
PC_plot <- DimPlot(SCG, reduction = "pca", group.by = "Antibody_classification.cleaned")
PC_heatmap <- DimHeatmap(SCG, dims = 1:30, cells = 500, balanced = TRUE)

# # NOTE: This process can take a long time for big datasets, comment out for expediency. More
# # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# # computation time
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# 
# JackStrawPlot(pbmc, dims = 1:20)

#other method for PCA determination
Elbowplot <- ElbowPlot(SCG, ndims = 20)

ggsave(Elbowplot, file = paste0(path, "Elbowplot.pdf"), width = 4.5, height = 6, units = "cm", dpi = 600, scale = 2.5)


# We decided to continue with PCA 1:18 based on the PCA visualizations
pcas <- 1:18

#Cluster the cells
SCG <- FindNeighbors(SCG, dims = pcas)
SCG <- FindClusters(SCG, resolution = 0.4)

# If you haven't installed UMAP, you can do SCG via reticulate::py_install(packages =
# 'umap-learn')
SCG <- RunUMAP(SCG, dims = pcas)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
UMAP_plot <- DimPlot(SCG, reduction = "umap", group.by = "ident", label = TRUE) + coord_fixed() + guides(color = FALSE) + ggtitle("Cluster")
ggsave(UMAP_plot, file = paste0(path, "UMAP.pdf"), width = 4.5*5, height = 6*5, units = "cm", dpi = 600, scale = 0.75)

cat(paste0(bold(green("\nThis many clusters are found:",  length(unique(SCG@meta.data$seurat_clusters)), "\n"))))


# Make the requested violin plots
# Xist expression per sample to double-check sex
Xist_plot <- VlnPlot(SCG, features = "Xist", group.by = "Antibody_classification.cleaned",
                     ncol = 1, pt.size = 0.1) + guides(fill = F) & theme(plot.title = element_text(face = "bold.italic"))
ggsave(Xist_plot, file = paste0(path, "Xist_violin_plot.pdf"), width = 4.5, height = 6, units = "cm", dpi = 600, scale = 2.5)


# Variety of manually chosen marker genes
selGenes_clusters <- c("Th", "S100b", "Pecam1", "Acta2", "Ptprc", "Dcn")
Marker_genes_cl <- FeaturePlot(SCG, features = selGenes_clusters, ncol = 3, pt.size = 0.1) + guides(fill = F) & theme(plot.title = element_text(face = "bold.italic"))
Marker_genes_cl
ggsave(Marker_genes_cl, file = paste0(path, "Markers_UMAP_plot.pdf"), width = 4.5*4, height = 12, units = "cm", dpi = 600, scale = 1.75)


selGenes_symp <- c("Th", "Dbh", "Snap25")
Marker_genes_sy <- FeaturePlot(SCG, features = selGenes_symp, ncol = 3, pt.size = 0.1) + guides(fill = F) & theme(plot.title = element_text(face = "bold.italic"))
Marker_genes_sy
ggsave(Marker_genes_sy, file = paste0(path, "Symp_markers_UMAP_plot.pdf"), width = 4.5*4, height = 6, units = "cm", dpi = 600, scale = 1.75)

# Distribution of samples inside of the UMAP
UMAP_plot_HTO <- DimPlot(SCG, reduction = "umap", group.by = "Antibody_classification.cleaned", label = F) + coord_fixed() + ggtitle("HTO Antibodies")
ggsave(UMAP_plot_HTO, file = paste0(path, "UMAP_HTO.pdf"), width = 4.5*5, height = 6*5, units = "cm", dpi = 600, scale = 0.75)



# Plot data and save as PDF (it looks bad, but it gets the job done)
pdf(file = paste0(path, "plots_combined.pdf"), 
    paper = "a4r", 
    onefile = T, width = 80, height = 180)

plot(violin_plot)
plot(scatter_features + scatter_mito)
plot(RidgePlot)
plot(HTO)
plot(plot1 + plot2)
plot(PC_plot)
plot(PC_heatmap)
plot(Elbowplot)
plot(UMAP_plot)
plot(UMAP_plot_HTO)
plot(Xist_plot)
plot(Marker_genes_cl)
plot(Marker_genes_sy)


quiet(dev.off())

cat(paste0(bold(green("\nDone: Plots have been painted!\n"))))

# End of script
timeEnd <- Sys.time()

cat(paste0(bold(green("Finished in: ", difftime(timeEnd, timeStart, units='mins'), " Minutes.", 
                      "\nThank you for analyzing with Ruben Methorst!\n\n\n", bold(silver("\t\tByebye\n\n\n\n\n"))))))

cat(paste0(yellow(italic(
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n",
  "The MIT License (MIT)\n",
  "Permission is hereby granted, free of charge, to any perSCGn obtaining a copy of this SCGftware and asSCGciated documentation files (the \"SCGftware\"), to deal in the SCGftware without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the SCGftware, and to permit perSCGns to whom the SCGftware is furnished to do SCG, subject to the following conditions:\n",
  "    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the SCGftware.\n",
  "THE SCGFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SCGFTWARE OR THE USE OR OTHER DEALINGS IN THE SCGFTWARE.\n",
  "Reference: http://openSCGurce.org.\n",
  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"))))


# Stop the log
sink()

## Ruben Methorst, 1998-2021  ##
