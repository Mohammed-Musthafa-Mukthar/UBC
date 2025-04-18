# Set your custom library path
custom_lib <- "/gpfs/fs7/aafc/phenocart/PhenomicsProjects/UFPSGPSCProject/4_Assets/R/custom_lib"
.libPaths(custom_lib)
dir.create(custom_lib, showWarnings = FALSE, recursive = TRUE)

# Set C++14 as required by some Bioconductor packages (e.g., fgsea)
Sys.setenv("CXX14" = "g++ -std=c++14")
Sys.setenv("CXX14STD" = "-std=c++14")
Sys.setenv("CXXFLAGS" = "-std=c++14")

# Install missing CRAN dependencies first
cran_pkgs <- c("ggplot2", "ggforce", "ggrepel", "viridis", "ggfun", "aplot",
               "ggnewscale", "shadowtext", "scatterpie", "qvalue")

install.packages(cran_pkgs,
                 repos = "https://cloud.r-project.org",
                 lib = custom_lib,
                 dependencies = TRUE)

# Ensure BiocManager is present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = custom_lib)
}

# Use BiocManager to install Bioconductor packages
bioc_pkgs <- c("fgsea", "DOSE", "ggraph", "ggtree", "enrichplot", "clusterProfiler", "org.Mm.eg.db")

BiocManager::install(bioc_pkgs,
                     lib = custom_lib,
                     ask = FALSE,
                     update = TRUE,
                     dependencies = TRUE)

# Final confirmation
cat("✅ Final installed packages:\n")
print(installed.packages(lib.loc = custom_lib)[, "Package"])
