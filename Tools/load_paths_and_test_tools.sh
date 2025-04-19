#!/bin/bash

# Export all binary paths
export PATH=/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/Apps/hisat2-2.2.1:$PATH
export PATH=/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/Apps/sratoolkit.3.2.1-ubuntu64/bin:$PATH
export PATH=/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/Apps/subread-2.1.0-Linux-x86_64/bin:$PATH
export LD_LIBRARY_PATH=/gpfs/fs7/aafc/phenocart/PhenomicsProjects/UFPSGPSCProject/4_Assets/R/icu:$LD_LIBRARY_PATH

echo "✅ PATHs updated."

# Test the tools
echo -e "\n🔍 Testing tool installations...\n"

echo "-- hisat2 --"
hisat2 --version

echo "-- prefetch --"
prefetch --version

echo "-- fasterq-dump --"
fasterq-dump --version

echo "-- featureCounts --"
featureCounts -v

echo "-- samtools --"
samtools --version
