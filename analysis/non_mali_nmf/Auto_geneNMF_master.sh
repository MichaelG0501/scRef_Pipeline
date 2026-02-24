#!/bin/bash
WD=/rds/general/project/tumourheterogeneity1/ephemeral/Auto_AG
cd $WD

for ct in macrophage fibroblast endothelial nk plasma cd4 cd8; do
  qsub Auto_geneNMF_${ct}.sh
  echo "Submitted Auto_geneNMF_${ct}.sh"
done
