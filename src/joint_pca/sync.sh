#!/bin/bash
rsync -avz -P lawless@jed.epfl.ch:/work/gr-fe/lawless/spss/exome/data/joint/pca_output/*eigen* ./data/joint/pca_output/
rsync -avz -P lawless@jed.epfl.ch:/work/gr-fe/lawless/spss/exome/data/joint/pca_output/phenotypes.csv ./data/joint/pca_output/

rsync -avz -P lawless@jed.epfl.ch:/work/gr-fe/lawless/spss/exome/src/joint/pca.R ./src/joint/
rsync -avz -P lawless@jed.epfl.ch:/work/gr-fe/lawless/spss/exome/src/joint/Rplots.pdf ./src/joint/



