# CRCmicrobiome
Gut microbiome-based machine learning model for colorectal cancer (CRC), including a large number of young/early onset CRC.
The related results were published in _Nature Communications_ https://doi.org/10.1038/s41467-024-47523-x. 

**Sample information**
- Public data: 1262 samples from 10 studies, which includes 660 CRC and 662 controls. The integrated data file was obtained from Beghini _et al_.
- Fudan study: 200 samples from a Fudan study Yang _et al_, which incudes 100 CRC (**50 samples from patients under age 50**) and 100 age-matched controls.
- Unpublished study: 460 newly sequencing CRC samples from patients recruited in a Guangzhou hospital (**167 samples from patients under age 50**).

**Data type**
- All samples were shotgun metagenomics.
- Speceis profiles were generated by MetaPhlAn3.
- Functional profiles (pathway & level4 ec) were generated by HUMAnN3.

**Machine learning model**
- random forest
- lasso logistic regression
  
**Essential R packages**
- mlr3: https://mlr3book.mlr-org.com/
- SIAMCAT: https://bioconductor.org/packages/release/bioc/html/SIAMCAT.html
- tested on R version 4.3.1 (2023-06-16) -- "Beagle Scouts"_
