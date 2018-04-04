# Rnmr1D

Rnmr1D is the main module in the NMRProcFlow web application (http://nmrprocflow.org) concerning the NMR spectra processing.

* Inside NMRProcFlow, Rnmr1D allows users to process their NMR spectra within a GUI application and thus the macro-command sequence coming from this process can be saved. 

* Outside NMRProcFlow Rnmr1D become an R package allowing users to replay  the macro-command sequence generated within NMRProcFlow. Moreover, without using NMRProcFlow, this package can also be used to replace any 'home-made script'  by a macro-command sequence.

## Installation of some dependencies

```R
source('http://bioconductor.org/biocLite.R');
biocLite('MassSpecWavelet'); biocLite('impute');
install.packages(c('doParallel', 'ptw', 'signal', 'speaq'), repos='http://cran.rstudio.com')
```

## Installation of the R package 

```R
require(devtools)
install_github("djacob65/Rnmr1D", dependencies = TRUE)
```

* Note: For Windows 7/10, it is highly recommended to proceed within R Studio so that the C++ code compilation will be correctly done.

## Quick tutorial


```R
library(doParallel)
library(Rnmr1D)

# Test with the provided example data
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "MMBBI_14P05")
CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
SAMPLEFILE <- file.path(data_dir, "Samples.txt")

# Detect the number of Cores
detectCores()

# Launch the pre-processing then the processing defined in the macro-command file
out <- Rnmr1D(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=detectCores())

# Have a look on returned data structure
ls(out)
ls(out$specMat)

# Stacked Plot with a perspective effect
plotSpecMat(out$specMat, ppm_lim=c(0.5,5))

# Overlaid Plot
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0)

# Get the data matrix 
outMat <- get_Buckets_dataset(out, norm_meth='CSN')

# Get the Signal/Noise Ratio (SNR) matrix 
outSNR <- get_SNR_dataset(out, c(10.2,10.5), ratio=TRUE)

# Get the bucket table
outBucket <- get_Buckets_table(out)

# Get the spectra data
spectra <- get_Spectra_Data(out)

```


