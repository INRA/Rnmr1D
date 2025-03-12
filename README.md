# Rnmr1D

Rnmr1D is the main module in the NMRProcFlow web application (http://nmrprocflow.org) concerning the NMR spectra processing.

* Inside NMRProcFlow, Rnmr1D allows users to process their NMR spectra within a GUI application and thus the macro-command sequence coming from this process can be saved. 

* Outside NMRProcFlow Rnmr1D become an R package allowing users to replay  the macro-command sequence generated within NMRProcFlow. Moreover, without using NMRProcFlow, this package can also be used to replace any 'home-made script'  by a macro-command sequence.

* See the [Macro-command Reference Guide](https://nmrprocflow.org/themes/pdf/Macrocommand.pdf) to have more details about macro-commands.

## Installation of some dependencies

* You may need to install a C++ compiler if not the case yet (see https://teuder.github.io/rcpp4everyone_en/020_install.html)

* Some R packages from Bioconductor:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("impute","MassSpecWavelet","pcaMethods"))
```

## Installation of the R package 

* Note for Windows 11: Before performing the installation within R GUI it may require to specify the Compiler binaries path in the PATH environment variable so that the C++ code compilation will be correctly done ( check with Sys.getenv("PATH") ). Otherwise, install [Rtools for Windows](https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html)

```R
if (!require("devtools"))
    install.packages("devtools", repos="https://cran.rstudio.com")
devtools::install_github("INRA/Rnmr1D", dependencies = TRUE)
```

## Example of use


```R
library(Rnmr1D)

# Test with the provided example data
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "CD_BBI_16P02")
CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
SAMPLEFILE <- file.path(data_dir, "Samples.txt")

# Detect the number of Cores
detectCores()

# Launch the pre-processing then the processing defined in the macro-command file
out <- Rnmr1D::doProcessing(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=detectCores())

# Have a look on returned data structure
ls(out)
ls(out$specMat)

### Stacked Plot with a perspective effect
dev.new()
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33)

### Overlaid Plot
dev.new()
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0, pY=0.1)

# Get the data matrix 
outMat <- getBucketsDataset(out, norm_meth='CSN')

# Get the Signal/Noise Ratio (SNR) matrix 
outSNR <- getSnrDataset(out, c(10.2,10.5), ratio=TRUE)

# Get the bucket table
outBucket <- getBucketsTable(out)

# Get the spectra data
spectra <- getSpectraData(out)
```

## See a more complete illustation within the vignette
```R
vignette("Rnmr1D")
```

You can access the vignette online : https://nmrprocflow.org/themes/html/Rnmr1D.html



