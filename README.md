# GPC algorithm


GPC R Codes includes R scripts to implement examples of Gibbs posterior inference on the Youden index cutoff for dementia diagnosis and diabetes diagnosis.

Additionally, the GPC_YI folder contains files needed to install an GPC_YI package with Rcpp codes for all  examples done in the above paper.

To install the GPC_YI package you will need:
  1.  R, version 3.3.4 was used to build GPC_YI.
  2.  Rtools.
  3.  The devtools package.
  
Steps:
  1.  Open Rgui.
  2.  Install the "devtools" package if you have not already using the dropdown "Packages -> Install package(s)...".  Load the devtools      package using the command library(devtools).  
  3.  Submit the command: install_github("nasyring/GPC_YI", subdir = "GPC_YI").
