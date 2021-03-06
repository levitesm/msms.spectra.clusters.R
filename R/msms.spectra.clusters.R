#' @details This package allows you to create spectral libraries from the .mgf files by clustering spectra according to their similarity. You need to prepare the MGF object form the .mgf file using the PrepareMgf function. Then you Create a library from mgf. After that you can add identifications for you library, and new MGFs to it, get the FDR table for the clusters, and check the similarity between different clusters (getClustersSimilarity) to look for PTMs(getModif).
#' @usage Say, you have a large .mgf file that you want to use for creating a spectral library for you organism called LibMgf.mgf. You have th id_dt data.table which contains the known identifications for spectra in LibMgf.mgf (table should contain Title, Sequence and Conf. The title of the spectrum in the .mgf file, it's sequence and a confidence of this identification. In case of conflict the sequence with a higher confidence will be used). You also have a file Mgf2.mgf which you would like to identify using your library.
#' First you prepare you LibMgf.mgf for clustering
#' libMgf<-Prepare_MGF('../LibMgf.mgf')
#' Then you create a Library
#' lib<-CreateLibFromMGF(libMgf)
#' You'll get a report for the clustering process after this step.
#' 
#' To add Identifications to your library you use:
#' lib<-addIDsToLibrary(lib,id_dt)
#' The IDs will be automatically extended to the members of the same cluster.
#' TO get the clustered MGF with its cluster members and identifications use:
#' mgf<-getFullMgf(lib)
#' 
#' To identify a new mgf with the help of your library you'll need to add its spectra to the clusters:
#' 
#' mgf2<-Prepare_MGF('../Mgf2.mgf')
#' lib<-add_MGF_to_LIBRARY(mgf2, lib)
#' After this step you'll get a report of the similarity of the new mgf to the library.
#' 
#' #' To get the added MGF and its identifications use:
#' dt<-getLastMgfFromLib(lib)
#' 
#' 
#' #' To see the table of all the spectra in the library with their cluster members, identifications and extended identifications and FDR for each cluster, use:
#' dt<-getFullMgf(lib)
#' 
#' In order to look for potential PTMs in your dataset, you need to look for pairs of similar clusters. 
#' simDt<-getClustersimilarity(lib)
#' 
#' Use getModif(simDt) to get the actual PTM masses.
"_PACKAGE"
#> [1] "_PACKAGE"