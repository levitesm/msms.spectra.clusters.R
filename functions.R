
#################################################
############# FUNCTIONS #########################
#################################################


################
#### Parsing MGF
#' A Parse MGF function
#'
#' This function allows you to read an .mgf file and convert it into a data.table
#' @param path A path to the .mgf file
#' @param mode The mode of MGF reading. Can be "PP" for a ProteinPilot MGF, or "centroid" for a CEnrtoid MGF. Defaults to centroid.
#' @keywords mgf 
#' @examples mgf_data_table <- Parse_MGF('C:/docs/mgf1.mgf')
Parse_MGF <- function(path, mode = 'centroid'){
  
  stopifnot(require(data.table))
  
  start_time <- Sys.time()
  require(data.table)
  require(stringi)
  
  tryCatch(path <- normalizePath(path), warning = function(w)warning(w), error = function(e)warning(e))
  
  print(paste0('Start reading at ', round(as.numeric(Sys.time() - start_time), 1)))
  
  mgf_dt <- fread(path, stringsAsFactors = F, header = F, colClasses = 'c', sep = '\n', blank.lines.skip = T,
                  verbose = T, showProgress = T, strip.white = T) # slow; don't know a faster way to read
  
  print(paste0('Finished reading at ', round(as.numeric(Sys.time() - start_time), 1)))
  
  
  mgf_dt<-mgf_dt[!grepl('###',V1)]
  
  line_ct <- mgf_dt[, .N]
  mgf_dt[, line_num := 1:line_ct]
  
  # -------- substr is faster than both str_sub and stri_sub
  mgf_dt[, is_ion := substr(V1, 1, 1) %in% as.character(0:9)] # According to Matrix Science website,
  # any line starting with numbers describes an ion,
  # any line starting with words describes spectrum or dataset
  
  print(paste0('Ions detected at ', round(as.numeric(Sys.time() - start_time), 1)))
  
  setkey(mgf_dt, is_ion)
  mgf_tags_dt <- setkey(mgf_dt[J(F), .(line_num, V1)], V1)
  mgf_dt <- setkey(mgf_dt[J(T), list(line_num, V1)], line_num)
  
  begins <- sort(mgf_tags_dt['BEGIN IONS', line_num])
  # ------ this is >5 times faster than begins <- sort(mgf_dt[V1 == 'BEGIN IONS', line_num]) because > 95% of mgf is ions, not tags
  
  ends <- sort(mgf_tags_dt['END IONS', line_num])
  
  print(paste0('Begins ends detected at ', round(as.numeric(Sys.time() - start_time), 1)))
  
  if (length(begins) == 0) {
    warning('No spectra found in your file (based on BEGIN IONS TAG).')
    return(NULL)
  }
  
  if (length(begins) != length(ends)) {
    warning('Problem with tags: number of BEGIN IONS tags is not the same as END IONS tags.')
    return(NULL)
  }
  
  lengths <- ends - begins + 1
  
  if (any(lengths <= 1)) {
    warning('Problem with tags: some BEGIN IONS and END IONS tags have wrong order.')
    return(NULL)
  }
  
  
  spec_ct <- length(begins)
  
  spec_index <- setkey(data.table(spec_num = rep.int(1L:spec_ct,
                                                     times = lengths),
                                  line_num = 1:line_ct), line_num)
  
  
  setkey(mgf_tags_dt, line_num)
  
  mgf_dt <- setkey(mgf_dt[spec_index, .(spec_num, V1), nomatch = 0], spec_num)
  mgf_tags_dt <- setkey(merge(mgf_tags_dt, spec_index, all = F)[, tag := sapply(stri_split(V1, fixed = '=',
                                                                                           n = 2), '[', 1)], tag)
  
  rm(spec_index)
  
  print(paste0('Tags detected at ', round(as.numeric(Sys.time() - start_time), 1)))
  
  if (mode != 'centroid'){
    mgf_tags_dt <-
      dcast.data.table(mgf_tags_dt[!J(c('BEGIN IONS', 'END IONS'))][spec_num != 0 & spec_num != spec_ct + 1], spec_num ~ tag,
                       value.var = 'V1')[, ':='(charge = as.integer(apply(stri_match_first_regex(CHARGE, '^CHARGE=([[:digit:]]+)([+-]?)$')[, 3:2],
                                                                          1, stri_join, collapse = '')),
                                                rt = as.integer(stri_match_first_regex(RTINSECONDS, '^RTINSECONDS=([[:digit:]]+)$')[, 2]),
                                                prec_mz = as.numeric(stri_match_first_regex(PEPMASS, '^PEPMASS=([[:digit:]]+\\.?[[:digit:]]*)$')[,2]),
                                                title = stri_match_first_regex(TITLE, '^TITLE=(.*)$')[, 2],
                                                CHARGE = NULL, PEPMASS = NULL, RTINSECONDS = NULL, TITLE = NULL
                       )]
    
    # Charge has sign, can be faster if done by CHARGE, but already very fast here
    # Time expected to be integer
    # Only expects mass to be like 100.111
    
  } else {
    mgf_tags_dt <-
      dcast.data.table(mgf_tags_dt[!J(c('BEGIN IONS', 'END IONS'))][spec_num != 0 & spec_num != spec_ct + 1], spec_num ~ tag,
                       value.var = 'V1')[, ':='(prec_mz = as.numeric(stri_match_first_regex(PEPMASS, '^PEPMASS=([[:digit:]]+\\.?[[:digit:]]*)$')[,2]),
                                                title = stri_match_first_regex(TITLE, '^TITLE=(.*)$')[, 2],
                                                PEPMASS = NULL, TITLE = NULL
                       )]
  }
  print(paste0('Tags parsed at ', round(as.numeric(Sys.time() - start_time), 1)))
  
  if (mode == 'PP' | mode == 'centroid') {
    mgf_tags_dt[, c('locus', 'file') := {
      tmp <- stri_match_first_regex(title, pattern = '^Locus:([[:digit:].]*) File:\"(.*)\"$')
      .(tmp[, 2], as.factor(tmp[, 3]))
    }]
  }
  
  print(paste0('Locuses parsed at ', round(as.numeric(Sys.time() - start_time), 1)))
  # ------
  if (mode != 'centroid') {
    res <- mgf_dt[, c('frg_mz', 'frg_int', 'frg_z', 'V1') :=
    {
      split <- stri_split_fixed(V1, pattern = ' ', n = 3, simplify = T);
      .(as.numeric(split[, 1]), as.numeric(split[, 2]), as.integer(split[, 3]), NULL)
    }][, {order <- order(frg_mz); .(frg_mz = list(frg_mz[order]), frg_int = list(frg_int[order]), frg_z = list(frg_z[order]))}, spec_num][mgf_tags_dt]
  } else {
    res <- mgf_dt[, c('frg_mz', 'frg_int', 'V1') :=
    {
      split <- stri_split_fixed(V1, pattern = ' ', n = 2, simplify = T);
      .(as.numeric(split[, 1]), as.numeric(split[, 2]), NULL)
    }][, {order <- order(frg_mz); .(frg_mz = list(frg_mz[order]), frg_int = list(frg_int[order]))}, spec_num][mgf_tags_dt]
  }
  
  # slow, but stri_split_fixed and n = 3 is 5 time faster then base strsplit; this is the slowest part (45% of time)
  
  print(paste0('Ions detected at ', round(as.numeric(Sys.time() - start_time), 1)))
  
  setkey(res, spec_num)
}##### END Parsing MGF



####################
### Noise estimation
#' A noise estimation function
#'
#' This function estimates noise level for a list of peak intensities. Supports Centroin and CWT peak pickings.
#' @param int A list of peak intensities form MSMS spectrum.
#' @keywords noise 
#' @examples noise <- noise_estimation(c(1,2,3,4,5))
noise_estimation <- function(int){
  
  l <- length(int)
  if (l < 20) {
    return(list(bw = NA_real_, noise = NA_real_))
  }
  
  si<-sort(int)[1:10]
  if(min(si)==max(si))
  {
    return(list(bw = NA_real_, noise = si[1]*4))### NOISE = Min_Noise * 4
  }
  
  diff_dens <- density(as.matrix(dist(int))[lower.tri(matrix(nrow = l, ncol = l))], bw = 0.01)
  bw <- diff_dens$x[which.max(diff_dens$y)]
  
  if (bw <= 0) {
    bw<-0.01
    
  }
  
  d <- density(int, bw = bw)
  der <- diff(d$y)/diff(d$x)
  start_search <- d$x[which.min(der)]
  noise <- d$x[2:length(d$x)][which(der > min(der)/10 & d$x[2:length(d$x)] > start_search)[1]]
  return(list(bw = bw, noise = noise))
}### END Noise estimation


###############################
### Correlation between 2 lists
#' A Parse MGF function
#'
#' This function calculates the correlation between two lists of peak masses. 
#' @param mass_list1 The first list of masses.
#' @param mass_list2 The secound list of masses.
#' @param D A maximal difference between two masses for them to be considered the same. Delaults to 0.04
#' @param mode If mode=1 the correlation score will be returned, else the number of correlationg peaks is returned. The correlation score is the percent of correlating peaks out of the mean length of two lists. Defaluts to 1.
#' @keywords correlation 
#' @examples corr <- Spec_list_corr(list1, list2)
Spec_list_corr<-function(mass_list1, mass_list2, D = 0.04, mode = 1)
{
  
  if(length(mass_list1)>length(mass_list2))
  {
    small_list<-mass_list2
    big_list<-mass_list1
  }
  else
  {
    small_list<-mass_list1
    big_list<-mass_list2
  }
  
  the_count<-0
  xj<-1
  for(i in 1: length(small_list))
  {
    mass1<-small_list[i]
    
    if(xj>length(big_list)){break}
    
    for(j in xj:length(big_list))
    {
      if(big_list[j]>=(mass1-D))
      {
        if(big_list[j]<=(mass1+D)){the_count<-the_count+1;xj<-j+1}
        break
      }
      else{xj<-j+1}
    }
  }
  
  
  score<-the_count/(length(small_list)+length(big_list))*200
  
  if (mode==1) return(score) else return(the_count)
  
}### END Correlation



#################
### Write MGF Centroid
#' A Write MGF Centroid function
#'
#' This function allows you to write an .mgf data.table into an .mgf file. MGF should be of a Centroid type.
#' @param mgf The mgf data.table object
#' @param path A path to the .mgf file
#' @keywords write 
#' @export
#' @examples write_mgf_centroid(mgf1,'C:/docs/mgf1.mgf')
write_mgf_centroid <- function(mgf, path) {
  
  stopifnot(require(data.table))
  
  for(n in 1:nrow(mgf)){
    write(x = 'BEGIN IONS', file = path, ncolumns = 1,append = T, sep = '\n');
    write(x = paste0('TITLE=', mgf[n, title]), file = path, ncolumns = 1,append = T, sep = '\n');
    
    if (is.finite(mgf[n, prec_mz])) {
      write(x = paste0('PEPMASS=', round(mgf[n, prec_mz], 5)), file = path, ncolumns = 1,append = T, sep = '\n')
    }
    
    if (mgf[n, frg_count] > 0) {
      
      
      fgmz<-mgf[n, frg_mz][[1]]
      fgint<-mgf[n, frg_int][[1]]
      
      vec<-c(1:length(fgmz))
      
      txt<-sapply(vec, function(i){
        return(paste0(round(fgmz[i], 4), ' ', round(fgint[i], 4),''))
      })
      
      write(x = txt, file = path, ncolumns = 1,append = T, sep = '')
    }
    write(x = 'END IONS', file = path, ncolumns = 1,append = T, sep = '\n');
  }
  
}### END Write MGF Centroid





#################
### Write MGF PP
#' A Write MGF PP function
#'
#' This function allows you to write an .mgf data.table into an .mgf file. MGF should be of a ProteinPilot type.
#' @param mgf The mgf data.table object
#' @param path A path to the .mgf file
#' @keywords write 
#' @export
#' @examples write_mgf_(mgf1,'C:/docs/mgf1.mgf')
write_mgf_PP <- function(mgf, path) {
  
  stopifnot(require(data.table))
  
  for(n in 1:nrow(mgf)){
    write(x = 'BEGIN IONS', file = path, ncolumns = 1,append = T, sep = '\n');
    write(x = paste0('TITLE=', mgf[n, title]), file = path, ncolumns = 1,append = T, sep = '\n');
    write(x = paste0('CHARGE=', mgf[n, charge],'+'), file = path, ncolumns = 1,append = T, sep = '\n');
    
    
    
    if (is.finite(mgf[n, prec_mz])) {
      write(x = paste0('PEPMASS=', round(mgf[n, prec_mz], 5)), file = path, ncolumns = 1,append = T, sep = '\n')
    }
    
    write(x = paste0('RTINSECONDS=', mgf[n, rt]), file = path, ncolumns = 1,append = T, sep = '\n');
    
    
    if (mgf[n, frg_count] > 0) {
      
      fgmz<-mgf[n, frg_mz][[1]]
      fgint<-mgf[n, frg_int][[1]]
      
      vec<-c(1:length(fgmz))
      
      txt<-sapply(vec, function(i){
        return(paste0(round(fgmz[i], 4), ' ', round(fgint[i], 4),''))
      })
      
      write(x = txt, file = path, ncolumns = 1,append = T, sep = '')
    }
    write(x = 'END IONS', file = path, ncolumns = 1,append = T, sep = '\n');
  }
  
}### END Write MGF PP



#################
###Y List from MGF
#' An Y list form MGF function
#'
#' This function recieves an mgf data.table object and extracts the mass lists from each spectrum accourding to the number of Y ions defined inside the MGF.
#' @param mgf The mgf data.table object
#' @keywords Y_list 
#' @examples mass_list<-Y_list_from_mgf(mgf1)
Y_list_from_mgf<-function(mgf)
{
  
  stopifnot(require(data.table))
  
  spectra <-
    mgf[, 
        .(frg_mz = {
          N <- min(number,frg_above_noise_count)
          if (is.na(N)) {
            N <- 20
          }
          l<-sort(head(
            frg_mz[[1]][order(frg_int[[1]], decreasing = T)],
            N))
          list(l)
        }
        
        ),
        title]
  
  
  mass_list <- spectra[, frg_mz]
  
  return(mass_list)
}### END Y List from MGF


#################
### Anti MGF
#' An Anti MGF function
#'
#' This function recieves an mgf data.table object and adds the decoy specta to it. Each spectrum gets an ANTI spectum which is used as a decoy for clustering algorithm.
#' @param mgf The mgf data.table object
#' @keywords anti 
#' @examples mgf_with_decoys<-addANTI(mgf1)
addANTI<-function(mgf){
  
  stopifnot(require(data.table))
  
  mycoMGF[,anti:=F]
  anti_mgf<-mycoMGF
  anti_mgf[,anti:=T]
  
  
  #NEW
  K_COUNT<-0
  R_COUNT<-0
  BAD1<-0
  BAD2<-0
  
  
  nnn<-c()
  nni<-c()
  for(i in 1: nrow(anti_mgf))
  {
    ls<-unlist(anti_mgf[i]$frg_mz)
    li<-unlist(anti_mgf[i]$frg_int)
    
    N <- min(anti_mgf[i]$number,anti_mgf[i]$frg_above_noise_count)
    if (is.na(N)) {
      
    }
   
    max_ls<-sort(head(
      ls[order(li, decreasing = T)],
      N))
    
    which_list<-which(ls %in% max_ls)
    
    max_li<-li[which_list]
    
    rest_ls<-ls[-which_list]
    rest_li<-li[-which_list]
    
    
    
    fan<-anti_mgf[i]$frg_mz[[1]]
    MASS<-anti_mgf[i]$prec_mz*anti_mgf[i]$charge-anti_mgf[i]$charge
    
    
    if(tail(fan[fan<(MASS-173)],n = 1)>(MASS-173.2) | (fan[fan>175][1]<175.2  & !(tail(fan[fan<(MASS-145)],n = 1)>(MASS-145.2)))) 
    {
      fix_max_ls<-abs(MASS-max_ls+2+18+156.10111)# Y and R
      fix_rest_ls<-abs(MASS-rest_ls+2-18-156.10111)# B and R
      R_COUNT<-R_COUNT+1
    }else
    {
      fix_max_ls<-abs(MASS-max_ls+2+18+128.09496)# Y and K
      fix_rest_ls<-abs(MASS-rest_ls+2-18-128.09496)# B and K
      K_COUNT<-K_COUNT+1
    }
   
    
    
   
    
    
    nls<-list(c(fix_max_ls,rest_ls))
    nli<-list(c(max_li,rest_li))
    
    dt<-data.table(x=unlist(nls),y=unlist(nli))
    dt<-dt[order(x)]
    
    nnn<-c(nnn,list(dt$x))
    nni<-c(nni,list(dt$y))
  }
  ####
  
  
  if(length(nnn)==1)## NEW FIX
  {
    anti_mgf$frg_mz<-list(nnn)
    anti_mgf$frg_int<-list(nni)     
  }else
  {
    anti_mgf$frg_mz<-nnn
    anti_mgf$frg_int<-nni
  }
  
  anti_mgf[,title:=paste0('ANTI__',title)]
  
  anti_mgf[, frg_above_noise := .(list(frg_mz[[1]][frg_int[[1]] > noise])), spec_num]
  
  
  mycoMGF[,anti:=F]
  
  fmg<-rbind(mycoMGF,anti_mgf)
  
  return(fmg)
}### END Anti




#################
### Add Spec to LIBRARY
#' An add Spec to library function
#' This function recieves a library object created by the CreateLibFromMGF function and a spectrum which is to be added to the library. It returns library with the new spectrum added.
#' @param library A spectral library for Spectrum addition
#' @param Spec A spectrum to be added to library. A single entity for an MGF object.
#' @param D_MS a maximal difference between precoursor masses of two spectra for them to be considered the same. Defaults to 0.04.
#' @param D_MSMS a maximal difference between fragmentation spectra peaks' masses from them to be considered the same. Delaults to 0.04.
#' @param MIN_PEAK_COUNT a minimal number of correlating peaks between two spactra for them to be consideredd the same. Defaults to 5.
#' @keywords addSpectrum 
#' @examples library1<-addSpecToLibrary(library1,mgf1[100])
addSpecToLibrary<-function(library,Spec, D_MS = 0.04, D_MSMS = 0.04, MIN_PEAK_COUNT = 5)
{
  
  stopifnot(require(data.table))
  
  clusters<-library[[1]]
  fullMgf<-library[[2]]
  CLUSTER_CUT<-library[[3]]
  mgf<-Spec
  
  mgf[,cluster_num:=0]
  mgf[,anti_num:=0]
  
  ### START Adding ANTI
  mgf<-addANTI(mgf)
  
  mgf_list<-Y_list_from_mgf(mgf)
  
  mgf_prec<-mgf$prec_mz
  mgf_charge<-mgf$charge
 
  
  for(i in 1:2){
    
  list1<-unlist(mgf_list[i])
  
  if(length(clusters)>0){
    
    #Regular

    scores<-rep(0,length(clusters))
    for(j in 1: length(clusters))
    {
      if(mgf_charge[i]==clusters[j][[1]]$cons_charge & abs(mgf_prec[i]-clusters[j][[1]]$cons_prec)<=D_MS ) # << PREC CUT >> MAYBE NUMBER!!!
      {
        list2<-unlist(clusters[j][[1]]$cons_spec)

        scores[j]<-Spec_list_corr(list1,list2,D = D_MSMS)

      }
      else next

    }
    
   
    
  }else{scores=0}#CLUSTERS EMPTY
  
  
  PC<-0
  if(max(scores)>CLUSTER_CUT)
  {
  X<-which.max(scores)
  PC<-Spec_list_corr(unlist(clusters[X][[1]]$cons_spec), list1, mode = 2)
  }
  
  if(max(scores)>CLUSTER_CUT & PC>MIN_PEAK_COUNT) # << CLUSTER CUT >>
  { #Add Spec to CLUSTER
    
    
    ### ADDITIONAL CHECK
    in_clust_score<-c(0)
    in_clust_count<-c(0)
    if(length(clusters[X][[1]]$SPECS)>0) ### Added minimal similarity
    {
      iii<-c(1:length(clusters[X][[1]]$SPECS))
      in_clust_score<-lapply(iii, function(ii){
        spc<-clusters[X][[1]]$SPECS[ii]
        return(Spec_list_corr(list1,unlist(spc),D = D_MSMS))
      })
      in_clust_score<-unlist(in_clust_score)
    
    
    if (min(in_clust_score)>CLUSTER_CUT) #Secound Part
    {
      iii<-c(1:length(clusters[X][[1]]$SPECS))
      in_clust_count<-lapply(iii, function(ii){
        spc<-clusters[X][[1]]$SPECS[ii]
        return(Spec_list_corr(list1,unlist(spc),D = D_MSMS,mode = 2))#>> PEAK COUNT
      })
      in_clust_count<-unlist(in_clust_count)
    }}
    
    if (min(in_clust_count)>MIN_PEAK_COUNT) #ADD to cluster //SIMILAR TO ALL
    {
      
      clusters[X][[1]]$NUM<-clusters[X][[1]]$NUM+1
      if(mgf[i]$anti)
      {
        clusters[X][[1]]$ANTI_COUNT<-clusters[X][[1]]$ANTI_COUNT+1
        mgf[1,anti_num:=X]
      }
      else
      {
        mgf[1,cluster_num:=X]
      }
      clusters[X][[1]]$TITLES<-c(clusters[X][[1]]$TITLES, mgf[i]$title)
      clusters[X][[1]]$SPECS<-append(clusters[X][[1]]$SPECS,list(list1))
      clusters[X][[1]]$MINS<-c(clusters[X][[1]]$MINS,min(in_clust_score))
      clusters[X][[1]]$MEANS<-c(clusters[X][[1]]$MEANS,mean(in_clust_score))  
      
      if(mgf[i]$sorter>clusters[X][[1]]$cons_sorter)#Changing the Cons Spec
      {
        clusters[X][[1]]$cons_spec<-list1
        clusters[X][[1]]$cons_title<-mgf[i]$title
        clusters[X][[1]]$cons_sorter<-mgf[i]$sorter
        clusters[X][[1]]$cons_prec<-mgf_prec[i]
      }
    }
    else
    { #Create new CLUSTER
      if(!mgf[i]$anti){
        mgf[1,cluster_num:=(length(clusters)+1)]
        new_clust<-list(cons_spec=list1, cons_prec=mgf_prec[i], cons_charge=mgf_charge[i], cons_title=mgf[i]$title, NUM=1, ANTI_COUNT=0, TITLES=list(mgf[i]$title) , SPECS=list(list1), MINS=c(100), MEANS=c(100), cons_sorter=mgf[i]$sorter )
        
        if(length(clusters)==0)
        {
          clusters<-list(new_clust)
        }
        else
        {
          clusters[length(clusters)+1]<-list(new_clust)
        }
        
        
      }
      
    }
  }
  else
  {#NEW CLUSTER
    if(!mgf[i]$anti){
      mgf[1,cluster_num:=(length(clusters)+1)]
      new_clust<-list(cons_spec=list1, cons_prec=mgf_prec[i], cons_charge=mgf_charge[i], cons_title=mgf[i]$title, NUM=1, ANTI_COUNT=0, TITLES=list(mgf[i]$title) , SPECS=list(list1), MINS=c(100), MEANS=c(100), cons_sorter=mgf[i]$sorter )
      
      if(length(clusters)==0)
      {
        clusters<-list(new_clust)
      }
      else
      {
        clusters[length(clusters)+1]<-list(new_clust)
      }
      
      
    }
    
  } 
  }
  
  mgf[,anti:=NULL]
  if(nrow(fullMgf)>0){fullMgf<-rbind(fullMgf,mgf[1], fill = T)}else{fullMgf<-mgf[1]}
  library<-list(clusters,fullMgf,CLUSTER_CUT)
  return(library)
}### END Add Spec to Library




#################
### Add MGF to Library
#' An add MGF to library function
#' This function recieves a library object created by the CreateLibFromMGF function and an MGF object created by PrepareMGF function which is to be added to the library. It returns library with the new mgf added.
#' @param mgf An mgf to be added to library.
#' @param library A spectral library for MGF addition
#' @param sort_type A way of sorting the MGF before addition. Can be: none - no sorting, quality - sorting be spectra quality form high to low, random - randomal sorting. Defaults to quality.
#' @keywords addMgf
#' @export
#' @examples library1<-add_MGF_to_LIBRARY(library1,mgf1)
add_MGF_to_LIBRARY<-function(mgf, library, sort_type='quality')
{
  
  stopifnot(require(data.table))
  
  print(paste0('>>> ',Sys.time(), '   >>> START CLUSTERING'))
  
  if(sort_type=='quality')
  {
    mgf<-mgf[order(sorter,decreasing = T)]
  }
  else if(sort_type=='none')
  {
    # Do nothing
  }
  else if(sort_type=='random')
  {
    mgf<-mgf[sample(1:nrow(mgf)),]
  }
  else
  {
    print('Sort_type not supported. Performing NO sort.')
  }
 
  print(paste0('>>> ',Sys.time(), '   >>> Adding ', nrow(mgf) , " Specs to library"))
  
  N1<-length(library[[1]])
  
  print(paste0('>>> ',Sys.time(), '   >>> Clusters in library - ', N1))
  
  for(i in 1:nrow(mgf))
  {
     
    library<-addSpecToLibrary(library,mgf[i])#,cl)
    if(i%%1000==0)print(i)
    
  }
  
  N2<-length(library[[1]])
  
  print(paste0('>>> ',Sys.time(), '   >>> CLusters after adding --- ', N2))
  
  print(paste0('>>> ',Sys.time(), '   >>> NEW Clusters --- ', (N2-N1)))
  
  clusters<-library[[1]]
  
  sum<-0
  
  for(i in (N1+1): length(clusters))
  {
    sum<-sum+clusters[i][[1]]$NUM
    
  }
  
  lib_num<-nrow(mgf)-sum
  
  print(paste0('>>> ',Sys.time(), '   >>> ', lib_num, ' spectra fell into the initial library. -- ' ,(lib_num/nrow(mgf)) ,'%'))
  
  
  print(paste0('>>> ',Sys.time(), '   >>> FINISH CLUSTERING'))
  
  library<-addFDRtoLibrary(library)
  
  library<-expendLibIDs(library)
  
  library[[4]]<-mgf$title
  
  return(library)
  
}### END Clusters




#################
### CreateLibFromMGF
#' A create librari form mgf function
#' This function recieves an MGF object created by PrepareMGF function and creates a spectral library out of it.
#' @param mgf An mgf to be added to library.
#' @param CLUSTER_CUT A minimal level of correlation between the two spactra from them to be considered the same. Defaults to 20.
#' @keywords createLibrary
#' @export
#' @examples library1<-CreateLibFromMGF(mgf1)
CreateLibFromMGF<-function(mgf, CLUSTER_CUT=20)
{
  
  stopifnot(require(data.table))
  
  print(paste0('>>> ',Sys.time(), '   >>> START CLUSTERING'))
  
  ### HAVE TO SORT by Quality
  mgf<-mgf[order(sorter,decreasing = T)]
  
  
  clusters<-list()
  m<-mgf[0:0]
 
  library<-list(clusters,m,CLUSTER_CUT)
  
  for(i in 1:nrow(mgf))
  {
    library<-addSpecToLibrary(library,mgf[i])#,cl)
    if(i%%1000==0)print(i)
  }
  
  print(paste0('>>> ',Sys.time(), '   >>> FINISH CLUSTERING'))
  
  library<-addFDRtoLibrary(library)
  library[[4]]<-mgf$title
  
  return(library)
  
}### END MGF to Clusters



#################
### Prepare MGF
#' A Prepare mgf function
#' This function recieves a path to the .mgf file and returnes an MGF data.table object. Mgf is parsed using Parse_mgf and nise_estimation functions, and prepared to clustering.
#' @param path A path to the .mgf file
#' @param mode A mode of mgf preparation. Can be PP for Protein PIlot type of mgf, and centroid for the centroid mgf type. Defaults to centroid.
#' @keywords prepareMgf
#' @export
#' @examples mgf1<-Prepare_MGF('c:/docs/mgf1.mgf')
Prepare_MGF<-function(path,mode='centroid')
{
  
  stopifnot(require(data.table))
  library(snowfall)
  
  print(paste0('>>> ',Sys.time(), '   >>> START Parsing MGF'))
  
  mgf <- Parse_MGF(path, mode = mode)
  
  print(paste0('>>> ',Sys.time(), '   >>> MGF Size = ', nrow(mgf)))
  
  ints <- mgf$frg_int
  
  print(paste0('>>> ',Sys.time(), '   >>> START Noise Estimation'))
  
  sfInit(parallel = T, cpus = 12)
  noise_bw <- sfLapply(ints, noise_estimation)
  
  sfStop()
  gc()
  
  noise_bw_corr <- lapply(noise_bw, function(x){if (length(x) == 1) {return(list(noise = NA_real_, bw = NA_real_))} else {return(x)}})
  
  noises <- sapply(noise_bw_corr, '[[', 'noise')
  bws <- sapply(noise_bw_corr, '[[', 'bw')
  
  mgf[, bw := bws]
  mgf[, noise := noises]
  
  rm(list = c('bws', 'noises', 'noise_bw'))
  
  mgf <- setkey(mgf, locus)
  
  ############################
  mgf<-mgf[!is.na(noise)]
  
  mgf[, spec_num:=c(1:nrow(mgf))]
  mgf[, frg_count := length(frg_mz[[1]]), spec_num]
  mgf[, frg_above_noise := .(list(frg_mz[[1]][frg_int[[1]] > noise])), spec_num]
  mgf[, frg_above_noise_count := length(frg_above_noise[[1]]), spec_num]
  mgf[, ratio := (mean(frg_int[[1]][frg_int[[1]]>noise])/mean(frg_int[[1]][frg_int[[1]]<=noise])), spec_num]
  
  mgf<-mgf[!is.na(ratio)]
  mgf[,sorter:=frg_above_noise_count*ratio, spec_num]
  
  if(mode=='PP')
  {
    mgf[is.na(charge),charge:=max(1,min(5,ceiling(max(unlist(frg_above_noise))/prec_mz))),spec_num]
  }
  if(mode=='centroid')
  {
    mgf[,charge:=max(1,min(5,ceiling(max(unlist(frg_above_noise))/prec_mz))),spec_num]
  }
  
  mgf[,number:=round(prec_mz*charge/111*2),spec_num] # Averagine mass = 111.117
  
  mgf<-mgf[number>10 & frg_above_noise_count>10] # CHOOSE GOOD
 
  print(paste0('>>> ',Sys.time(), '   >>> MGF Size GOOD SPECS = ', nrow(mgf)))
  
  print(paste0('>>> ',Sys.time(), '   >>> FINISH Prepare MGF'))
  
  return(mgf)
  
  
} ### END PRepare MGF


#################
### Get FDR Table
#' A get FDR table function
#' This function recieves a spectral library created be CreateLIBfromMGF function and returns a table of FDR values by cluster. 
#' @param library A spectral library
#' @keywords FDR
#' @examples FDRtable<-getFDRtable(lib1)
getFDRtable<-function(library)
{
  
  stopifnot(require(data.table))

  clusters<-library[[1]]
  
  nums<-c(1:length(clusters))
  dt<-data.table(cluster_num=nums, FullCount=0,AntiCount=0)
  for(i in 1: length(clusters))
  {
    dt[i]$FullCount=clusters[i][[1]]$NUM
    dt[i]$AntiCount=clusters[i][[1]]$ANTI_COUNT
  }
  
  dt[,FDR:=((AntiCount/FullCount*(FullCount-AntiCount)+AntiCount)/FullCount)] 
  return(dt)
}### END Get FDR



#################
### Add FDR 
#' A get FDR table function
#' This function recieves a spectral library created be CreateLIBfromMGF function and adds FDR to it. FDR is calculated using the getFDRtable function.
#' @param library A spectral library
#' @keywords FDR
#' @examples lib1<-addFDRtoLibrary(lib1)
addFDRtoLibrary<-function(library)
{
  
  stopifnot(require(data.table))
  
  dt<-getFDRtable(library)
  
  l_mgf<-library[[2]]
  
  if(!("FDR" %in% names(l_mgf)))
    {
    l_mgf[,FDR:=NULL]
    l_mgf[,FullCount:=NULL]
    l_mgf[,AntiCount:=NULL]
    
    }
  
    
  
  n_mgf<-merge(l_mgf, dt, by='cluster_num', all.x=T)
 
  
  library[[2]]<-n_mgf
  return(library)
}### END Add FDR



#################
### Add IDs
#' Add IDs to library function
#' This function recieves a spectral library created be CreateLIBfromMGF function and a data.table with spectra identifications, and adds id's to the library. expenedLibIDs function is also called.
#' @param library A spectral library
#' @param id_dt A data.table of spectra IDs. SHould contain: Title, Sequence, Conf
#' @keywords IDs
#' @export
#' @examples lib1<-addIDsToLibrary(lib1)
addIDsToLibrary<-function(library, id_dt)
{
  
  stopifnot(require(data.table))

  l_mgf<-library[[2]]
  title<-id_dt$Title
  sequence<-id_dt$Sequence
  conf<-id_dt$Conf
  dt<-data.table(title=title, sequence=sequence, conf=conf)
  n_mgf<-merge(l_mgf,dt,by = 'title', all.x = T)
  
  library[[2]]<-n_mgf
  library<-expendLibIDs(library)
  return(library)
  
}### END Add IDs


#################
### Expend IDs
#' expend Lib IDs  function
#' This function recieves a spectral library with added Ids and expends them to the containing clusters.
#' @param library A spectral library
#' @keywords Expend
#' @examples lib1<-expendLibIDs(lib1)
expendLibIDs<-function(library)
{
  
  stopifnot(require(data.table))
  
  l_mgf<-library[[2]]
  if(!("sequence" %in% names(l_mgf)))return(library)# If there are no IDs, there's nothing to expend
  if(!("Expanded_sequence" %in% names(l_mgf))) l_mgf[,Expanded_sequence:=NULL]
  
  
  n_mgf<-merge(l_mgf, l_mgf[, .SD[which.max(conf), .(Expanded_sequence=sequence)], cluster_num], by = 'cluster_num', all = TRUE)
  
   
  library[[2]]<-n_mgf
  return(library)
  
}### END Expend IDs



#################
### Get Last mgf
#' get last MGF from Lib function
#' This function recieves a spectral library and returns an mgf data.table object of the mgf last added to the library.
#' @param library A spectral library
#' @keywords last
#' @export
#' @examples mgf1<-getLastMgfFromLib(lib1)
getLastMgfFromLib<-function(library)
{
  
  stopifnot(require(data.table))
 
return(library[[2]][title %in% library[[4]]])

  
}### END get last mgf



#################
### Get Modif
#' get Modif function
#' This function recieves clusters similarity table from the getClustersimilarity() function and returns a table of potential PTMs.
#' @param corr_dt A table of spectra similarity
#' @param D Iteration of modofication mass. Defaults to 0.01.
#' @keywords Modifications
#' @export
#' @examples modif1<-getModif(corr_dt1)
getModif<-function(corr_dt, D=0.01)
{
  
  stopifnot(require(data.table))
  
  max<-max(corr_dt$MassDiff)
  
  dt<-data.table(diff=numeric(), count=integer())
  
  i<-0
  while (i<max) {
    X<-nrow(corr_dt[abs(MassDiff)>i & abs(MassDiff)<=(i+D) & ChargeI==ChargeJ & Corr>C])
    if(X>0) dt<-rbindlist(list(dt,list(i,X)))
    i<-i+D
    
  }
  
  return(dt)
  
}### END get Modif




#################
### Get similarity 
#' get clusters similarity function
#' This function recieves a spectral library and returns a data.table of clusters pairs who's representative spectra are similar to each other.
#' @param library A spectral library
#' @param N Number of similar cluster pairs that need to be found. Defaults to 1000.
#' @param COR A minimal correlation between cluster representatives needed. Defaults to 50.
#' @param SIZE A minimal size of clusters considered. Defaults to 5.
#' @keywords similarity
#' @export
#' @examples corr_dt<-getClustersimilarity(lib1)
getClustersimilarity<-function(library, N=1000, COR=50, SIZE=5)
{
  
  stopifnot(require(data.table))
  
  cls<-library[[1]]
  
  nums<-rep(0,L)
  for(i in 1:L) nums[i]<-cls[[i]]$NUM
  big_nums<-which(nums>=SIZE)
  
  dt<-data.table(I=integer(), J=integer(), ChargeI=integer(), ChargeJ=integer(), MassDiff=numeric(), Corr=numeric())
  L<-length(cls)
  
  while (nrow(dt)<N) {
    i<-sample(big_nums,1)
    j<-sample(big_nums,1)
    
    I<-cls[i][[1]]
    J<-cls[j][[1]]
    diff<-(I$cons_prec-J$cons_prec)
    
    if(i!=j & I$cons_charge==J$cons_charge & abs(diff)>1)
    {
      list1<-I$cons_spec
      list2<-J$cons_spec
      cor<-Spec_list_corr(list1,list2)
      if(cor>=COR)
      {
        dt<-unique(rbindlist(list(dt,list(i,j,I$cons_charge,J$cons_charge,diff,cor))))
   
      }
    }
    
  }
  
  return(dt)
  
  
}### END get similarity 
