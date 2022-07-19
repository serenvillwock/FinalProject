


## Implement phenotypic selection on simulated population:
Pphenotypes <- AlphaSimR::pheno(mtPhenos)
newPhenosP <- mtPhenos
traceP <- data.frame()

for (i in 1:nCycles){

  print(paste0("phenotypic selection cycle", i)) #verbose checkpoint

  # Calculate selection index
  indexP <- Pphenotypes %*% selIdxCoef

  # Select best lines
  whichToSelectP <- order(indexP, decreasing=T)[1:nToSelect]
  Pselected_varieties <- newPhenosP[whichToSelectP]

  # Make random crosses among selected lines
  newPopP <- randCross(Pselected_varieties, nCrosses=nProgeny, nProgeny=6)

  # Phenotype the progeny
  newPhenosP = setPheno(newPopP, varE=errCov, simParam = SP)
  Pphenotypes <- AlphaSimR::pheno(newPhenosP)


  # Save phenotypes
  to_save <- data.frame(Pphenotypes,
                        cycle = rep(i, length(whichToSelectP)))
  traceP <- rbind(traceP, to_save)

} # end phenotypic selection




## Implement multitrait genomic selection on simulated population:
newProgenyPop <- mtPhenos
MTGphenotypes <- AlphaSimR::pheno(mtPhenos)
traceMTG <- data.frame()

for (i in 1:nCycles){

  print(paste0("MT genomic selection cycle", i)) #verbose checkpoint

  ## Calculate genomic relationship matrix
  source("./code/calcRelationshipMatrices.R")
  snpMat <- pullSnpGeno(newProgenyPop)
  snpRelMat <- calcGenomicRelationshipMatrix(snpMat)


  ## Model the data with sommer, predict best lines
  # create dataframe to pass to the mixed model software
  dataf=data.frame(mean=as.factor(rep(1,length(simdata[,1]))),variety=as.factor(simdata[,1]),
                   DM_phenotype=as.double(simdata[,2]),TC_phenotype=as.double(simdata[,3]),
                   trueGV1=as.double(simdata[,4]),trueGV2=as.double(simdata[,5]))

  # Adding a small value to the diagonal of the G matrix to ensure it can be inverted
  Gsnp=snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1]))
  #Gsnpi=solve(snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1])))

  # Map elements in the relationship matrix to the phenotypes
  rownames(Gsnp)=levels(dataf$variety)
  colnames(Gsnp)=levels(dataf$variety)


  # Fit multi-trait model

  MTM <- mmer(cbind(DM_phenotype, TC_phenotype) ~ 1,
              #give a vector of two traits and sommer knows you are running a multitrait model
              random= ~ vs(variety, Gtc=unsm(2), Gu=Gsnp) ,
              rcov= ~ vs(units, Gtc=unsm(2)), #Gtc = unstructured 2x2 matrix
              data=dataf, verbose = TRUE)

  # Extract BLUPs
  DM_BLUPs <- unlist(randef(MTM)[[1]][1])
  TC_BLUPs <- unlist(randef(MTM)[[1]][2])
  MTG_BLUPs <- cbind(DM_BLUPs, TC_BLUPs)

  # Calculate selection index for multitrait model:
  MTGindex <- as.matrix(MTG_BLUPs) %*% selIdxCoef
  MTGwhichToSelect <- order(MTGindex, decreasing=T)[1:nToSelect]

  # Select best lines
  MTGselected_varieties <- newProgenyPop[MTGwhichToSelect]

  # Make random crosses among selected lines
  newProgenyPop <- randCross(MTGselected_varieties, nCrosses=nProgeny, nProgeny=6)

  # Phenotype the progeny
  newProgenyPopPheno = setPheno(newProgenyPop, varE=errCov, simParam = SP)

  MTGphenotypes <- AlphaSimR::pheno(newProgenyPopPheno)

  # Save phenotypes of selected lines
  to_save <- data.frame(MTGphenotypes,
                        cycle = rep(i, length(MTGwhichToSelect)))
  traceMTG <- rbind(traceMTG, to_save)

} # end multi-trait genomic selection





## Implement single trait genomic selection on simulated population:
newProgenyPop <- mtPhenos
STGphenotypes <- AlphaSimR::pheno(mtPhenos)
traceSTG <- data.frame()

for (i in 1:nCycles){

  print(paste0("ST genomic selection cycle", i)) #verbose checkpoint

  ## Calculate genomic relationship matrix
  source("./code/calcRelationshipMatrices.R")
  snpMat <- pullSnpGeno(newProgenyPop)
  snpRelMat <- calcGenomicRelationshipMatrix(snpMat)


  ## Model the data with sommer, predict best lines
  # create dataframe to pass to the mixed model software
  dataf=data.frame(mean=as.factor(rep(1,length(simdata[,1]))),variety=as.factor(simdata[,1]),
                   DM_phenotype=as.double(simdata[,2]),TC_phenotype=as.double(simdata[,3]),
                   trueGV1=as.double(simdata[,4]),trueGV2=as.double(simdata[,5]))

  # Adding a small value to the diagonal of the G matrix to ensure it can be inverted
  Gsnp=snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1]))
  #Gsnpi=solve(snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1])))

  # Map elements in the relationship matrix to the phenotypes
  rownames(Gsnp)=levels(dataf$variety)
  colnames(Gsnp)=levels(dataf$variety)


  #Single trait models
  DM_S <- mmer(DM_phenotype ~ 1,
               random= ~ vs(variety, Gu=Gsnp),
               rcov= ~ units,
               data=dataf, verbose = TRUE)

  TC_S <- mmer(TC_phenotype ~ 1,
               random= ~ vs(variety, Gu=Gsnp),
               rcov= ~ units,
               data=dataf, verbose = TRUE)


  # Extract BLUPs
  DM_BLUPs <- unlist(randef(DM_S)[[1]])
  TC_BLUPs <- unlist(randef(TC_S)[[1]])
  STG_BLUPs <- cbind(DM_BLUPs, TC_BLUPs)

  # Calculate selection index for multitrait model:
  STGindex <- as.matrix(STG_BLUPs) %*% selIdxCoef
  STGwhichToSelect <- order(STGindex, decreasing=T)[1:nToSelect]

  # Select best lines
  STGselected_varieties <- newProgenyPop[STGwhichToSelect]

  # Make random crosses among selected lines
  newProgenyPop <- randCross(STGselected_varieties, nCrosses=nProgeny, nProgeny=6)

  # Phenotype the progeny
  newProgenyPopPheno = setPheno(newProgenyPop, varE=errCov, simParam = SP)

  STGphenotypes <- AlphaSimR::pheno(newProgenyPopPheno)

  # Save phenotypes of selected lines
  to_save <- data.frame(STGphenotypes,
                        cycle = rep(i, length(STGwhichToSelect)))
  traceSTG <- rbind(traceSTG, to_save)

} # end single-trait genomic selection


return(list(traceP, traceMTG, traceSTG))


} # end implement_selection function
