
## Function to implement phenotypic & genomic selection on simulated population with given trait parameters ##

implement_selection <- function(DM_TC_adcorr = 0.2,
                        DM_errvar = 3,
                        TC_errvar = 0.1,
                        DM_TC_ercorr = 0.2,
                        nCycles = 5){

  ## Species trait and population parameters
  historicalNe <- 200
  nFounders <- historicalNe
  nChr <- 18
  segSites <- 4400 #ICGMC 2015
  nQTL <- 80
  nSNP <- 1000

  ## Specify trait means and variances
  founderHaps <- runMacs2(nInd=nFounders, nChr=nChr, segSites=segSites, Ne=historicalNe, histNe=NULL, histGen=NULL)
  SP <- SimParam$new(founderHaps)

  traitMeans <- c(DM=28.6, TC=3.8) #DM= dry matter; TC= total carotenoids/g fresh weight

  # Specify the additive variance and correlation: 1 on the diagonal
  addVar <- c(DM=16.52, TC=1)
  addCor <- matrix(c(1, DM_TC_adcorr, DM_TC_adcorr, 1), nrow=2) #additive correlation matrix
  # Specify the error correlation and calculate error covariance
  errVar <- c(DM=DM_errvar, TC=TC_errvar)
  errCor <- DM_TC_ercorr; errCov <- errCor*prod(sqrt(errVar))
  errCov <- matrix(c(errVar[1], errCov, errCov, errVar[2]), nrow=2) #error covariance matrix

  SP$addTraitA(nQtlPerChr=nQTL, mean=traitMeans, var=addVar, corA=addCor)
  SP$addSnpChip(nSnpPerChr=nSNP)

  # Create a founder population
  founders <- newPop(founderHaps, simParam=SP)
  nProgeny = nFounders
  progenyPop <- randCross(founders, nCrosses=nProgeny, nProgeny=6)
  mtPhenos = setPheno(progenyPop,varE=errCov)

  #Create a matrix with phenotypes and genotypes
  simdata=cbind(c(1:length(progenyPop@id)), mtPhenos@pheno,mtPhenos@gv)


  #Calculate selection index
  econVal <- c(DM = 1, TC = 0.8)
  addCov <- addCor * sqrt( addVar[1] * addVar[2]  )
  phenCov <-  addCov + errCov
  selIdxCoef <- solve(phenCov) %*% addCov %*% econVal  #solve inverts the matrix
  nToSelect <- nProgeny*0.1 #select top 10%




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
