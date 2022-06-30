
## Function to generate a founder population
generate_founders <- function(DM_TC_adcorr = -0.2,
                              DM_TC_ercorr = -0.1,
                                DM_errvar = 41.08,
                                TC_errvar = 4.35){

  ## Species trait and population parameters
  historicalNe <- 200
  nFounders <- historicalNe
  nChr <- 18
  segSites <- 4400 #ICGMC 2015
  nQTL <- 80
  nSNP <- 1000

  ## Specify trait means and variances
  founderHaps <- runMacs2(nInd=nFounders, nChr=nChr, segSites=segSites,
                          Ne=historicalNe, histNe=NULL, histGen=NULL)
  SP <<- SimParam$new(founderHaps)

  traitMeans <- c(DM=24.12, TC=6.53) #DM= dry matter; TC= total carotenoids/g fresh weight
  # Parkes et al. 2020

  # Specify the additive variance and correlation: 1 on the diagonal
  addVar <- c(DM=18.44, TC=1.36) #Parkes 2020
  addCor <- matrix(c(1, DM_TC_adcorr, DM_TC_adcorr, 1), nrow=2) #additive correlation matrix
  #default correlation: -0.2
  # Specify the error correlation and calculate error covariance

  errVar <- c(DM=DM_errvar, TC=TC_errvar)
  #default error variances: DM: 41.08, TC: 4.35 (Parkes 2020 dominance + error variation)
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
  simdata=cbind(c(1:length(progenyPop@id)), mtPhenos@pheno, mtPhenos@gv)

  return(list(mtPhenos, as.data.frame(simdata)))

}



## Function to implement phenotypic & genomic selection on simulated population with given trait parameters ##
# Input: `progeny` is a list output from the generate_founders() function
# with the first element containing the progeny of the founding AlphaSim population
# and the second element containing a dataframe with progeny's phenotypes and genotypes

implement_pheno_selection <- function(progeny, multitrait=TRUE, nCycles = 5){

  pop_data <- progeny[[1]]
  pheno_data <- progeny[[2]]
  colnames(pheno_data) <- c("variety", "DM_pheno", "TC_pheno", "DM_gv", "TC_gv")
  pheno_data$variety <- as.factor(pheno_data$variety)

  ## Generate a relationship matrix from a three-column pedigree
  source("./code/calcRelationshipMatrices.R")
  ped <- cbind(pop_data@id, pop_data@father, pop_data@mother)
  parents <- unique(c(pop_data@father, pop_data@mother))
  founders <- cbind(parents, rep(0,length(parents)), rep(0, length(parents)))
  pedigree <- rbind(founders, ped)
  pedigree_named <- convertNamesToRows(pedigree)
  CC_mat <- pedigreeToCCmatrix(pedigree_named)

  # Calculate A matrix
  if(i == 1){
    #remove unphenotyped founders from CC_mat for first cycle
    CC_mat_prog <- CC_mat[-c(1:(length(parents))), -c(1:(length(parents)))]

  } else {CC_mat_prog <- CC_mat} #otherwise leave it as is

  # Map elements in the relationship matrix to the phenotypes
  rownames(CC_mat_prog) <- levels(pheno_data$variety)
  colnames(CC_mat_prog) = levels(pheno_data$variety)

  A_mat <- CC_mat_prog * 2


  ## Estimate variance parameters with a linear model
  # using the A matrix as the covariance structure for variety

  if (multitrait ==TRUE){

    MTmodel <- mmer(cbind(DM_pheno, TC_pheno) ~ 1,
                   random= ~ vs(variety, Gtc=unsm(2), Gu=A_mat),
                   rcov= ~ vs(units, Gtc=unsm(2)),
                          data = pheno_data, verbose = TRUE)


  } else {

    DMmodel <- mmer(DM_pheno ~ 1,
                    random= ~ vs(variety, Gtc=unsm(2), Gu=CC_mat),
                    rcov= ~ vs(units, Gtc=unsm(2)),
                    data = pheno_data, verbose = TRUE)

    TCmodel <- mmer(TC_pheno ~ 1,
                    random= ~ vs(variety, Gtc=unsm(2), Gu=CC_mat),
                    rcov= ~ vs(units, Gtc=unsm(2)),
                    data = pheno_data, verbose = TRUE)
    }





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
