
## Function to generate a founder population
generate_founders <- function(DM_TC_adcorr = -0.2,
                              DM_TC_ercorr = -0.1,
                                DM_errvar = 41.08,
                                TC_errvar = 4.35,
                              FounderNum = 200,
                              numProgeny = 6){

  ## Species trait and population parameters
  historicalNe <- 200
  nFounders <- FounderNum
  nChr <- 18
  segSites <- 4400 #ICGMC 2015
  nQTL <- 80
  nSNP <- 1000

  ## Specify trait means and variances
  founderHaps <- runMacs2(nInd=nFounders, nChr=nChr, segSites=segSites,
                          Ne=historicalNe, histNe=NULL, histGen=NULL)
  SP <- SimParam$new(founderHaps)

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
  saveRDS(errCov, "./data/errCov.RDS")

  SP$addTraitA(nQtlPerChr=nQTL, mean=traitMeans, var=addVar, corA=addCor)
  SP$addSnpChip(nSnpPerChr=nSNP)

  #save to global environment
  saveRDS(SP, file="./data/SP.RDS")
  SP <<- SP

  # Create a founder population
  founders <- newPop(founderHaps, simParam=SP)
  nProgeny <- numProgeny
  progenyPop <- randCross(founders, nCrosses=nFounders, nProgeny=nProgeny)
  mtPhenos = setPheno(progenyPop,varE=errCov)

  #Create a matrix with phenotypes and genotypes
  simdata=cbind(progenyPop@id, mtPhenos@pheno, mtPhenos@gv)

  return(list(mtPhenos, as.data.frame(simdata)))

}



## Function to implement phenotypic selection on simulated population with given trait parameters ##
# Input: `progeny` is a list output from the generate_founders() function
# with the first element containing the progeny of the founding AlphaSim population
# and the second element containing a dataframe with progeny's phenotypes and genotypes

implement_pheno_selection <- function(progeny, multitrait=TRUE, nCycles = 5,
                                      selection_proportion=0.1, nProgeny=6, nFounders=200){

  pop_data <- progeny[[1]]
  pheno_data <- progeny[[2]]


  colnames(pheno_data) <- c("variety", "DM_pheno", "TC_pheno", "DM_gv", "TC_gv")
  pheno_data$variety <- as.factor(pheno_data$variety)

  #DM_pheno <- as.numeric(pheno_data$DM_pheno)
  #TC_pheno <- as.numeric(pheno_data$TC_pheno)
  pheno_data$DM_pheno <- as.numeric(pheno_data$DM_pheno)
  pheno_data$TC_pheno <- as.numeric(pheno_data$TC_pheno)

  # set up pedigree
  ped_original <- cbind(pop_data@id, pop_data@father, pop_data@mother)
  parents <- unique(c(pop_data@father, pop_data@mother))
  founders <- cbind(parents, rep(0,length(parents)), rep(0, length(parents)))
  pedigree <- rbind(founders, ped_original)


  # Set up a dataframe to store data from each cycle
  traceP <- data.frame()

  # Save initial phenotypes
  to_save <- data.frame(cbind(pheno_data, cycle = rep(i, nrow(pheno_data))))
  traceP <- rbind(traceP, to_save)

  for(i in (1:nCycles)){
    print(i)

    ## Generate a relationship matrix from a three-column pedigree
    source("./code/calcRelationshipMatrices.R")
    if(i != 1){ #pedigree was already set up above for first cycle
      ped <- cbind(pop_data@id, pop_data@father, pop_data@mother)
      pedigree <- rbind(pedigree, ped)}

    pedigree_named <- convertNamesToRows(pedigree)
    CC_mat <- pedigreeToCCmatrix(pedigree_named)


    # Subset the CC_matrix for just this generation by cutting off the parental rows and columns
    CC_end <- dim(CC_mat)[1]
    CC_start <- CC_end - (nFounders*nProgeny) + 1
    CC_mat_prog <- CC_mat[CC_start:CC_end, CC_start:CC_end]


    # Map elements in the relationship matrix to the phenotypes
    rownames(CC_mat_prog) <- pheno_data$variety
    colnames(CC_mat_prog) <- pheno_data$variety


    # Calculate A matrix (additive relationship = twice coefficient of coancestry)
    A_mat <- CC_mat_prog * 2
    Amatfilename <- paste0("../FinalProject_RDS_backup/Amat_cycle", i, ".RDS")
    saveRDS(A_mat, file = Amatfilename)
    print("A matrix calculated")


    ## Estimate variance parameters with a linear model
    # using the A matrix as the covariance structure for variety
    if (multitrait == TRUE){

      # Fit model
      MTmodel <- mmer(cbind(DM_pheno, TC_pheno) ~ 1,
                     random= ~ vsr(variety, Gtc=unsm(2), Gu=A_mat),
                     rcov= ~ vsr(units, Gtc=unsm(2)),
                            data = pheno_data, verbose = TRUE)

      MTfilename <- paste0("data/MTmodel_cycle", i, ".RDS")
      saveRDS(MTmodel, file = MTfilename)
      print("model fit done")


      # Predict variety effects
      aMT <- sommer::predict.mmer(MTmodel, classify=c("variety"))
      pMT = aMT$pvals
      psMT = pMT[order(pMT$trait, pMT$variety),]

      # save checkpoint
      psMTfilename <- paste0("../FinalProject_RDS_backup/psMT_cycle", i, ".RDS")
      saveRDS(psMT, file = psMTfilename)
      print("predictions done")

      # reformat predictions
      DM_preds <- psMT[psMT$trait == "DM_pheno",]
      TC_preds <- psMT[psMT$trait == "TC_pheno",]
      all_preds <- merge(DM_preds, TC_preds, by="variety")


      # Extract variance/covariance estimates
      addCov_est <- MTmodel$sigma$`u:variety`
      errCov_est <- MTmodel$sigma$`u:units`
      phenCov_est <-  addCov_est + errCov_est


      # Calculate selection index
      econVal <- c(DM = 1, TC = 0.8)
      selIdxCoef <- solve(phenCov_est) %*% addCov_est %*% econVal  #solve inverts the matrix
      all_preds$selIdx <- as.matrix(all_preds[,c(3,6)])  %*%  selIdxCoef


      # Make selections #it is selecting some parents / maybe variety names are off
      nToSelect <- nFounders
      selected <- all_preds[order(all_preds$selIdx, decreasing=T),"variety"][1:nToSelect]
      MT_Pselected_varieties <- pop_data[as.character(selected)]
      print("selections done")

      # Make crosses, unless it's the last cycle
      if(i != nCycles){

        new_pop <- randCross(MT_Pselected_varieties, nCrosses=nFounders, nProgeny=nProgeny)
        new_pop_phenod <- setPheno(new_pop, varE=errCov, simParam=SP)

        # Set up matrix for next cycle
        pop_data <- new_pop
        pheno_data <- as.data.frame(cbind(new_pop_phenod@id, new_pop_phenod@pheno, new_pop_phenod@gv))

        colnames(pheno_data) <- c("variety", "DM_pheno", "TC_pheno", "DM_gv", "TC_gv")
        pheno_data$variety <- as.factor(pheno_data$variety)
        pheno_data$DM_pheno <- as.numeric(pheno_data$DM_pheno)
        pheno_data$TC_pheno <- as.numeric(pheno_data$TC_pheno)

        #save data
        to_save <- data.frame(cbind(pheno_data, cycle = rep(i, nrow(pheno_data))))
        traceP <- rbind(traceP, to_save)
        print(paste0("cycle ", i, " saved!"))
        saveRDS(traceP, file="../FinalProject_RDS_backup/lasttraceP.RDS")

      } # end crosses
  } #end multitrait model predictions

    else { #single trait model:

      # set up pedigree
      ped_original <- cbind(pop_data@id, pop_data@father, pop_data@mother)
      parents <- unique(c(pop_data@father, pop_data@mother))
      founders <- cbind(parents, rep(0,length(parents)), rep(0, length(parents)))
      pedigree <- rbind(founders, ped_original)

      pedigree_named <- convertNamesToRows(pedigree)
      CC_mat <- pedigreeToCCmatrix(pedigree_named)

      # Subset the CC_matix for just this generation
      CC_end <- dim(CC_mat)[1]
      CC_start <- CC_end - (nFounders*nProgeny) + 1
      CC_mat_prog <- CC_mat[CC_start:CC_end, CC_start:CC_end]

      # Map elements in the relationship matrix to the phenotypes
      rownames(CC_mat_prog) <- levels(pheno_data$variety)
      colnames(CC_mat_prog) <- levels(pheno_data$variety)


      # Calculate A matrix (additive relationship = twice coefficient of coancestry)
      A_mat <- CC_mat_prog * 2
      Amatfilename <- paste0("../FinalProject_RDS_backup/Amat_ST_cycle", i, ".RDS")
      saveRDS(A_mat, file = Amatfilename)
      print("A matrix calculated")

      DMmodel <- mmer(DM_pheno ~ 1,
                      random= ~ vsr(variety, Gu=A_mat),
                      rcov= ~ units,
                      data = pheno_data, verbose = TRUE)

      DMfilename <- paste0("data/DMmodel_cycle", i, ".RDS")
      saveRDS(DMmodel, file = DMfilename)
      print("DM model fit")


      TCmodel <- mmer(TC_pheno ~ 1,
                      random= ~ vsr(variety, Gu=A_mat),
                      rcov= ~ units,
                      data = pheno_data, verbose = TRUE)

      TCfilename <- paste0("data/TCmodel_cycle", i, ".RDS")
      saveRDS(TCmodel, file = TCfilename)
      print("TC model fit")


      # Predict variety effects for DM
      aDM <- sommer::predict.mmer(DMmodel, classify=c("variety"))
      pDM = aDM$pvals
      psDM = pDM[order(pDM$trait,pDM$variety),]

      # save checkpoint
      psDMfilename <- paste0("../FinalProject_RDS_backup/psDM_cycle", i, ".RDS")
      saveRDS(psDM, file = psDMfilename)
      print("DM predictions done")

      # Predict variety effects for TC
      aTC <- sommer::predict.mmer(TCmodel, classify=c("variety"))
      pTC = aTC$pvals
      psTC = pTC[order(pTC$trait,pTC$variety),]

      # save checkpoint
      psTCfilename <- paste0("../FinalProject_RDS_backup/psTC_cycle", i, ".RDS")
      saveRDS(psTC, file = psTCfilename)
      print("TC predictions done")

      # merge predictions
      #DM_preds <- psDM[psMT$trait == "DM_pheno",]
      #TC_preds <- psMT[psMT$trait == "TC_pheno",]
      all_preds <- merge(psDM, psTC, by="variety")


      ## Estimate variance/covariance parameters:
      # (Additive covariance = cov of the breeding values?)
      addCov_est <- cov(cbind(all_preds$predicted.value.x, all_preds$predicted.value.y))
      # (Error covariance = cov of the residuals?)
      errCov_est <- cov(cbind(DMmodel$residuals, TCmodel$residuals))
      phenCov_est <-  addCov_est + errCov_est


      # Calculate selection index
      econVal <- c(DM = 1, TC = 0.8)
      selIdxCoef <- solve(phenCov_est) %*% addCov_est %*% econVal  #solve inverts the matrix
      all_preds$selIdx <- as.matrix(all_preds[,c(3,6)])  %*%  selIdxCoef


      # Make selections
      nToSelect <- nrow(all_preds)*selection_proportion #select top 10%
      selected <- all_preds[order(all_preds$selIdx, decreasing=T),"variety"][1:nToSelect]
      ST_Pselected_varieties <- pop_data[as.numeric(selected)]
      print("selections done")

    # Make crosses, unless it's the last cycle
      if(i != nCycles){

        new_pop <- randCross(ST_Pselected_varieties, nCrosses=nFounders, nProgeny=nProgeny)
        new_pop_phenod = setPheno(new_pop, varE=errCov)

        # Set up matrix for next cycle
        pop_data <- new_pop
        pheno_data <- as.data.frame(cbind(new_pop_phenod@id, new_pop_phenod@pheno, new_pop_phenod@gv))

        colnames(pheno_data) <- c("variety", "DM_pheno", "TC_pheno", "DM_gv", "TC_gv")
        pheno_data$variety <- as.factor(pheno_data$variety)
        pheno_data$DM_pheno <- as.numeric(pheno_data$DM_pheno)
        pheno_data$TC_pheno <- as.numeric(pheno_data$TC_pheno)

        #save data
        to_save <- data.frame(cbind(pheno_data, cycle = rep(i, nrow(pheno_data))))
        traceP <- rbind(traceP, to_save)
        saveRDS(traceP, file="../FinalProject_RDS_backup/lasttraceP.RDS")
        print(paste0("cycle ", i, " saved!"))

      } #end crosses
    } #end single trait predictions
  } #end cycle

  return(traceP)

} #end function








## Function to implement genomic selection on simulated population with given trait parameters ##
# Input: `progeny` is a list output from the generate_founders() function
# with the first element containing the progeny of the founding AlphaSim population
# and the second element containing a dataframe with progeny's phenotypes and genotypes

implement_geno_selection <- function(progeny, multitrait=TRUE, nCycles = 5,
                                      selection_proportion=0.1, nProgeny=6,
                                     nFounders=200){

  pop_data <- progeny[[1]]
  pheno_data <- progeny[[2]]

  colnames(pheno_data) <- c("variety", "DM_pheno", "TC_pheno", "DM_gv", "TC_gv")
  pheno_data$variety <- as.factor(pheno_data$variety)

  # Set up a dataframe to store data from each cycle
  traceG <- data.frame()

  # Save initial population phenotypes
  to_save <- data.frame(cbind(pheno_data, cycle = rep(i, nrow(pheno_data))))
  traceG <- rbind(traceG, to_save)


  for (i in 1:nCycles){

    ## Calculate genomic relationship matrix
    source("./code/calcRelationshipMatrices.R")
    snpMat <- pullSnpGeno(pop_data)
    snpRelMat <- calcGenomicRelationshipMatrix(snpMat)

    if (multitrait == TRUE){
      print(paste0("starting cycle ", i)) #checkpoint

      # Calculate genomic relationship matrix
      snpMat <- pullSnpGeno(pop_data)
      snpRelMat <- calcGenomicRelationshipMatrix(snpMat)
      print("genomic relationship matrix calculated")

      # Adding a small value to the diagonal of the G matrix to ensure it can be inverted
      Gsnp=snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1]))
      #Gsnpi=solve(snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1])))

      # Map elements in the relationship matrix to the phenotypes
      rownames(Gsnp)=levels(pheno_data$variety)
      colnames(Gsnp)=levels(pheno_data$variety)


      # Fit multi-trait model

      MTM <- mmer(cbind(DM_phenotype, TC_phenotype) ~ 1,
                  #give a vector of two traits and sommer knows you are running a multitrait model
                  random= ~ vs(variety, Gtc=unsm(2), Gu=Gsnp) ,
                  rcov= ~ vs(units, Gtc=unsm(2)), #Gtc = unstructured 2x2 matrix
                  data=pheno_data, verbose = TRUE)
      print("model fit done")

      # Extract BLUPs
      DM_BLUPs <- unlist(randef(MTM)[[1]][1])
      TC_BLUPs <- unlist(randef(MTM)[[1]][2])
      MTG_BLUPs <- cbind(DM_BLUPs, TC_BLUPs)

      # Extract variance/covariance estimates:
      addCov_est <- MTM$sigma #fix this
      errCov_est <- MTM$sigma #fix this
      phenCov_est <-  addCov_est + errCov_est

      # Calculate selection index
      econVal <- c(DM = 1, TC = 0.8)
      selIdxCoef <- solve(phenCov_est) %*% addCov_est %*% econVal  #solve inverts the matrix
      all_preds$selIdx <- as.matrix(all_preds[,c(3,6)])  %*%  selIdxCoef


      # Make selections
      nToSelect <- nrow(all_preds)*selection_proportion #select top 10%
      selected <- all_preds[order(all_preds$selIdx, decreasing=T),"variety"][1:nToSelect]
      ST_Pselected_varieties <- pop_data[as.numeric(selected)]
      print("selections done")

      # Select best lines
      MTGselected_varieties <- newProgenyPop[MTGwhichToSelect]

      # Make random crosses among selected lines
      pop_data <- randCross(MTGselected_varieties, nCrosses=nProgeny, nProgeny=6)

      # Phenotype the progeny
      pop_data_phenod = setPheno(pop_data, varE=errCov, simParam = SP)

      pheno_data <- as.data.frame(cbind(as.factor(pop_data_phenod@id),
                                        pop_data_phenod@pheno,
                                        pop_data_phenod@gv))
      colnames(pheno_data) <- c("variety", "DM_pheno", "TC_pheno", "DM_gv", "TC_gv")

      # Save phenotypes of selected lines
      to_save <- data.frame(pheno_data,
                            cycle = rep(i, length(pheno_data)))
      traceG <- rbind(traceG, to_save)

      } # end multi-trait genomic selection

    else{

      print(paste0("ST genomic selection cycle", i)) #verbose checkpoint

      # Calculate genomic relationship matrix
      snpMat <- pullSnpGeno(pop_data)
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
      pop_data <- randCross(STGselected_varieties, nCrosses=nProgeny, nProgeny=6)

      # Phenotype the progeny
      pop_data_phenod = setPheno(pop_data, varE=errCov, simParam = SP)

      pheno_data <- as.data.frame(cbind(as.factor(pop_data_phenod@id),
                                        pop_data_phenod@pheno,
                                        pop_data_phenod@gv))
      colnames(pheno_data) <- c("variety", "DM_pheno", "TC_pheno", "DM_gv", "TC_gv")

      # Save phenotypes of selected lines
      to_save <- data.frame(pheno_data,
                            cycle = rep(i, length(pheno_data)))
      traceG <- rbind(traceG, to_save)

    } #end single-trait genomic selection

  } #end cycle

  return(traceG)

} #end function








