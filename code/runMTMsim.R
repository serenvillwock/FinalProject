
## function to simulate data on two traits with different heritabilities and different genetic correlations between traits

runMTMsim <- function(T1h2, T2h2, addgenCor) {


  ### Species and population parameters
  historicalNe <- 50
  nFounders <- historicalNe
  nChr <- 7
  segSites <- 1400
  nQTL <- 50
  nSNP <- 500


  ### Initiate simulation
  #Specify Trait means and variances

  founderHaps <- runMacs2(nInd=nFounders, nChr=nChr, segSites=segSites, Ne=historicalNe, histNe=NULL, histGen=NULL)
  SP <- SimParam$new(founderHaps)


  traitMeans <- c(T1=31, T2=16)
  # Specify the additive variance and correlation: 1 on the diagonal
  addVar <- c(T1=1.5, T2=1.7) #trait 1 and trait 2
  addCor <- matrix(c(1, addgenCor, addgenCor, 1), nrow=2) #additive correlation matrix
  # Specify the error correlation and calculate error covariance
  ##### #calculate error var using the trait heritability
  errVar <- c(T1=  ((addVar[1] - T1h2*addVar[1])/ T1h2)  , T2=  ((addVar[2] - T2h2*addVar[2])/ T2h2) )
  #####
  errCor <- 0.5; errCov <- errCor*prod(sqrt(errVar))
  errCov <- matrix(c(errVar[1], errCov, errCov, errVar[2]), nrow=2) #error covariance matrix

  SP$addTraitA(nQtlPerChr=nQTL, mean=traitMeans, var=addVar, corA=addCor)
  SP$addSnpChip(nSnpPerChr=nSNP)

  # Create a new population of founders
  founders <- newPop(founderHaps, simParam=SP)
  nProgeny=nFounders
  progenyPop <- randCross(founders, nCrosses=nProgeny, nProgeny=6)
  mtPhenos=setPheno(progenyPop,varE=errCov)
  #Creating a matrix with phenotypes and genotypes
  simdata=cbind(c(1:length(progenyPop@id)), mtPhenos@pheno,mtPhenos@gv)


  ## Calculate genomic relationship matrix

  source("./code/calcRelationshipMatrices.R")
  snpMat <- pullSnpGeno(progenyPop)
  snpRelMat <- calcGenomicRelationshipMatrix(snpMat)




  ##### Analyze data using sommer
  library(sommer)

  # create dataframe to pass to the mixed model software
  dataf=data.frame(mean=as.factor(rep(1,length(simdata[,1]))),variety=as.factor(simdata[,1]),phenotype1=as.double(simdata[,2]),phenotype2=as.double(simdata[,3]),trueGV1=as.double(simdata[,4]),trueGV2=as.double(simdata[,5]))

  # Adding a small value to the diagonal of the G matrix to ensure it can be inverted
  Gsnp=snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1]))
  #Gsnpi=solve(snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1])))

  # Map elements in the relationship matrix to the phenotypes
  rownames(Gsnp)=levels(dataf$variety)
  colnames(Gsnp)=levels(dataf$variety)

  # Instructing asreml if the G matrix is already inverted
  attr(Gsnp, "INVERSE")=FALSE

  #running multi-trait and univariate models in sommer
  #MTM

  MTM <- mmer(cbind(phenotype1, phenotype2) ~ 1, #give a vector of two traits and sommer knows you are running a multitrait model
              random= ~ vs(variety, Gtc=unsm(2), Gu=Gsnp) ,
              rcov= ~ vs(units, Gtc=unsm(2)), #Gtc = unstructured 2x2 matrix; we don't know what those covariances are or if there is a pattern
              data=dataf, verbose = FALSE)
  summary(MTM)

  #running predict function to get variety effects for trait 1 and trait 2
  aMT<-sommer::predict.mmer(MTM,classify=c("variety"))
  pMT=aMT$pvals
  psMT=pMT[order(pMT$trait,pMT$variety),]

  #Single trait models
  T1 <- mmer(phenotype1 ~ 1,
             random= ~ vs(variety, Gu=Gsnp) ,
             rcov= ~ units,
             data=dataf, verbose = FALSE)
  summary(T1)

  aT1<-sommer::predict.mmer(T1,classify=c("variety"))
  pT1=aT1$pvals
  psT1=pT1[order(pT1$variety),]

  T2 <- mmer(phenotype2 ~ 1,
             random= ~ vs(variety, Gu=Gsnp) ,
             rcov= ~ units,
             data=dataf, verbose = FALSE)
  summary(T2)

  aT2<-sommer::predict.mmer(T2,classify=c("variety"))
  pT2=aT2$pvals
  psT2=pT2[order(pT2$variety),]


  ### Return correlations with true genetic values
  return( c(   cor(psMT$predicted.value[1:300],dataf$trueGV1), #correlation for trait 1 predictions MTM
               cor(psMT$predicted.value[301:600],dataf$trueGV2), #correlation for trait 2 predictions MTM
               cor(psT1$predicted.value,dataf$trueGV1), #correlation for trait 1 predictions single trait
               cor(psT2$predicted.value,dataf$trueGV2) #correlation for trait 2 predictions single trait
               ) )

}
