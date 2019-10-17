#' @export
#' @import pbapply


ffdGroupEvidenceFFBS <- function(ffdGroup, covariates, m0=0, Cova=100,
                                 delta = 0.95, S0 = 1, n0 = 1, N1 = FALSE, Nsimu1=100, Cutpos=30, r1, Test, mask, Ncores = NULL){
  
  if(N1==FALSE){N1 = dim(covariates)[1]}
  
  #TAKING THE POSITIONS FROM THE 4D IMAGE WITH NON-NULL VALUES 
  posiffd <- which(mask[,,] != 0, arr.ind = TRUE)
  
  
  
  
  ffd.out <- pbapply::pbapply(posiffd, 1, ffdGroupVoxelFFBS, ffdGroup, covariates, m0, Cova,
                              delta, S0, n0, N1, Nsimu1, r1, Test, Cutpos, cl = Ncores)
  
  #number of tests from the output of ffdIndividualVoxelFFBS  (Joint, marginal and LTT)
  Ntest <- 3
  vol.evidence <- list()
  
  
  for(k in 1:(Ntest)){
    vol.evidence[[k]] <- array(0, c(dim(covariates)[2], dim(mask)[1:3]))
  }
  
  
  for(ii in 1:dim(posiffd)[1]){
    vol.evidence[[1]][ ,posiffd[ii,1], posiffd[ii,2], posiffd[ii,3]] <- ffd.out[[ii]]$EvidenceJoint
    vol.evidence[[2]][ ,posiffd[ii,1], posiffd[ii,2], posiffd[ii,3]] <- ffd.out[[ii]]$EvidenceMargin
    vol.evidence[[3]][ ,posiffd[ii,1], posiffd[ii,2], posiffd[ii,3]] <- ffd.out[[ii]]$EvidenLTT
  }
  
  return(vol.evidence)
  
  
  
}











