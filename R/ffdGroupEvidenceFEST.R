#' @export
#' @import pbapply


ffdGroupEvidenceFEST <- function(ffdGroup, covariates, m0=0, Cova=100,
                                delta = 0.95, S0 = 1, n0 = 1, N1 = FALSE, Nsimu1=100, Cutpos=30, r1, Test, mask, Ncores = NULL){
  
  if(N1==FALSE){N1 = dim(covariates)[1]}
  
  #TAKING THE POSITIONS FROM THE 4D IMAGE WITH NON-NULL VALUES 
  posiffd <- which(mask[,,] != 0, arr.ind = TRUE)
  
  
  
  if(Test == "LTT"){
    ffd.out <- pbapply::pbapply(posiffd, 1, ffdGroupVoxelFEST, ffdGroup, covariates, m0, Cova,
                                delta, S0, n0, N1, Nsimu1, r1, Test, Cutpos, cl = Ncores)
    
    vol.evidence <- list()
    
    for(k in 1:(dim(covariates)[2])){
      vol.evidence[[k]] <- array(0, dim(mask))
    }
    
    
    for(j in 1:(dim(covariates)[2])){
      for(i in 1:dim(posiffd)[1]){
        vol.evidence[[j]][posiffd[i,1], posiffd[i,2], posiffd[i,3]] <- ffd.out[j, i]
      }
    }
    
    return(vol.evidence)
    
  }
  if(Test == "Joint"){
    
    ffd.out = pbapply::pbapply(posiffd, 1, ffdGroupVoxelFEST, ffdGroup, covariates, m0, Cova,
                               delta, S0, n0, N1, Nsimu1, r1, Test, Cutpos, cl = Ncores)
    
    Ntest <- 2
    pAux <- dim(covariates)[2]
    vol.evidence <- list()
    
    
    for(k in 1:(dim(covariates)[2]*Ntest)){
      vol.evidence[[k]] <- array(0, dim(mask)[1:3])
    }
    
    
    for(j in 1:dim(covariates)[2]){
      for(ii in 1:dim(posiffd)[1]){
        vol.evidence[[j]][posiffd[ii,1], posiffd[ii,2], posiffd[ii,3]] <- ffd.out[[ii]]$EvidenceJoint[j]
        vol.evidence[[Ntest+j]][posiffd[ii,1], posiffd[ii,2], posiffd[ii,3]] <- ffd.out[[ii]]$EvidenceMargin[j]
      }
    }
    
    
    return(vol.evidence)
    
  }
  
  
}
    
  









