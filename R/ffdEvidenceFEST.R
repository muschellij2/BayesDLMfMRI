#' @export
#' @import pbapply

ffdEvidenceFEST = function(ffdc, covariates, m0=0, Cova=100,
         delta=0.95, S0=1, n0=1, N1=FALSE, Nsimu1 = 100, Cutpos1=30, r1 = 1, perVol = 0.10, Test = "LTT", Ncores = NULL){
  
  
  if(N1==FALSE){N1=dim(ffdc)[4]}
  
  #TAKING THE POSITIONS FROM THE 4D IMAGE WITH NON-NULL VALUES 
  posiffd1 <- which(ffdc[,,,1] != 0, arr.ind = TRUE)
  
  
  if(Test == "LTT"){
    #COMPUTING THE EVIDENCE FOR BRAIN ACTIVATION: VOXEL-WISE ANALYSIS
    ffd.out = pbapply::pbapply(posiffd1, 1, ffdIndividualVoxelFEST, covariates, ffdc, m0, Cova,
                    delta, S0, n0, N1, Nsimu1, Cutpos1, Min.vol = perVol*max(ffdc), r1, Test, cl = Ncores)
    
    
    vol.evidence <- list()
    
    
    for(k in 1:(dim(covariates)[2])){
      vol.evidence[[k]] <- array(0, dim(ffdc)[1:3])
    }
    
    
    for(j in 1:dim(covariates)[2]){
      for(ii in 1:dim(posiffd1)[1]){
        vol.evidence[[j]][posiffd1[ii,1], posiffd1[ii,2], posiffd1[ii,3]] <- ffd.out[j, ii]
      }
    }
    
    return(vol.evidence) 
  }
  if(Test == "JointTest"){
    #COMPUTING THE EVIDENCE FOR BRAIN ACTIVATION: VOXEL-WISE ANALYSIS
    ffd.out = pbapply::pbapply(posiffd1, 1, ffdIndividualVoxelFEST, covariates, ffdc, m0, Cova,
                    delta, S0, n0, N1, Nsimu1, Cutpos1, Min.vol = perVol*max(ffdc), r1, Test, cl = Ncores)
    #number of tests from the output of ffdIndividualVoxelLTT  (Joint and marginal)
    Ntest <- 2
    pAux <- dim(covariates)[2]
    vol.evidence <- list()
    
    
    for(k in 1:(dim(covariates)[2]*Ntest)){
      vol.evidence[[k]] <- array(0, dim(ffdc)[1:3])
    }
    
    
    for(j in 1:dim(covariates)[2]){
     for(ii in 1:dim(posiffd1)[1]){
       vol.evidence[[j]][posiffd1[ii,1], posiffd1[ii,2], posiffd1[ii,3]] <- ffd.out[[ii]]$EvidenceJoint[j]
       vol.evidence[[Ntest+j]][posiffd1[ii,1], posiffd1[ii,2], posiffd1[ii,3]] <- ffd.out[[ii]]$EvidenceMargin[j]
      }
    }
    
    return(vol.evidence)
  }
  
  
  
}

  

  
