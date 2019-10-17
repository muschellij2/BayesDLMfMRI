#' @export
#' @import pbapply

ffdEvidenceFFBS = function(ffdc, covariates, m0=0, Cova=100,
                           delta=0.95, S0=1, n0=1, N1=FALSE, Nsimu1 = 100, Cutpos1=30, r1 = 1, perVol = 0.10, Ncores = NULL){
  
  
  if(N1==FALSE){N1=dim(ffdc)[4]}
  
  #TAKING THE POSITIONS FROM THE 4D IMAGE WITH NON-NULL VALUES 
  posiffd1 <- which(ffdc[,,,1] != 0, arr.ind = TRUE)
  
  

    #COMPUTING THE EVIDENCE FOR BRAIN ACTIVATION: VOXEL-WISE ANALYSIS
    ffd.out = pbapply::pbapply(posiffd1, 1, ffdIndividualVoxelFFBS, covariates, ffdc, m0, Cova,
                    delta, S0, n0, N1, Nsimu1, Cutpos1, Min.vol = perVol*max(ffdc), r1, cl = Ncores)
    #number of tests from the output of ffdIndividualVoxelFFBS  (Joint, marginal and LTT)
    Ntest <- 3
    vol.evidence <- list()
    
    
    for(k in 1:(Ntest)){
      vol.evidence[[k]] <- array(0, c(dim(covariates)[2], dim(ffdc)[1:3]))
    }
    
    
      for(ii in 1:dim(posiffd1)[1]){
        vol.evidence[[1]][ ,posiffd1[ii,1], posiffd1[ii,2], posiffd1[ii,3]] <- ffd.out[[ii]]$EvidenceJoint
        vol.evidence[[2]][ ,posiffd1[ii,1], posiffd1[ii,2], posiffd1[ii,3]] <- ffd.out[[ii]]$EvidenceMargin
        vol.evidence[[3]][ ,posiffd1[ii,1], posiffd1[ii,2], posiffd1[ii,3]] <- ffd.out[[ii]]$EvidenLTT
      }
    
    return(vol.evidence)
  
  
  
}



