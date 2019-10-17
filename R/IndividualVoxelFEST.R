#' @export


IndividualVoxelFEST <- function(posi.ffd, covariates, ffdc, m0, Cova, delta, S0, n0, N1, Nsimu1, Cutpos1, Min.vol, r1, Test){
  
  
  if(r1 == 0){
    
    posi <- distanceNeighbors(posi.refer = as.vector(posi.ffd), r1)
    
    #BOLD RESPONSE SERIES IN THE CLUSTER RELATED TO posi
    series.def <- sapply(1:(dim(posi)[1]), function(k){ffdc[posi[k,1], posi[k,2], posi[k,3], ]})
    #CHEKING THE THRESHOLD Min.vol FOR THE MAIN TS: JUST TEMPORAL SERIES ABOVE THE THRESHOLD, DISCARD TEMPORAL SERIES WITH NON-SIGNIFICANT SIGNAL
    
    if(min(series.def[,1]) < Min.vol){if(Test=="LTT"){return( rep(NA, dim(covariates)[2]) )}
      if(Test=="JointTest"){return( list(EvidenceJoint = rep(NA, dim(covariates)[2]), 
                                         EvidenceMargin = rep(NA, dim(covariates)[2]) ) )}}else{
                                           
                                           series.def <- matrix((series.def - mean(series.def))/sd(series.def), ncol=1)                                     
                                           m01   <- matrix(rep(m0, dim(covariates)[2]*dim(series.def)[2]), ncol=1)
                                           Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                           S01 <- diag(rep(S0,dim(series.def)[2]))
                                           #DISCOUNT FACTORS MATRIX
                                           delta1 <- sqrt(delta)
                                           Beta1 <-  diag(1/c(rep(delta1, dim(covariates)[2])))
                                           
                                           if(Test=="LTT"){
                                             return(NA)
                                           }
                                           
                                           if(Test=="JointTest"){
                                             res <- .Individual_FunctionalMultiTest(ffd1 = as.matrix(series.def), Cova = as.matrix(covariates), m0In = m01, c0In = Cova1, 
                                                                                    S0In = S01, beta0In = Beta1, nt0In = n0, NIn = N1, Nsimu = Nsimu1, CUTpos = Cutpos1)
                                             
                                             return(list(EvidenceJoint = rep(NA, dim(covariates)[2]), EvidenceMargin = res$EvidenMarginal))
                                             
                                             
                                           }
                                           
                                           
                                         }
  }else{
    
    #THIS LINE RETURN THE POSITIONS OF EACH VOXEL INSIDE THE CLUSTER GIVEN THE DISTANCE r1
    posi1 <- distanceNeighbors(posi.refer = as.vector(posi.ffd), r1)
    
    
    aux.pos <- dim(ffdc)[1:3]
    #GOING THROUGH EACH ROW AND CHECKING IF ANY POSITION IS OUTSIDE THE BOUNDS
    row_sub1 = apply(posi1, 1, function(row, x1){0 < row & row<=x1}, x1=aux.pos)
    posi <- posi1[apply(t(row_sub1), 1, sum)==3, ]
    #BOLD RESPONSE SERIES FOR THE CLUSTER RELATED TO posi
    series <- sapply(1:(dim(posi)[1]), function(k){ffdc[posi[k,1], posi[k,2], posi[k,3], ]})
    #CHEKING THE THRESHOLD Min.vol FOR THE MAIN TS: JUST TEMPORAL SERIES ABOVE THE THRESHOLD, DISCARD TEMPORAL SERIES WITH NON-SIGNIFICANT SIGNAL
    if(min(series[,1])<Min.vol){if(Test=="LTT"){return( rep(NA, dim(covariates)[2]) )}
      if(Test=="JointTest"){return( list(EvidenceJoint = rep(NA, dim(covariates)[2]), 
                                         EvidenceMargin = rep(NA, dim(covariates)[2])) )}}else{
                                           # IDENTIFYING AND REMOVING TEMPORAL SERIES INSIDE THE CLUSTER WITH ZERO VALUES
                                           zero.series <- unique(which(series==0, arr.ind = TRUE)[,2]) 
                                           if(length(zero.series)==0){series.def <- series}else{series.def <- series[,-(zero.series)]}
                                           #CHECKING THE SIZE OF THE CLUSTER: q=1 or q>1
                                           #is.vector(series.def)==TRUE THEN q=1 OTHERWISE q>1
                                           if(is.vector(series.def)){ 
                                             series.def <- matrix((series.def - mean(series.def))/sd(series.def), ncol=1)
                                             #PRIOR HYPERPARAMETERS FOR q1=1
                                             m01   <- matrix(rep(m0, dim(covariates)[2]*dim(series.def)[2]), ncol=1)
                                             Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                             S01 <- diag(rep(S0,dim(series.def)[2]))
                                             #DISCOUNT FACTORS MATRIX
                                             delta1 <- sqrt(delta)
                                             Beta1 <-  diag(1/c(rep(delta1, dim(covariates)[2])))}else{
                                               series.def <- apply(series.def, 2, function(x){(x-mean(x))/sd(x)})
                                               #PRIOR HYPERPARAMETERS FOR q1>1
                                               m01   <- matrix(rep(m0, dim(covariates)[2]*dim(series.def)[2]), ncol=dim(series.def)[2])
                                               Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                               S01 <- diag(rep(S0,dim(series.def)[2]))
                                               delta1 <- sqrt(delta)
                                               #DISCOUNT FACTORS MATRIX
                                               Beta1 <-  diag(1/c(rep(delta1, dim(covariates)[2])))
                                               
                                             }
                                           
                                           
                                           
                                           if(Test=="LTT"){
                                             res <- .Individual_FunctionalTestLT(ffd1 = as.matrix(series.def), Cova = as.matrix(covariates), m0In = m01, c0In = Cova1, S0In = S01,
                                                                                 beta0In = Beta1, nt0In = n0, NIn = N1, Nsimu = Nsimu1, CUTpos = Cutpos1)
                                             
                                             #EVIDENCE OF ACTIVATION FOR A SINGLE VOXEL TAKING INTO ACCOUNT THE INFORMATION OF THE ENTIRE CLUSTER OF SIZE q
                                             return(res)
                                           }
                                           
                                           if(Test=="JointTest"){
                                             res <- .Individual_FunctionalMultiTest(ffd1 = as.matrix(series.def), Cova = as.matrix(covariates), m0In = m01, c0In = Cova1, 
                                                                                    S0In = S01, beta0In = Beta1, nt0In = n0, NIn = N1, Nsimu = Nsimu1, CUTpos = Cutpos1)
                                             
                                             #EVIDENCE OF ACTIVATION FOR A SINGLE VOXEL TAKING INTO ACCOUNT THE INFORMATION OF THE ENTIRE CLUSTER OF SIZE q
                                             return(res)
                                             
                                             
                                           }
                                           
                                           
                                         }
  }
  
}
#END FUNCTION



