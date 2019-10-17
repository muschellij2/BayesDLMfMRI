#' @export

GroupSingleVoxelFEST <- function(posi.ffd, DatabaseGroup, covariates, m0, Cova, delta, S0, n0, N1, Nsimu1, r1, Test, Cutpos){
  
  if(N1==FALSE){N1 = dim(covariates)[1]}
  
  if(r1 == 0){
    
    posi <- distanceNeighbors(posi.refer = posi.ffd, r1)
    Ngroup <- length(DatabaseGroup)
    #BOLD RESPONSE SERIES FOR A SPECIFIC CLUSTER
    series.group = NULL 
    #system.time(
    for(i in 1:Ngroup){
      ffd.c <- DatabaseGroup[[i]]
      series <- sapply(1:dim(posi)[1], function(ii){ffd.c[posi[ii,1], posi[ii,2], posi[ii,3], ]})
      series.group <- rbind(series.group, series)
    }
    
    #case for single voxels
    if(any(series.group==0)){ 
      if(Test=="LTT"){return(NA)}
      if(Test=="Joint"){return( NA )}}else{
                                       
                                       if(Test=="Joint"){
                                         Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                         delta1<- sqrt(delta)
                                         Beta1 <- diag(1/c(rep(delta1, dim(covariates)[2])))
                                         res   <- .Group_FunctionalMultiTest(ffd1 = series.group, Cova = covariates, m0In = m0, c0In = Cova1, S0In = S0, 
                                                                             beta0In = Beta1, nt0In = n0, flag1 = 0, NIn = N1, NS = Ngroup, Nsimu = Nsimu1, CUTpos = Cutpos)
                                         
                                         return(res)
                                         
                                       }
                                       
                                       if(Test=="LTT"){return(NA)}}}else{
                                         
                                         
                                         
                                         Ngroup <- length(DatabaseGroup)
                                         #THIS LINE RETURN THE POSITIONS OF EACH VOXEL INSIDE THE CLUSTER GIVEN THE DISTANCE r1
                                         posi1 <- distanceNeighbors(posi.refer = posi.ffd, r1)
                                         aux.pos <- dim(DatabaseGroup[[1]])[1:3]
                                         row_sub1 <- apply(posi1, 1, function(row, x1){0 < row & row<=x1}, x1=aux.pos)
                                         
                                         posi <- posi1[apply(t(row_sub1), 1, sum)==3, ]
                                         
                                         
                                         
                                         #BOLD RESPONSE SERIES FOR A SPECIFIC CLUSTER
                                         series.group = NULL 
                                         for(i in 1:Ngroup){
                                           ffd.c <- DatabaseGroup[[i]]
                                           series <- sapply(1:dim(posi)[1], function(ii){ffd.c[posi[ii,1], posi[ii,2], posi[ii,3], ]})
                                           series.group <- rbind(series.group, series)
                                         }
                                         
                                         
                                         if(any(series.group[,1]==0)){if(Test=="LTT"){return(NA)}
                                           if(Test=="Joint"){return(NA)}}else{
                                                                            flag <- any(series.group==0)
                                                                            
                                                                            Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                                                            delta1 <- sqrt(delta)
                                                                            Beta1 <- diag(1/c(rep(delta1, dim(covariates)[2])))
                                                                            
                                                                            if(Test=="LTT"){
                                                                              
                                                                              res <- .Gruop_FunctionalTestLT(series.group, covariates, m0, Cova1, S0, Beta1, n0, sum(flag), N1, Ngroup, Nsimu1, Cutpos)  
                                                                              return(res)
                                                                            }
                                                                            
                                                                            if(Test=="Joint"){
                                                                              
                                                                              res   <- .Group_FunctionalMultiTest(ffd1 = series.group, Cova = covariates, m0In = m0, c0In = Cova1, S0In = S0, 
                                                                                                              beta0In = Beta1, nt0In = n0, flag1 = sum(flag), NIn = N1, NS = Ngroup, Nsimu = Nsimu1, CUTpos = Cutpos)
                                                                              
                                                                              return(res)
                                                                              
                                                                            }
                                                                            
                                                                            
                                                                          }  
                                         
                                       }
  
  
  
  
}
