

ffdAverageGroupVoxelFEST <- function(posi.ffd, DatabaseGroup, covariates, m0, Cova, delta, S0, n0, N1, Nsimu1, r1, Cutpos){
  
  
  if(r1 == 0){
    
    posi <- distanceNeighbors(posi.refer = posi.ffd, r1)
    #posi <- rbind(c(23, 32, 45))
    Ngroup <- length(DatabaseGroup)
    #Ngroup <- 3
    #BOLD RESPONSE SERIES FOR A SPECIFIC CLUSTER
    series.group = NULL 
    for(i in 1:Ngroup){
      ffd.c <- DatabaseGroup[[i]]
      #ffd.c <- readNIfTI(paste("~/Documents/TEsis/sub-7.feat/std.nii.gz", sep=""), reorient = FALSE)
      series <- sapply(1:dim(posi)[1], function(ii){ffd.c[posi[ii,1], posi[ii,2], posi[ii,3], ]})
      series.group <- cbind(series.group, series)
    }
    
    #series.group0 = read.csv("seriesgroup2.csv", sep = ",")
    #str(series.group0)
    #series.group = cbind(series.group0[1:310,1], series.group0[311:620,1], series.group0[621:930,1])
    #covariates = as.matrix(read.csv("~/Documents/TEsis/sub-7.feat/covariables.csv", header=FALSE, sep=""))
    
      if(any(series.group==0)){return( list(EvidenceLTT = rep(NA, dim(covariates)[2]), 
                                            EvidenceMargin = rep(NA, dim(covariates)[2])) )}else{
                                             
                                             #m0=0; Cova = 100; S0 =1; delta = 0.95; n0 = 1; N1 =310; Nsimu1 = 100; Cutpos1 = 30
                                             
                                             m01   <- matrix(rep(m0, dim(covariates)[2]), ncol=1)
                                             Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                             S01 <- diag(rep(S0,1))
                                             #DISCOUNT FACTORS MATRIX
                                             delta1 <- sqrt(delta)
                                             Beta1 <-  diag(1/c(rep(delta1, dim(covariates)[2])))
                                             
                                             OnlineThetaM = array(0, c(dim(covariates)[2], N1-Cutpos, Nsimu1))
                                           
                                             
                                            for(j in 1:Ngroup){
                                              x = as.matrix(series.group[,j])
                                              series.def <- matrix((x - mean(x))/sd(x), ncol=1)
                                              res = .Individual_FunctionalMultiTest(ffd1 = series.def, Cova = covariates, m0In = m01, c0In = Cova1, 
                                                                                   S0In = S01, beta0In = Beta1, nt0In = n0, NIn = N1, Nsimu = Nsimu1, CUTpos = Cutpos)
                                              OnlineThetaM = OnlineThetaM + res$Online_theta
                                             
                                            }
                                             
                                             OnlineThetaM = OnlineThetaM/Ngroup
                                             evidence = 1 - sapply(1:dim(covariates)[2],function(j){mean(sapply(1:(Nsimu1), function(i) any(OnlineThetaM[j,,i]<0)))})
                                             
                                             
                                             return( list(EvidenceLTT = rep(NA, dim(covariates)[2]), 
                                                          EvidenceMargin = evidence) )
                                             
                                             }
      
    }else{Ngroup <- length(DatabaseGroup)
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
         
         # IDENTIFYING AND REMOVING TEMPORAL SERIES INSIDE THE CLUSTER WITH ZERO VALUES
         #series.group = read.csv("seriesgroup2.csv", sep = ",")
         #str(series.group0)
         #series.group = cbind(series.group0[1:310,1], series.group0[311:620,1], series.group0[621:930,1])
         #covariates = as.matrix(read.csv("~/Documents/TEsis/sub-7.feat/covariables.csv", header=FALSE, sep=""))
         
         if(any(series.group[,1]==0)){return( list(EvidenceLTT = rep(NA, dim(covariates)[2]), 
                                                   EvidenceMargin = rep(NA, dim(covariates)[2])) )}else{

                                            
                                           zero.series <- unique(which(series.group==0, arr.ind = TRUE)[,2]) 
                                           if(length(zero.series)==0){series.group <- series.group}else{series.group <- as.matrix(series.group[,-(zero.series)])}                  
                                                     
                                            
                                           if(dim(series.group)[2]>1){
                                                    
                                            m01   <- matrix(rep(m0, dim(covariates)[2]*dim(series.group)[2]), ncol=dim(series.group)[2])        
                                            Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                            S01 <- diag(rep(S0, dim(series.group)[2]))
                                            delta1 <- sqrt(delta)
                                            Beta1 <- diag(1/c(rep(delta1, dim(covariates)[2])))
                                            
                                            OnlineThetaM = array(0, c(dim(covariates)[2], N1-Cutpos, Nsimu1))
                                            OnlineThetaL = array(0, c(N1-Cutpos, dim(covariates)[2], Nsimu1))
                                            
                                            for(j in 1:Ngroup){
                                              
                                              x = as.matrix(series.group[(1+(j-1)*dim(covariates)[1]):(j*dim(covariates)[1]),])
                                              series.def <- apply(x, 2, function(y){(y-mean(y))/sd(y)})
                                              
                                              res = .Individual_FunctionalMultiTest(ffd1 = series.def, Cova = covariates, m0In = m01, c0In = Cova1, 
                                                                                   S0In = S01, beta0In = Beta1, nt0In = n0, NIn = N1, Nsimu = Nsimu1, CUTpos = Cutpos)
                                              
                                              res2 = .Individual_FunctionalTestLT(ffd1 = series.def, Cova = covariates, m0In = m01, c0In = Cova1, 
                                                                                   S0In = S01, beta0In = Beta1, nt0In = n0, NIn = N1, Nsimu = Nsimu1, CUTpos = Cutpos)
                                              
                                              OnlineThetaM = OnlineThetaM + res$Online_theta
                                              OnlineThetaL = OnlineThetaL + res2$Online_theta
                                            }
                                            
                                            OnlineThetaM = OnlineThetaM/Ngroup
                                            OnlineThetaL = OnlineThetaL/Ngroup
                                            
                                            evidenceM = 1 - sapply(1:dim(covariates)[2],function(j){mean(sapply(1:(Nsimu1), function(i) any(OnlineThetaM[j,,i]<0)))})
                                            evidenceL = 1 - sapply(1:dim(covariates)[2],function(j){mean(sapply(1:(Nsimu1), function(i) any(OnlineThetaL[,j,i]<0)))})
                                            
                                            return( list(EvidenceLTT = evidenceL, 
                                                         EvidenceMargin = evidenceM) )
                                            
                                           }else{
                                             
                                             m01   <- matrix(rep(m0, dim(covariates)[2]*dim(series.group)[2]), ncol=dim(series.group)[2])        
                                             Cova1 <- diag(rep(Cova, dim(covariates)[2]))
                                             S01 <- diag(rep(S0, dim(series.group)[2]))
                                             delta1 <- sqrt(delta)
                                             Beta1 <- diag(1/c(rep(delta1, dim(covariates)[2])))
                                             
                                             OnlineThetaM = array(0, c(dim(covariates)[2], N1-Cutpos, Nsimu1))
                                             
                                             
                                             for(j in 1:Ngroup){
                                               
                                               x = as.matrix(series.group[(1+(j-1)*dim(covariates)[1]):(j*dim(covariates)[1]),])
                                               series.def <- apply(x, 2, function(y){(y-mean(y))/sd(y)})
                                               
                                               res = .Individual_FunctionalMultiTest(ffd1 = series.def, Cova = covariates, m0In = m01, c0In = Cova1, 
                                                                                     S0In = S01, beta0In = Beta1, nt0In = n0, NIn = N1, Nsimu = Nsimu1, CUTpos = Cutpos)
                                               
                                               
                                               OnlineThetaM = OnlineThetaM + res$Online_theta
                                               
                                             }
                                             
                                             OnlineThetaM = OnlineThetaM/Ngroup
                                             
                                             
                                             evidenceM = 1 - sapply(1:dim(covariates)[2],function(j){mean(sapply(1:(Nsimu1), function(i) any(OnlineThetaM[j,,i]<0)))})
                                             evidenceL = rep(NA, dim(covariates)[2])
                                             
                                             return( list(EvidenceLTT = evidenceL, 
                                                          EvidenceMargin = evidenceM) )
                                             
                                             
                                           }

                                            
  
                                            
                                            
                                          }
                                         
                                         
                                         
                                         
                                         
                                  
                                         
        }
  
  
  
  
}
