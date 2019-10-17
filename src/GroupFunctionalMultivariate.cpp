// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
# include <RcppArmadillo.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]

Rcpp::List Group_FunctionalMultiTest(const arma::mat & ffd1, const arma::mat & Cova, const double m0In,
                                          const arma::mat & c0In, const double S0In, const arma::mat & beta0In,
                                          const double nt0In, const int flag1, const int NIn, const int NS, const int Nsimu, const int CUTpos){

  arma::mat Y1 = ffd1;
  arma::mat X1 = Cova;
  double m01 = m0In;
  arma::mat c00 = c0In;
  double S01 = S0In;
  arma::mat beta0 = beta0In;
  double nt00 = nt0In;
  int flag = flag1;
  int  N   = NIn;
  int  Ns  = NS;
  //int  j  = Rcpp::as<int>(Serie);


  //NUMEBER OF CURVES DRAWING FROM OUR ALGORITM. -FOR THE FUNCTIONAL TEST-
  int   nsimu = Nsimu;
  //NUMBER OF VOXELS THAT CONFORM THE CLUSTER
  int  q1 = Y1.n_cols;
  //NUMBER OF COVARIATES IN THE MODEL
  int  p1 = X1.n_cols;
  //CUT POSITION FOR THE FUNCTIONAL TEST
  int cutpos = CUTpos;

  
  //DECLARING VECORS TO ALLOCATE THE THETA CURVES FOR THE FUNCTIONAL TEST
  arma::mat  Theta1(p1, N-cutpos);
  arma::cube ThetaOut(p1, N-cutpos, nsimu);
  
  //DECLARING LOGICAL VECTOR TO COMPUTE THE FUNCTIONAL EVIDENCE
  arma::mat res1(p1, nsimu);
  arma::mat res2(p1, nsimu);
  arma::mat Aux_bin(p1, N-cutpos);
  arma::mat Aux_bin2(p1, N-cutpos);
  
  //VECTORS TO ALLOCATE THE APPROXIMATION ERROR (MEASURE OF FITNESS)
  arma::rowvec Efit(N-cutpos);
  arma::rowvec Efit2(nsimu);


  //TAKING OFF THE COLUMNS WITH NULL VALUES (flag==1 MEANS THAT THERE IS ANY COLUMN WITH ZERO VALUES)
  if(flag==1){
    for(int i=0; i<q1; i++){
      if(any(Y1.col(i)==0)){Y1.shed_col(i);
        i=i-1;
        q1=q1-1;}
    }
  }

  //DECLARING ARRAYS TO ALLOCATE THE SIMULATED BOLD RESPONSE
  arma::cube ysimu(nsimu, q1, N-cutpos);
  arma::mat Ysimu(nsimu, q1);

  //CREATING THE MATRIX S0 WITH DIMENSIONS ACCORDING TO THE VALUE OF q1
  arma::mat S00(q1, q1);
  S00.eye();
  S00 = S01*S00;
  //CREATING THE MATRIX m0 WITH DIMENSIONS ACCORDING TO THE VALUE OF q1
  arma::mat m00(p1, q1);
  m00.ones();
  m00 = m01*m00;


  //COPY THE PRIOR HYPERPARAMETERS FOR THE APPLICATION OF OUR ALGORIHTM
  arma::mat m01c = m00;
  arma::mat c01c = c00;
  arma::mat S01c = S00;
  double nt01c   = nt00;

  //DEFINING AUXILIAR VECTORS AND MATRICES
  arma::uvec colAux(q1);
  arma::uvec rowAux(N);
  arma::mat Y_final(N, q1);


  //DEFINING AUXILIAR MATRICES TO ALLOCATE THE MATRICES RESULTING FROM THE GROUP ANALYSIS
  arma::colvec nt0aux(Ns);
  arma::mat m0aux(p1*Ns, q1);
  arma::mat c0aux(p1*Ns, p1);
  arma::mat S0aux(q1*Ns, q1);
  arma::mat c0aux1(p1, p1);
  arma::mat S02(p1*q1, p1*q1);
  arma::mat m011c(p1, q1);
  arma::mat c011c(p1, p1);
  arma::mat S011c(q1, q1);

  //FILL THE MATRICES AND VECTOR HYPERPARAMETERS WITH THE PRIOR VALUES
  for(int j=0; j<Ns; j++){
    nt0aux[j] = nt00;
    m0aux.submat(j*p1, 0, p1*(j+1)-1, q1-1) = m00;
    c0aux.submat(j*p1, 0, p1*(j+1)-1, p1-1) = c00;
    S0aux.submat(j*q1, 0, q1*(j+1)-1, q1-1) = S00;
  }



  //FILLING ONE OF THE AUXILIARY VECTORS WITH THE POSITIONS RELATED TO THE CLUSTER SIZE (q1)
  for(int k=0; k<q1; k++){colAux[k] = k;}


  for(int i = 0; i<N; i++){


    //FITTING THE MULTIVARIATE DLM FOR EVERY SUBJECT IN THE SAMPLE (Ns SUBJECTS)
    for(int j=0; j<Ns; j++){

      double    nt0 = nt0aux[j];
      arma::mat m0  = m0aux.submat(j*p1, 0, p1*(j+1)-1, q1-1);
      arma::mat c0  = c0aux.submat(j*p1, 0, p1*(j+1)-1, p1-1);
      arma::mat S0  = S0aux.submat(j*q1, 0, q1*(j+1)-1, q1-1);

      //FILLING ONE OF THE AUXILIARY VECTORS WITH THE POSITIONS RELATED TO THE BOLD RESPONSES OBSERVATIONS FOR THE SUBJECT j (N OBSERVATIONS).
      for(int k2=0 ; k2< N; k2++){rowAux[k2] = j*N + k2;}
      //TAKING THE APPROPIATE SUBMATRIX (OR NEIGHBORHOOD VOXEL SERIES FOR THE SUBJECT j) USING THE AUXILIARY VECTORS
      Y_final   = Y1.submat(rowAux, colAux);
      //STANDARDIZING THE DOLB RESPONSES
      for(int ii=0; ii<q1; ii++){Y_final.col(ii)   = (Y_final.col(ii) - arma::mean(Y_final.col(ii)))/arma::stddev(Y_final.col(ii));}


      arma::rowvec F1 = X1.row(i);
      arma::mat    a1  = m0;
      arma::mat    R1  = beta0*c0*beta0;
      arma::rowvec f1  = F1*m0;
      double       Q1  = arma::as_scalar(F1*c0*F1.t()) + 1.0;
      arma::rowvec e1  = Y_final.row(i) - f1;
      arma::colvec A1  = R1*F1.t()/Q1;
      int          n1  = nt0 + 1;
      arma::mat    S1  = (nt0*S0 + e1.t()*e1/Q1)/n1;
      arma::mat    C1  = R1 - A1*A1.t()*Q1;
      arma::mat    m1  = a1 + A1*e1;

      nt0aux[j] = n1;
      m0aux.submat(j*p1, 0, p1*(j+1)-1, q1-1)  = m1;
      c0aux.submat(j*p1, 0, p1*(j+1)-1, p1-1)  = C1;
      S0aux.submat(j*q1, 0, q1*(j+1)-1, q1-1)  = S1;


    }

    arma::rowvec F1 = X1.row(i);


    //THE POSTERIOR PARAMETERS AT TIME cutpos ARE STORED AND WILL BE USED AS PRIOR HYPERPARAMETERS ON THE SIMULATION
    //OF THE ON-LINE THETA TRAJECTORIES
    if(i==cutpos-1){S011c.zeros();
      c011c.zeros();
      m011c.zeros();


      for(int jj=0; jj<Ns; jj++){
        m011c = m011c + m0aux.submat(jj*p1, 0, p1*(jj+1)-1, q1-1);
        c011c = c011c + c0aux.submat(jj*p1, 0, p1*(jj+1)-1, p1-1);
        S011c = S011c + S0aux.submat(jj*q1, 0, q1*(jj+1)-1, q1-1);
      }

      m011c = m011c/Ns;
      c011c = c011c/pow(Ns, 2);
      S011c = S011c/pow(Ns, 2);


      m01c     = m011c;
      c01c     = c011c;
      S01c     = S011c;
      nt01c    = nt0aux[1];
    }


    //IN THIS BLOCK ARE SIMULATED nsimu BOLD RESPONSES
    if(i>=(cutpos)){

      S011c.zeros();
      c011c.zeros();
      m011c.zeros();


      for(int jj=0; jj<Ns; jj++){
        m011c = m011c + m0aux.submat(jj*p1, 0, p1*(jj+1)-1, q1-1);
        c011c = c011c + c0aux.submat(jj*p1, 0, p1*(jj+1)-1, p1-1);
        S011c = S011c + S0aux.submat(jj*q1, 0, q1*(jj+1)-1, q1-1);
      }

      m011c = m011c/Ns;
      c011c = c011c/pow(Ns, 2);
      S011c = S011c/pow(Ns, 2);


      for(int kkk=0; kkk<nsimu; kkk++){
        arma::mat Z1 = arma::randn(p1, q1);
        arma::mat thetaSimu  = m011c  + arma::chol(c011c)*Z1*arma::chol(S011c);
        arma::mat Z2 = arma::randn(1.0, q1);
        Ysimu.row(kkk) = F1*thetaSimu +  Z2*arma::chol(S011c);
      }

      ysimu.slice(i-cutpos) = Ysimu;






    }





  }
  //END BLOCK




  //BEGIN(SIMULATE THE ON-LINE THETA TRAJECTORIES)
  for(int j=0; j < nsimu; j++){
    arma::mat m01  = m01c;
    arma::mat c01  = c01c;
    arma::mat S01  = S01c;
    double    nt01 = nt01c;

    for(int i = cutpos; i<N; i++){
      arma::rowvec F1 = X1.row(i);
      arma::mat    a11  = m01;
      arma::mat    R11  = beta0*c01*beta0;
      arma::rowvec f11  = F1*m01;
      double       Q11  = arma::as_scalar(F1*c01*F1.t()) + 1.0;
      arma::rowvec Ysim(q1);
      for(int k=0; k<q1; k++){Ysim(k) = ysimu(j, k, i-cutpos);}
      arma::rowvec e11  = Ysim - f11;
      arma::colvec A11  = R11*F1.t()/Q11;
      int          n11  = nt01 + 1;
      arma::mat    S11  = (nt01*S01 + e11.t()*e11/Q11)/n11;
      arma::mat    C11  = R11 - A11*A11.t()*Q11;
      arma::mat    m11  = a11 + A11*e11;

      nt01 = n11;
      m01  = m11;
      c01  = C11;
      S01  = S11;

      Efit(i-cutpos) = norm(Y1.row(i) - f11)/norm(Y1.row(i));
      
      Theta1.col(i-cutpos) = m01.col(0);
      
      for(int kk=0; kk<p1; kk++){
        if(m01(kk, 0)<0){Aux_bin2(kk, i-cutpos) = 0;}else{Aux_bin2(kk, i-cutpos) = 1;}
        if(any(m01.row(kk)<0)){Aux_bin(kk, i-cutpos) = 0;}else{Aux_bin(kk, i-cutpos) = 1;}
      }
      
    }
    
    Efit2(j) = median(Efit);
    
    
    for(int jj=0; jj<p1; jj++){if(sum(Aux_bin.row(jj))==(N-cutpos)){res1(jj, j) = 1;}else{res1(jj, j) = 0;}}
    for(int jj=0; jj<p1; jj++){if(sum(Aux_bin2.row(jj))==(N-cutpos)){res2(jj, j) = 1;}else{res2(jj, j) = 0;}}
    ThetaOut.slice(j) = Theta1;
    
  }
  //END(SIMULATE THE ON-LINE THETA TRAJECTORIES)
  
  //COMPUTATION OF THE FITNESS MEASURE
  double FitnessV = median(Efit2);
  
  
  //COMPUTE THE FUNCTIONAL EVIDENCE FOR A POSSIBLE BRAIN ACTIVATION
  arma::rowvec evidenTheta1(p1);
  arma::rowvec evidenTheta2(p1);
  for(int i=0; i<p1; i++){evidenTheta1(i) = mean(res1.row(i));
    evidenTheta2(i) = mean(res2.row(i));}
  
  return Rcpp::List::create(Rcpp::Named("EvidenMultivariate") = evidenTheta1,
                            Rcpp::Named("EvidenMarginal") = evidenTheta2,
                            Rcpp::Named("Online_theta") = ThetaOut,
                            Rcpp::Named("Y_simu") = ysimu,
                            Rcpp::Named("FitnessV") =FitnessV);
}
