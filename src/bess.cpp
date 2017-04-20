#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List bess_lm(Eigen::MatrixXd X, Eigen::VectorXd y, int T, int max_steps, Eigen::VectorXd beta0) {
  int n=X.rows();
  int p=X.cols();
  double max_T=0.0;
  std::vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k;
  }
  std::vector<int>I(p-T);
  std::vector<int>A(T);
  std::vector<int>J(p-T);
  std::vector<int>B(T);
  Eigen::VectorXd beta=beta0;
  Eigen::VectorXd d=(X.transpose()*(y-X*beta)) /double(n);
  Eigen::VectorXd bd=beta+d;
  bd=bd.cwiseAbs();
  for(int k=0;k<=T-1;k++) {
    max_T=bd.maxCoeff(&A[k]);
    bd(A[k])=0;
  }
  std::sort (A.begin(),A.end());
  std::set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  for(int l=1;l<=max_steps;l++) {
   Eigen::MatrixXd X_A(n,T);
    for(int mm=0;mm<=T-1;mm++) {
      X_A.col(mm)=X.col(A[mm]);
    }
    Eigen::MatrixXd X_I(n,p-T);
    for(int mm=0;mm<=p-T-1;mm++) {
      X_I.col(mm)=X.col(I[mm]);
    }
    Eigen::VectorXd beta_A=X_A.colPivHouseholderQr().solve(y);
    for(int mm=0;mm<=T-1;mm++) {
      d(A[mm])=0.0;
      beta(A[mm])=beta_A(mm);
    }
    Eigen::VectorXd d_I=(X_I.transpose()*(y-X_A*beta_A))/double(n);
    for(int mm=0;mm<=p-T-1;mm++) {
      beta(I[mm])=0.0;
      d(I[mm])=d_I(mm);
    }
    bd=beta+d;
    bd=bd.cwiseAbs();
    for(int k=0;k<=T-1;k++) {
      max_T=bd.maxCoeff(&B[k]);
      bd(B[k])=0;
    }
    std::sort (B.begin(),B.end());
    std::set_difference(E.begin(),E.end(), B.begin(),B.end(),J.begin());
    if(A==B) { break;} else {
      A=B;
      I=J;
    }
   }
   double mse=(y-X*beta).array().square().sum();
   mse=mse/(2*n);
   //return beta;
   return List::create(Named("beta")=beta,Named("max_T")=max_T,Named("mse")=mse);
}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List get_A(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd beta, double coef0,int T, Eigen::VectorXd B) {
  double max_T=0.0;
  int n=X.rows();
  int m=X.cols();
  std::vector<int>E(m);
  for(int k=0;k<=m-1;k++) {
    E[k]=k+1;
  }
  std::vector<int>A(T);
  std::vector<int>I(m-T);
  Eigen::VectorXd one(n);
  Eigen::VectorXd coef(n);
  for(int i=0;i<=n-1;i++) {
    coef(i)=coef0;
  }
  one = Eigen::VectorXd::Ones(n);
  Eigen::VectorXd xbeta_exp=(X*beta+coef).array().exp();
  for(int i=0;i<=n-1;i++) {
    if(xbeta_exp(i)>1.0e30)
      xbeta_exp(i)=1.0e30;
  }
  Eigen::VectorXd one_xbeta_exp=one+xbeta_exp;
  Eigen::VectorXd p=xbeta_exp.cwiseQuotient(one_xbeta_exp);
  Eigen::VectorXd l1=-X.adjoint()*(y-p);
  Eigen::MatrixXd X2=X.adjoint();
  X2=X2.array().square();
  Eigen::VectorXd l2=X2*(p.cwiseProduct(one-p));
  Eigen::VectorXd d=-l1.cwiseQuotient(l2);
  if(B.size()<m) {
    for(int k=0;k<=B.size()-1;k++) {
      d(B(k-1))=0;
    }
  }
  Eigen::VectorXd bd=beta+d;
  bd=bd.cwiseAbs();
  bd=bd.cwiseProduct(l2.cwiseSqrt());
  for(int k=0;k<=T-1;k++) {
    max_T=bd.maxCoeff(&A[k]);
    bd(A[k])=0;
    A[k]=A[k]+1;
  }
  std::sort (A.begin(),A.end());
  std::set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  return List::create(Named("p")=p,Named("A")=A,Named("I")=I,Named("max_T")=max_T);
}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List get_A2(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd beta, double coef0, double max_T, Eigen::VectorXd B) {
  double max_t=0.0;
  int n=X.rows();
  int m=X.cols();
  std::vector<int>E(m);
  for(int k=0;k<=m-1;k++) {
    E[k]=k+1;
  }
  std::vector<int>A(m);
  Eigen::VectorXd one(n);
  one = Eigen::VectorXd::Ones(n);
  Eigen::VectorXd coef(n);
  for(int i=0;i<=n-1;i++) {
    coef(i)=coef0;
  }
  Eigen::VectorXd xbeta_exp=(X*beta+coef).array().exp();
  for(int i=0;i<=n-1;i++) {
    if(xbeta_exp(i)>1.0e30)
      xbeta_exp(i)=1.0e30;
  }
  Eigen::VectorXd one_xbeta_exp=one+xbeta_exp;
  Eigen::VectorXd p=xbeta_exp.cwiseQuotient(one_xbeta_exp);
  Eigen::VectorXd l1=-X.adjoint()*(y-p);
  Eigen::MatrixXd X2=X.adjoint();
  X2=X2.array().square();
  Eigen::VectorXd l2=X2*(p.cwiseProduct(one-p));
  Eigen::VectorXd d=-l1.cwiseQuotient(l2);
  if(B.size()<m) {
    for(int k=0;k<=B.size()-1;k++) {
      d(B(k-1))=0;
    }
  }
  Eigen::VectorXd bd=beta+d;
  bd=bd.cwiseAbs();
  bd=bd.cwiseProduct(l2.cwiseSqrt());
  int mark=0;
  for(int k=0;k<=m-1;k++) {
    if(bd(k)>=max_T) {
      A[mark]=k+1;
      mark=mark+1;
    }
  }
  Eigen::VectorXd bd_A(mark);
  for(int k=0;k<=mark-1;k++) {
    bd_A[k]=bd(A[k]-1);
  }
  max_t=bd_A.minCoeff();
  A.resize(mark);
  std::sort(A.begin(),A.end());
  std::vector<int>I(m-mark);
  std::set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  return List::create(Named("p")=p,Named("A")=A,Named("I")=I,Named("max_T")=max_t);
}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List getcox_A(Eigen::MatrixXd X, Eigen::MatrixXd y, Eigen::VectorXd beta, int T, Eigen::VectorXd B, Eigen::VectorXd status) {
  double max_T=0.0;
  int n=X.rows();
  int m=X.cols();
  std::vector<int>E(m);
  for(int k=0;k<=m-1;k++) {
    E[k]=k+1;
  }
  std::vector<int>A(T);
  std::vector<int>I(m-T);
  Eigen::VectorXd l1(m);
  Eigen::VectorXd l2(m);
  Eigen::VectorXd theta=X*beta;
  theta=theta.array().exp();
  Eigen::VectorXd cum_theta=Eigen::VectorXd::Zero(n);
  cum_theta(n-1)=theta(n-1);
  for(int k=n-2;k>=0;k--) {
    cum_theta(k)=cum_theta(k+1)+theta(k);
  }
  Eigen::MatrixXd xtheta(n,m);
  for(int k=0;k<=m-1;k++) {
    xtheta.col(k)=theta.cwiseProduct(X.col(k));
  }
  Eigen::MatrixXd x2theta(n,m);
  for(int k=0;k<=m-1;k++) {
    x2theta.col(k)=X.col(k).cwiseProduct(xtheta.col(k));
  }
  for(int k=n-2;k>=0;k--) {
    xtheta.row(k)=xtheta.row(k+1)+xtheta.row(k);
  }
  for(int k=n-2;k>=0;k--) {
    x2theta.row(k)=x2theta.row(k+1)+x2theta.row(k);
  }
  for(int k=0;k<=m-1;k++) {
    xtheta.col(k)=xtheta.col(k).cwiseQuotient(cum_theta);
  }
  for(int k=0;k<=m-1;k++) {
    x2theta.col(k)=x2theta.col(k).cwiseQuotient(cum_theta);
  }
  x2theta=x2theta.array()-xtheta.array().square().array();
  xtheta=X.array()-xtheta.array();
  if(status.size()>0) {
    for(int k=0;k<=status.size()-1;k++) {
      xtheta.row(status(k)-1)=Eigen::VectorXd::Zero(m);
      x2theta.row(status(k)-1)=Eigen::VectorXd::Zero(m);
    }
  }
  l1=-xtheta.colwise().sum();
  l2=x2theta.colwise().sum();
  Eigen::VectorXd d=-l1.cwiseQuotient(l2);
  if(B.size()<m) {
    for(int k=0;k<=B.size()-1;k++) {
      d(B(k-1))=0;
    }
  }
  Eigen::VectorXd bd=beta+d;
  bd=bd.cwiseAbs();
  bd=bd.cwiseProduct(l2.cwiseSqrt());
  for(int k=0;k<=T-1;k++) {
    max_T=bd.maxCoeff(&A[k]);
    bd(A[k])=0;
    A[k]=A[k]+1;
  }
  std::sort (A.begin(),A.end());
  std::set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  return List::create(Named("A")=A,Named("I")=I,Named("max_T")=max_T);
}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List getcox_A2(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd beta, double max_T, Eigen::VectorXd B, Eigen::VectorXd status) {
  double max_t=0.0;
  int n=X.rows();
  int m=X.cols();
  std::vector<int>E(m);
  for(int k=0;k<=m-1;k++) {
    E[k]=k+1;
  }
  std::vector<int>A(m);
  Eigen::VectorXd one(n);
  one = Eigen::VectorXd::Ones(n);
  Eigen::VectorXd l1(m);
  Eigen::VectorXd l2(m);
  Eigen::VectorXd theta=X*beta;
  theta=theta.array().exp();
  Eigen::VectorXd cum_theta=Eigen::VectorXd::Zero(n);
  cum_theta(n-1)=theta(n-1);
  for(int k=n-2;k>=0;k--) {
    cum_theta(k)=cum_theta(k+1)+theta(k);
  }
  Eigen::MatrixXd xtheta(n,m);
  for(int k=0;k<=m-1;k++) {
    xtheta.col(k)=theta.cwiseProduct(X.col(k));
  }
  Eigen::MatrixXd x2theta(n,m);
  for(int k=0;k<=m-1;k++) {
    x2theta.col(k)=X.col(k).cwiseProduct(xtheta.col(k));
  }
  for(int k=n-2;k>=0;k--) {
    xtheta.row(k)=xtheta.row(k+1)+xtheta.row(k);
  }
  for(int k=n-2;k>=0;k--) {
    x2theta.row(k)=x2theta.row(k+1)+x2theta.row(k);
  }
  for(int k=0;k<=m-1;k++) {
    xtheta.col(k)=xtheta.col(k).cwiseQuotient(cum_theta);
  }
  for(int k=0;k<=m-1;k++) {
    x2theta.col(k)=x2theta.col(k).cwiseQuotient(cum_theta);
  }
  x2theta=x2theta.array()-xtheta.array().square().array();
  xtheta=X.array()-xtheta.array();
  if(status.size()>0) {
    for(int k=0;k<=status.size()-1;k++) {
      xtheta.row(status(k)-1)=Eigen::VectorXd::Zero(m);
      x2theta.row(status(k)-1)=Eigen::VectorXd::Zero(m);
    }
  }
  l1=-xtheta.colwise().sum();
  l2=x2theta.colwise().sum();
  Eigen::VectorXd d=-l1.cwiseQuotient(l2);
  if(B.size()<m) {
    for(int k=0;k<=B.size()-1;k++) {
      d(B(k-1))=0;
    }
  }
  Eigen::VectorXd bd=beta+d;
  bd=bd.cwiseAbs();
  bd=bd.cwiseProduct(l2.cwiseSqrt());
  int mark=0;
  for(int k=0;k<=m-1;k++) {
    if(bd(k)>=max_T) {
      A[mark]=k+1;
      mark=mark+1;
    }
  }
  Eigen::VectorXd bd_A(mark);
  for(int k=0;k<=mark-1;k++) {
    bd_A[k]=bd(A[k]-1);
  }
  max_t=bd_A.minCoeff();
  A.resize(mark);
  std::sort(A.begin(),A.end());
  std::vector<int>I(m-mark);
  std::set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  return List::create(Named("A")=A,Named("I")=I,Named("max_T")=max_t);
}

// [[Rcpp::export]]
List abess_lm(Eigen::MatrixXd X, Eigen::VectorXd y, int T0L, int T0R, int K_max, int max_steps, Eigen::VectorXd beta0, double eps) {
  //initial
  int n=X.rows();
  int p=X.cols();
  List fit_L1,fit_L,fit_R,fit_M,fit_MR;
  List fit_ML;
  int TL;
  TL=T0L;
  int TM;
  int TR;
  TR=T0R;
  double mseL,mseM,mseR,mseML,mseMR;
  Eigen::VectorXd mse_all(K_max+2);
  Eigen::VectorXd lambda(K_max+2);
  Eigen::VectorXd beta0R(p);
  beta0R=beta0;
  Eigen::VectorXd beta0L(p);
  beta0L=beta0;
  Eigen::VectorXd betaL(p);
  Eigen::MatrixXd beta_all=Eigen::MatrixXd::Zero(p,K_max+2);
  Eigen::VectorXd Ths(K_max+2);
  fit_L1=bess_lm(X,y,TL,max_steps,beta0L);
  double nullmse=y.array().square().sum();
  nullmse=nullmse/(2*n);
  fit_L=fit_L1;
  fit_R=bess_lm(X,y,TR,max_steps,beta0R);
  mse_all(0)=fit_L[2];
  mse_all(1)=fit_R[2];
  lambda(0)=fit_L[1];
  lambda(1)=fit_R[1];
  Ths(0)=T0L;
  Ths(1)=T0R;
  Eigen::VectorXd beta0M(p);
  betaL=fit_L[0];
  beta0M=fit_R[0];
  beta_all.col(0)=betaL;
  beta_all.col(1)=beta0M;
  int k;
  for(k=1;k<=K_max;k++) {
    TM=round(TL+0.618*(TR-TL));
    Ths(k+1)=TM;
    fit_M=bess_lm(X,y,TM,max_steps,beta0M);
    mseL=fit_L[2];
    mseM=fit_M[2];
    mseR=fit_R[2];
    mse_all(k+1)=fit_M[2];
    lambda(k+1)=fit_M[1];
    beta0M=fit_M[0];
    beta_all.col(k+1)=beta0M;
    printf("%d-th iteration s.left: %d s.split: %d s.right: %d\n",k,TL,TM,TR);
    if(abs(mseL-mseM)/abs(nullmse*(TM-TL))> eps && abs(mseR-mseM)/abs(nullmse*(TM-TR))< eps) {
      TR=TM;
      fit_R=fit_M;
    } else if(abs(mseL-mseM)/abs(nullmse*(TM-TL))> eps && abs(mseR-mseM)/abs(nullmse*(TM-TR))> eps) {
      TL=TM;
      fit_L=fit_M;
    } else {
      TR=TM;
      fit_R=fit_M;
      TL=T0L;
      fit_L=fit_L1;
    }
    if(TR-TL==1) break;
    fit_ML=bess_lm(X,y,TM-1,max_steps,beta0M);
    fit_MR=bess_lm(X,y,TM+1,max_steps,beta0M);
    mseML=fit_ML[2];
    mseMR=fit_MR[2];
    if(abs(mseML-mseM)/abs(nullmse)> eps && abs(mseMR-mseM)/abs(nullmse)< eps) break;
  }
  return List::create(Named("beta")=beta_all,Named("mse")=mse_all,Named("nullmse")=nullmse, Named("lambda")=lambda,Named("k")=k+2,Named("s.list")=Ths);
}
