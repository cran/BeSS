#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List bess_lm(Eigen::MatrixXd X, Eigen::VectorXd y, int T, int max_steps, Eigen::VectorXd beta0) {
  int n=X.rows();
  int p=X.cols();
  double max_T=0.0;
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k;
  }
  vector<int>I(p-T);
  vector<int>A(T);
  vector<int>J(p-T);
  vector<int>B(T);
  Eigen::VectorXd beta=beta0;
  Eigen::VectorXd d=(X.transpose()*(y-X*beta)) /double(n);
  Eigen::VectorXd bd=beta+d;
  bd=bd.cwiseAbs();
  for(int k=0;k<=T-1;k++) {
    max_T=bd.maxCoeff(&A[k]);
    bd(A[k])=0;
  }
  sort (A.begin(),A.end());
  set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
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
    sort (B.begin(),B.end());
    set_difference(E.begin(),E.end(), B.begin(),B.end(),J.begin());
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
List get_A(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd beta, double coef0,int T, Eigen::VectorXd B, Eigen::VectorXd weights) {
  double max_T=0.0;
  int n=X.rows();
  int p=X.cols();
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k+1;
  }
  vector<int>A(T);
  vector<int>I(p-T);
  Eigen::VectorXd A_out(T);
  Eigen::VectorXd I_out(p-T);
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
  Eigen::VectorXd pr=xbeta_exp.cwiseQuotient(one_xbeta_exp);
  Eigen::VectorXd l1=-X.adjoint()*((y-pr).cwiseProduct(weights));
  Eigen::MatrixXd X2=X.adjoint();
  X2=X2.array().square();
  Eigen::VectorXd l2=X2*((pr.cwiseProduct(one-pr)).cwiseProduct(weights));
  Eigen::VectorXd d=-l1.cwiseQuotient(l2);
  if(B.size()<p) {
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
  sort (A.begin(),A.end());
  set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  for(int i=0;i<T;i++)
	  A_out(i) = A[i];
  for(int i=0;i<p-T;i++)
	  I_out(i) = I[i];
  return List::create(Named("p")=pr,Named("A")=A_out,Named("I")=I_out,Named("max_T")=max_T);
}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List getcox_A(Eigen::MatrixXd X, Eigen::MatrixXd y, Eigen::VectorXd beta, int T, Eigen::VectorXd B, Eigen::VectorXd status, Eigen::VectorXd weights) {
  double max_T=0.0;
  int n=X.rows();
  int p=X.cols();
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k+1;
  }
  vector<int>A(T);
  vector<int>I(p-T);
  Eigen::VectorXd A_out(T);
  Eigen::VectorXd I_out(p-T);
  Eigen::VectorXd l1(p);
  Eigen::VectorXd l2(p);
  Eigen::VectorXd theta=X*beta;
  theta=weights.array()*theta.array().exp();
  Eigen::VectorXd cum_theta=Eigen::VectorXd::Zero(n);
  cum_theta(n-1)=theta(n-1);
  for(int k=n-2;k>=0;k--) {
    cum_theta(k)=cum_theta(k+1)+theta(k);
  }
  Eigen::MatrixXd xtheta(n,p);
  for(int k=0;k<=p-1;k++) {
    xtheta.col(k)=theta.cwiseProduct(X.col(k));
  }
  Eigen::MatrixXd x2theta(n,p);
  for(int k=0;k<=p-1;k++) {
    x2theta.col(k)=X.col(k).cwiseProduct(xtheta.col(k));
  }
  for(int k=n-2;k>=0;k--) {
    xtheta.row(k)=xtheta.row(k+1)+xtheta.row(k);
  }
  for(int k=n-2;k>=0;k--) {
    x2theta.row(k)=x2theta.row(k+1)+x2theta.row(k);
  }
  for(int k=0;k<=p-1;k++) {
    xtheta.col(k)=xtheta.col(k).cwiseQuotient(cum_theta);
  }
  for(int k=0;k<=p-1;k++) {
    x2theta.col(k)=x2theta.col(k).cwiseQuotient(cum_theta);
  }
  x2theta=x2theta.array()-xtheta.array().square().array();
  xtheta=X.array()-xtheta.array();
  if(status.size()>0) {
    for(int k=0;k<=status.size()-1;k++) {
      xtheta.row(status(k)-1)=Eigen::VectorXd::Zero(p);
      x2theta.row(status(k)-1)=Eigen::VectorXd::Zero(p);
    }
  }
  l1=-xtheta.adjoint()*weights;
  l2=x2theta.adjoint()*weights;
  Eigen::VectorXd d=-l1.cwiseQuotient(l2);
  if(B.size()<p) {
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
  sort (A.begin(),A.end());
  set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  for(int i=0;i<T;i++)
	  A_out(i) = A[i];
  for(int i=0;i<p-T;i++)
	  I_out(i) = I[i];
  return List::create(Named("A")=A_out,Named("I")=I_out,Named("max_T")=max_T);
}
