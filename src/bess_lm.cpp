#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include "normalize.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;
void bess_lm_pdas(Eigen::MatrixXd& X, Eigen::VectorXd& y, int T0, int max_steps, Eigen::VectorXd& beta, Eigen::VectorXi& A_out, int& l) {
  int n = X.rows();
  int p = X.cols();
  vector<int>A(T0);
  vector<int>B(T0);
  Eigen::MatrixXd X_A(n, T0);
  Eigen::VectorXd beta_A(T0);
  Eigen::VectorXd res = (y-X*beta)/double(n);
  Eigen::VectorXd d(p);
  for(int i=0;i<p;i++){
    d(i) = res.dot(X.col(i));
  }
  Eigen::VectorXd bd = beta+d;
  bd = bd.cwiseAbs();
  for(int k=0;k<T0;k++) {             //update A
    bd.maxCoeff(&A[k]);
    bd(A[k]) = 0.0;
  }
  sort(A.begin(),A.end());
  for(l=1;l<=max_steps;l++) {
    for(int mm=0;mm<T0;mm++) {
      X_A.col(mm) = X.col(A[mm]);
    }
    beta_A = X_A.colPivHouseholderQr().solve(y);  //update beta_A
    beta = Eigen::VectorXd::Zero(p);
    for(int mm=0;mm<T0;mm++) {
      beta(A[mm]) = beta_A(mm);
    }
    res = (y-X_A*beta_A)/double(n);
    for(int mm=0;mm<p;mm++){     //update d_I
      bd(mm) = res.dot(X.col(mm));
    }
    for(int mm=0;mm<T0;mm++) {
      bd(A[mm]) = beta_A(mm);
    }
    bd = bd.cwiseAbs();
    for(int k=0;k<T0;k++) {
      bd.maxCoeff(&B[k]);
      bd(B[k]) = 0.0;
    }
    sort(B.begin(),B.end());
    if(A==B) break;
    else A = B;
  }
  for(int i=0;i<T0;i++){
    A_out(i) = A[i] + 1;
  }
}
// [[Rcpp::export]]
List bess_lm(Eigen::MatrixXd& X, Eigen::VectorXd& y, int T0, int max_steps, Eigen::VectorXd& beta, Eigen::VectorXd& weights, bool normal = true){
  int n = X.rows();
  int p = X.cols();
  int l;
  double mse;
  double nullmse;
  double aic;
  double bic;
  double gic;
  double coef0 = 0.0;
  Eigen::VectorXd meanx(p);
  Eigen::VectorXd normx(p);
  Eigen::VectorXi A_out(T0);
  double meany = 0.0;
  if(normal){
    Normalize(X, y, weights, meanx, meany, normx);
  }
  for(int i=0;i<n;i++){
    X.row(i) = X.row(i)*sqrt(weights(i));
    y(i) = y(i)*sqrt(weights(i));
  }
  bess_lm_pdas(X, y, T0, max_steps, beta, A_out, l);
  mse = (y-X*beta).squaredNorm()/double(n);
  nullmse = y.squaredNorm()/double(n);
  aic = double(n)*log(mse)+2.0*T0;
  bic = double(n)*log(mse)+log(double(n))*T0;
  gic = double(n)*log(mse)+log(double(p))*log(log(double(n)))*T0;
  if(normal){
    beta = sqrt(double(n))*beta.cwiseQuotient(normx);
    coef0 = meany - beta.dot(meanx);
  }
  return List::create(Named("beta")=beta, Named("coef0")=coef0, Named("mse")=mse, Named("nullmse")=nullmse, Named("aic")=aic, Named("bic")=bic, Named("gic")=gic, Named("A")=A_out);
}
// [[Rcpp::export]]
List bess_lms(Eigen::MatrixXd& X, Eigen::VectorXd& y, Eigen::VectorXi& T_list, int max_steps, Eigen::VectorXd& beta0, Eigen::VectorXd& weights, bool warm_start = false, bool normal = true){
  int n = X.rows();
  int p = X.cols();
  int l;
  int m = T_list.size();
  double nullmse;
  Eigen::VectorXd mse(m);
  Eigen::VectorXd aic(m);
  Eigen::VectorXd bic(m);
  Eigen::VectorXd gic(m);
  Eigen::VectorXd coef0(m);
  Eigen::VectorXd meanx(p);
  Eigen::VectorXd normx(p);
  Eigen::VectorXi A_out(p);
  Eigen::VectorXd beta = beta0;
  Eigen::MatrixXd beta_out(p, m);
  double meany = 0.0;
  if(normal){
    Normalize(X, y, weights, meanx, meany, normx);
  }
  int i = 0;
  for(i=0;i<n;i++){
    X.row(i) = X.row(i)*sqrt(weights(i));
    y(i) = y(i)*sqrt(weights(i));
  }
  for(i=0;i<m;i++){
    bess_lm_pdas(X, y, T_list(i), max_steps, beta, A_out, l);
    beta_out.col(i) = beta;
    if(!warm_start) beta = beta0;
    mse(i) = (y-X*beta_out.col(i)).squaredNorm()/double(n);
    aic(i) = double(n)*log(mse(i))+2.0*T_list(i);
    bic(i) = double(n)*log(mse(i))+log(double(n))*T_list(i);
    gic(i) = double(n)*log(mse(i))+log(double(p))*log(log(double(n)))*T_list(i);
  }
  nullmse = y.squaredNorm()/double(n);
  // cout<<m<<endl;
  if(normal){
    for(i=0;i<m;i++){
      beta_out.col(i) = sqrt(double(n))*beta_out.col(i).cwiseQuotient(normx);
      coef0(i) = meany - beta_out.col(i).dot(meanx);
    }
  }
  return List::create(Named("beta")=beta_out, Named("coef0")=coef0, Named("mse")=mse, Named("nullmse")=nullmse, Named("aic")=aic, Named("bic")=bic, Named("gic")=gic);
}
// [[Rcpp::export]]
List bess_lm_gs(Eigen::MatrixXd& X, Eigen::VectorXd& y, int s_min, int s_max, int K_max, int max_steps, double epsilon, Eigen::VectorXd& beta0, Eigen::VectorXd& weights, bool warm_start = false, bool normal = true){
  int n = X.rows();
  int p = X.cols();
  int l;
  int k;
  int sL = s_min;
  int sR = s_max;
  int sM;
  double nullmse;
  double mseL1;
  double mseL;
  double mseR;
  double mseM;
  double mseML;
  double mseMR;
  Eigen::VectorXi T_list(K_max+2);
  Eigen::VectorXd mse(K_max+2);
  Eigen::VectorXd aic(K_max+2);
  Eigen::VectorXd bic(K_max+2);
  Eigen::VectorXd gic(K_max+2);
  Eigen::VectorXd coef0(K_max+2);
  Eigen::VectorXd meanx(p);
  Eigen::VectorXd normx(p);
  Eigen::VectorXi A_out(p);
  Eigen::VectorXd beta_0R = beta0;
  Eigen::VectorXd beta_0L = beta0;
  Eigen::VectorXd beta_0M(p);
  Eigen::VectorXd beta_0MR(p);
  Eigen::VectorXd beta_0ML(p);
  Eigen::MatrixXd beta_out(p, K_max+2);
  double meany = 0.0;
  if(normal){
    Normalize(X, y, weights, meanx, meany, normx);
  }
  for(int i=0;i<n;i++){
    X.row(i) = X.row(i)*sqrt(weights(i));
    y(i) = y(i)*sqrt(weights(i));
  }
  bess_lm_pdas(X, y, sL, max_steps, beta_0L, A_out, l);
  beta_out.col(0) = beta_0L;
  if(warm_start) beta_0R = beta_0L;
  bess_lm_pdas(X, y, sR, max_steps, beta_0R, A_out, l);
  beta_out.col(1) = beta_0R;
  //cout<<beta_out.leftCols(2)<<endl;
  if(warm_start) beta_0M = beta_0R; else beta_0M = beta0;
  beta_0M = beta_0R;
  mse(0) = (y-X*beta_out.col(0)).squaredNorm()/double(n);
  mse(1) = (y-X*beta_out.col(1)).squaredNorm()/double(n);
  mseL = mse(0);
  mseL1 = mseL;
  mseR = mse(1);
  nullmse = y.squaredNorm()/double(n);
  T_list(0) = sL;
  T_list(1) = sR;
  for(k=2;k<=K_max+1;k++){
    sM = round(sL + 0.618*(sR - sL));
    T_list(k) = sM;
    bess_lm_pdas(X, y, sM, max_steps, beta_0M, A_out, l);
    //cout<<sL<<" "<<sM<<" "<<sR<<endl;
    beta_out.col(k) = beta_0M;
    if(!warm_start) beta_0M = beta0;
    mse(k) = (y-X*beta_out.col(k)).squaredNorm()/double(n);
    mseM = mse(k);
    if((abs(mseL - mseM)/abs(nullmse*(sM - sL)) > epsilon) && (abs(mseR - mseM)/abs(nullmse*(sM - sR)) < epsilon)) {
      sR = sM;
      mseR = mseM;
    } else if((abs(mseL - mseM)/abs(nullmse*(sM - sL)) > epsilon) && (abs(mseR - mseM)/abs(nullmse*(sM - sR)) > epsilon)) {
      sL = sM;
      mseL = mseM;
    } else {
      sR = sM;
      mseR = mseM;
      sL = s_min;
      mseL = mseL1;
    }
    if(sR - sL == 1) break;
    if(warm_start) {
      beta_0ML = beta_0M;
      beta_0MR = beta_0M;
    } else {
      beta_0ML = beta0;
      beta_0MR = beta0;
    }
    bess_lm_pdas(X, y, sM-1, max_steps, beta_0ML, A_out, l);
    bess_lm_pdas(X, y, sM+1, max_steps, beta_0MR, A_out, l);
    mseML = (y-X*beta_0ML).squaredNorm()/double(n);
    mseMR = (y-X*beta_0MR).squaredNorm()/double(n);
    if((abs(mseML - mseM)/nullmse > epsilon) && (2*abs(mseMR - mseM)/nullmse < epsilon)) break;
  }
  if(k>K_max+1) k = K_max+1;
  if(normal){
    for(int kk=0;kk<=k;kk++){
      beta_out.col(kk) = sqrt(double(n))*beta_out.col(kk).cwiseQuotient(normx);
      coef0(kk) = meany - beta_out.col(kk).dot(meanx);
    }
  }
  for(int i=0;i<=k;i++){
    aic(i) = double(n)*log(mse(i))+2.0*T_list(i);
    bic(i) = double(n)*log(mse(i))+log(double(n))*T_list(i);
    gic(i) = double(n)*log(mse(i))+log(double(p))*log(log(double(n)))*T_list(i);
  }
  return List::create(Named("beta")=beta_out.leftCols(k+1), Named("coef0")=coef0.head(k+1), Named("s_list")=T_list.head(k+1), Named("mse")=mse.head(k+1), Named("nullmse")=nullmse, Named("aic")=aic.head(k+1), Named("bic")=bic.head(k+1), Named("gic")=gic.head(k+1));
}
