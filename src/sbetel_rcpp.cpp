// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

NumericMatrix smooth_rcpp(NumericMatrix gMat, int bw, int td) {
  
  double scale = (1.0/(2.0*bw+1.0));
  double scale_inv = pow(scale, -0.5);
  scale = pow(scale, 0.5);
  
  //Possible dummy prior observations are not to be smoothed but 
  //only inflated by beta^(-1)
  int tt = gMat.nrow();
  int nc = gMat.ncol();
  NumericMatrix gMat_data(tt-td,nc);
  NumericMatrix gMat_new(tt, nc);
  if(td > 0) {
    gMat_data = gMat(Range(td,tt-1),_);
    for(int i = 0; i < td; i++) {
      gMat_new(i,_) = gMat(i,_) * scale_inv;
    }
  } else {
    gMat_data = gMat;
  }
  int t = gMat_data.nrow();
  
  for(int i = 1; i < t+1; i++) {
    
    int maxx;
    int minn;
    if(i-t > -bw) {
      maxx = i-t;
    } else {
      maxx = -bw;
    }
    if(i-1 < bw) {
      minn = i-1;
    } else {
      minn = bw;
    }
    int to = (i-maxx)-1;
    int from = (i-minn)-1;
    
    for(int j = 0; j < nc; j++) {
      
      NumericVector subvec(to-from+1);
      int count = 0;
      for(int h = from; h < to+1; h++) {
        subvec(count) = gMat_data(h,j);
        count ++;
      }
      
      gMat_new(i-1+td,j) = sum(subvec) * scale;
      
    }
    
  }
  
  return gMat_new;
}

NumericVector laGrangian_rcpp(NumericMatrix gMat, int itermax, double tol = 1e-10) {
  
  int nc = gMat.ncol();
  int t = gMat.nrow();
  
  NumericVector eta(nc);
  arma::vec eta_arma = as<arma::vec>(eta);
  arma::vec eta_last = eta_arma;
  
  NumericMatrix temp(nc, nc);
  
  for(int i = 0; i < itermax; i++) {
    
    //First derivative
    NumericMatrix temp1(t,nc);
    for(int j = 0; j < t; j++) {
      NumericVector gMat_row = gMat(j,_);
      arma::vec gMat_row_arma = as<arma::vec>(gMat_row);
      temp1(j,_) = as<NumericVector>(wrap(gMat_row_arma * arma::exp(eta_arma.t() * gMat_row_arma)));
    }
    NumericVector d1(nc);
    for(int j = 0; j < nc; j++) {
      NumericVector temp1_col = temp1(_,j);
      arma::vec temp1_col_arma = as<arma::vec>(temp1_col);
      d1[j] = as<double>(wrap(arma::mean(temp1_col_arma)));
    }
    
    //Second derivative (Jacobian)
    NumericMatrix temp2(nc*nc,t);
    for(int j = 0; j < t; j++) {
      NumericVector gMat_row = gMat(j,_);
      arma::vec gMat_row_arma = as<arma::vec>(gMat_row);
      temp2(_,j) = as<NumericVector>(wrap(arma::vectorise(gMat_row_arma * arma::exp(eta_arma.t() * gMat_row_arma) * gMat_row_arma.t())));
    }
    NumericVector temp22(nc*nc);
    for(int j = 0; j < nc*nc; j++) {
      NumericVector temp2_row = temp2(j,_);
      arma::vec temp2_row_arma = as<arma::vec>(temp2_row);
      temp22[j] = as<double>(wrap(arma::mean(temp2_row_arma)));
    }
    NumericMatrix J(nc, nc);
    for(int j = 0; j < nc; j++) {
      int from = j*nc;
      int to = j*nc+nc-1;
      J(_,j) = temp22[Range(from, to)];
    }
    
    //Newton step
    arma::mat J_arma = as<arma::mat>(J);
    arma::vec d1_arma = as<arma::vec>(d1);
    arma::vec direction = arma::solve(J_arma, d1_arma, arma::solve_opts::fast);
    eta_arma = eta_arma - direction;
    
    double diff_norm_arma = arma::norm(eta_arma - eta_last);
    if(diff_norm_arma < tol) {
      break;
    }
    eta_last = eta_arma;
    
    //If no convergence is achieved, empty vector is returned
    if(i == (itermax - 1)) {
      NumericVector empty;
      return empty;
    }
    
    //For choosing efficient tol
    //Rcout << diff_norm_arma << "\n";
    //Rcout << i << "\n";
  }
  
  eta = as<NumericVector>(wrap(eta_arma));
  return eta;
}

// [[Rcpp::export]]
double etel_rcpp(NumericVector th, 
                 Function g, 
                 NumericMatrix y,
                 int bw,
                 int td,
                 int itermax,
                 List args) {
  
  NumericMatrix gMat = g(th, y, args);
  NumericMatrix gMat_smooth = smooth_rcpp(gMat, bw, td);
  NumericVector lambdaHat = laGrangian_rcpp(gMat_smooth, itermax);
  
  if(lambdaHat.size() == 0) {
    return R_NegInf;
  }
  
  arma::vec lambdaHat_arma = as<arma::vec>(lambdaHat);
  arma::mat gMat_smooth_arma = as<arma::mat>(gMat_smooth);
  arma::vec p1 = gMat_smooth_arma * lambdaHat_arma;
  
  //LogSumExp
  double maxx = max(as<NumericVector>(wrap(p1)));
  arma::vec expp = arma::exp(p1 - maxx);
  double summ = arma::sum(expp);
  double p2 = log(summ) + maxx;
  
  double ret = as<double>(wrap(arma::sum(p1 - p2)));
  return ret;
}

