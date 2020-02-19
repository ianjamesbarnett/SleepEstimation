#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double pi = 3.14159265358979323846264338327950288419716939937510582097494;

// [[Rcpp::export]]
arma::vec hermite_C(int points, double z) {
  double p1 = 1/pow(pi,0.4);
  double p2 = 0;
  double p3;
  for(int j=1;j<=points;j++){
    p3 = p2;
    p2 = p1;
    p1 = z * sqrt(2.0/j) * p2 - sqrt((j - 1.0)/j) * p3;
  }
  double pp = sqrt(2 * points) * p2;
  arma::vec retval(2);
  retval[0]=p1;
  retval[1]=pp;
  return retval;
}

// [[Rcpp::export]]
arma::vec test_C(){
  arma::vec x(10);
  return x;
}

// [[Rcpp::export]]
arma::mat gausshermite_C(int points, int iterlim = 50) {
  arma::vec x(points);
  arma::vec w(points);
  double m = (points + 1.0)/2;
  double z;
  double z1;
  arma::vec p(2);
  for(int i=1;i<=m;i++){
    if (i == 1){
      z= sqrt(2.0 * points + 1) - 2 * pow(2.0 * points + 1,-1.0/6);
    }else if (i == 2){
      z=z - sqrt(points*1.0)/z;
    }else if (i == 3 || i == 4){
      z=1.9 * z - 0.9 * x(i - 2-1);
    }else{
      z=2 * z - x(i - 2-1);
    } 
    for(int j=1;j<=iterlim;j++) {
      z1 = z;
      p = hermite_C(points, z);
      z = z1 - p(0)/p(1);
      if(std::abs(z - z1) <= 0.000000000000001){
        break;
      }
    }
    x(points + 1 - i -1) = -z;
    x(i -1) = z;
    w(i-1) = 2.0/pow(p(1),2);
    w(points + 1 - i -1) = w(i-1);
  }
  arma::mat r;
  r.reshape(points,2);
  r.col(0)=x*sqrt(2.0);
  r.col(1)=w/sum(w);
//colnames(r) = c("Points", "Weights");
  return r;
}

// [[Rcpp::export]]
double IndLik_C(double t_init,double wt,double xs,double xw,double lambda_s,double lambda_w,double mu_s,double mu_w){
  double denom;
  double numer;
  if(t_init+wt>mu_s+24){
    if(t_init<xs){
      denom=1-exp(-lambda_w*(xs-t_init))+exp(-lambda_s*(xs-t_init))-exp(-lambda_s*(xw-t_init))+exp(-lambda_w*(xw-t_init))-exp(-lambda_w*(mu_s+24-t_init))+(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
      numer=(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
      return numer/denom;
    }else if(t_init<xw){
      denom=1-exp(-lambda_s*(xw-t_init))+exp(-lambda_w*(xw-t_init))-exp(-lambda_w*(mu_s+24-t_init))+(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
      numer=(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
      return numer/denom;
    }else{
      denom=1-exp(-lambda_w*(mu_s+24-t_init))+(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
      numer=(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
      return numer/denom;
    }
  }
  if(t_init<xs){
    denom=1-exp(-lambda_w*(xs-t_init))+exp(-lambda_s*(xs-t_init))-exp(-lambda_s*(xw-t_init))+exp(-lambda_w*(xw-t_init))-exp(-lambda_w*(mu_s+24-t_init))+(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
    if(t_init+wt<xs){
      return lambda_w*exp(-lambda_w*(wt))/denom;
    }else if(t_init+wt<xw){
      return lambda_s*exp(-lambda_w*(xs-t_init)-lambda_s*(t_init+wt-xs))/denom;
    }else{
      return lambda_w*exp(-lambda_w*(xs-t_init)-lambda_s*(xw-xs)-lambda_w*(t_init+wt-xw))/denom;
    }
  }else if(t_init<xw){
    denom=1-exp(-lambda_s*(xw-t_init))+exp(-lambda_w*(xw-t_init))-exp(-lambda_w*(mu_s+24-t_init))+(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
    if(t_init+wt<xw){
      return lambda_s*exp(-lambda_s*(wt))/denom;
    }else{
      return lambda_w*exp(-lambda_s*(xw-t_init)-lambda_w*(t_init+wt-xw))/denom;   
    }
  }else{
    denom=1-exp(-lambda_w*(mu_s+24-t_init))+(1/(1-exp(-lambda_s*24)))*(exp(-lambda_s*(mu_s+24-t_init))-exp(-lambda_s*(mu_w+24-t_init)))+(1/(1-exp(-lambda_w*24)))*(exp(-lambda_s*(mu_w+24-t_init))-exp(-lambda_w*(mu_w+48-t_init)));
    return lambda_w*exp(-lambda_w*(wt))/denom;
  }
}


// [[Rcpp::export]]
double JointLikSum_C(arma::mat mat,double x_s,double x_w,double lambda_s,double lambda_w,double mu_s,double mu_w,bool loglik=true) {
  double t1;
  int nrows=mat.n_rows;
  if(loglik){
    t1=0;
    for(int i=1;i<=nrows;i++){
      t1=t1+log(IndLik_C(mat(i-1,0),mat(i-1,1)-mat(i-1,0),x_s,x_w,lambda_s,lambda_w,mu_s,mu_w));
    }
  }else{
    t1=1;
    for(int i=1;i<=nrows;i++){
      t1=t1*IndLik_C(mat(i-1,0),mat(i-1,1)-mat(i-1,0),x_s,x_w,lambda_s,lambda_w,mu_s,mu_w);
    }
  }
  return t1;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#hermite_C(10,2.2)
#gausshermite_C(10)
#mat=matrix(c(15.00000, 15.21286, 15.24486, 15.30929, 15.64363, 16.31569, 16.61824, 16.67921, 16.93295, 17.22514, 17.61230, 17.62596, 18.25076, 18.40756,
#              18.48203, 18.99600, 19.06261, 19.08162, 19.68675, 20.02501, 20.09454, 20.13489, 20.49266, 20.56043, 20.58651, 20.59758, 20.64402, 20.78747,
#              20.81538, 21.38978, 21.42816, 21.60945, 21.66603, 22.23404, 22.29023, 22.93119, 22.95101, 22.95973, 23.27916, 23.67772, 23.88266, 24.06001,
#              24.15701, 24.36459, 24.52154, 24.60729, 24.85133, 25.05457, 27.64829, 28.86993, 31.80698, 32.58085, 32.61966, 33.13295, 33.48600, 33.55530,
#              33.56787, 33.80169, 33.95074, 34.41976, 34.83178, 35.56793, 35.70479, 35.73453, 35.81238, 35.88981, 35.89139, 36.06039, 36.62606, 36.66174,
#              36.76703, 37.15986, 37.28550, 37.33921, 37.46255, 37.53565, 37.98056, 38.52197, 38.82134, 38.83057, 38.93535, 38.97066, 15.21286, 15.24486,
#              15.30929, 15.64363, 16.31569, 16.61824, 16.67921, 16.93295, 17.22514, 17.61230, 17.62596, 18.25076, 18.40756, 18.48203, 18.99600, 19.06261,
#              19.08162, 19.68675, 20.02501, 20.09454, 20.13489, 20.49266, 20.56043, 20.58651, 20.59758, 20.64402, 20.78747, 20.81538, 21.38978, 21.42816,
#              21.60945, 21.66603, 22.23404, 22.29023, 22.93119, 22.95101, 22.95973, 23.27916, 23.67772, 23.88266, 24.06001, 24.15701, 24.36459, 24.52154,
#              24.60729, 24.85133, 25.05457, 27.64829, 28.86993, 31.80698, 32.58085, 32.61966, 33.13295, 33.48600, 33.55530, 33.56787, 33.80169, 33.95074,
#              34.41976, 34.83178, 35.56793, 35.70479, 35.73453, 35.81238, 35.88981, 35.89139, 36.06039, 36.62606, 36.66174, 36.76703, 37.15986, 37.28550,
#              37.33921, 37.46255, 37.53565, 37.98056, 38.52197, 38.82134, 38.83057, 38.93535, 38.97066, 39.09946),byrow=F,ncol=2)
#JointLik(mat,mu_s=25.75,mu_w=31.75,sigma_s=0.25,sigma_w=0.25,rho=0,lambda_s=1.18,lambda_w=4.6,24.53,32.96,loglik=FALSE,incl_rho=TRUE)

#start_time <- Sys.time()
#for(i in 1:100){
#  JointLikSum_C(mat,24.53,32.96,1.18,4.6,25.75,31.75,FALSE)
#}
#end_time <- Sys.time()
#end_time-start_time
#start_time <- Sys.time()
#for(i in 1:100){
#  JointLikSum(mat,24.53,32.96,1.18,4.6,25.75,31.75,FALSE)
#}
#end_time <- Sys.time()
#end_time-start_time
*/
