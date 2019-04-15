#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "idptSTM.h"
#include "PoissonRegression.h"

template <class T1, class T2>
std::multimap<T2, T1> swapPairs(std::map<T1, T2> m) {
  std::multimap<T2, T1> m1;
  
  for (auto&& item : m) {
    m1.emplace(item.second, item.first);
  }
  
  return m1;
};

void gen_model(arma::uvec& v, double& OldCrt,
               std::map< std::vector<int>, int>& CandidateModels, 
               std::map< std::vector<int>, double>& CriterionRecord,
               const arma::vec& theta, const int& p, const arma::vec& skip, 
               const arma::mat& x, const arma::mat& y, const int maxit){
  
  double crt, lik; 
  auto it = skip.cbegin();
  
  for(int pp = 0; pp != p; ++pp){
    if(!skip.size() || it == skip.cend() || *it != pp){
      v(pp) = 1 - v(pp);
      // get crt for alternative v
      auto v_key_iter = CriterionRecord.find(arma::conv_to<std::vector<int> >::from(v)); 
      if(v_key_iter == CriterionRecord.end()){
        lik = Poisson_Newton(x, y, theta, v, maxit);
        crt = 2*lik + log(y.size()) * sum(v);
        CriterionRecord.insert(std::make_pair(arma::conv_to<std::vector<int> >::from(v), crt));
      }else{
        crt = v_key_iter->second;
      }
      
      // generate model
      double s = exp(crt - OldCrt), 
        prob = Rcpp::traits::is_infinite<REALSXP>(s) ? 1 :( v(pp) ? 1/(1 + s) : s/(1 + s));
      double v_pp = R::rbinom(1, prob);
      if(v_pp == v(pp)){ // keep updated v, update OldCriterion
        OldCrt = crt; 
      }else{  // get back to previous v
        v(pp) = 1 - v(pp);
      }
    }else{ 
      ++it;}
  }
  
  // override CandidateModels
  ++CandidateModels[arma::conv_to<std::vector<int> >::from(v)];  
  
}

// [[Rcpp::export]]
Rcpp::List logGLMselect_cpp(const arma::vec& y, const arma::mat& x, const int maxit, 
                            const arma::vec& skip, const int ModelCnt, bool Message){
  
  // full model
  int p = x.n_cols, tmp_maxit = maxit; bool happy = true;
  arma::vec theta(p, arma::fill::zeros); theta(0) = sum(y)/x.n_rows;
  double lik = Poisson_Newton(x, y, theta, tmp_maxit, happy); //override theta
  double crt = 2*lik + log(y.n_rows) * p;
  if(Message) Rcpp::Rcout << "glm initial converged in " << maxit - tmp_maxit << " iterations" << std::endl;
  
  // model selection
  arma::uvec v(p, arma::fill::ones);
  
  std::map< std::vector<int>, double> CriterionRecord;
  CriterionRecord.insert(std::make_pair(arma::conv_to<std::vector<int> >::from(v), crt));
  std::map< std::vector<int>, int> CandidateModels;
  
  for(int iteration = 0; iteration != ModelCnt; ++iteration){
    gen_model(v, crt, CandidateModels, CriterionRecord, theta, p, skip, x, y, maxit); // override v_k, crt, CandidateModels, CriterionRecord
  }
  
  std::multimap<int, std::vector<int> > Generated_Models = swapPairs(CandidateModels);
  auto iter = Generated_Models.rbegin();
  std::vector<int> selected = iter->second;
  
  if(Message){
    int count = 0;
    while(count < 5 && iter != Generated_Models.rend()){
      int tmp_cnt = iter->first;
      for(auto pos = Generated_Models.equal_range(tmp_cnt); count < 5 && pos.first != pos.second; ++pos.first){
        std::vector<int> tmp_v = pos.first->second;
        Rcpp::Rcout << "Model ";
        for(auto c : tmp_v){
          Rcpp::Rcout << c;
        }
        Rcpp::Rcout << " appeared " << tmp_cnt << " times" << std::endl; 
        ++count; 
      }
      ++iter;
    }
  }
  
  //record
  v = arma::conv_to<arma::uvec>::from(selected);
  lik = Poisson_Newton(x, y, theta, v, maxit); // overide tmp_theta
  arma::mat hessian(sum(v), sum(v)); 
  get_hessian(x, y, theta, v, hessian); //override hessian
  
  return Rcpp::List::create(Rcpp::Named("likelihood") = lik, 
                            Rcpp::Named("parameters") = theta, 
                            Rcpp::Named("v") = v, 
                            Rcpp::Named("se") = arma::diagvec(arma::inv_sympd(hessian)));
}



// [[Rcpp::export]]
Rcpp::List idpSTModelSelection_cpp(const arma::mat& dat, int n, const int maxit,int ModelCnt, bool Message){  
  
  auto in_data = data_indpt(dat, n);
  
  arma::mat x = in_data["covariates"];
  arma::mat y = in_data["response"];
  int K = in_data["K"];
  
  arma::mat theta((K + 1), K), se((K + 1), K), v(K, K); 
  arma::vec skip(1); skip(0) = 0; 
  double lik = 0; 
  
  for(int k = 0; k != K; ++k){
    
    if(Message) Rcpp::Rcout << "Sensitivity of group "<< k + 1<< ": " << std::endl; 
    
    auto tmp = logGLMselect_cpp(y.col(k), x, maxit, skip, ModelCnt, Message);
    
    if(Message) Rcpp::Rcout << std::endl;
   
    //record
    arma::uvec tmp_v = tmp["v"]; arma::vec tmp_se(K + 1, arma::fill::zeros);
    
    theta.col(k) = Rcpp::as<arma::vec>(tmp["parameters"]);
    tmp_se.elem(arma::find(tmp_v)) = Rcpp::as<arma::vec>(tmp["se"]);
    se.col(k) = tmp_se;
    v.col(k) = (Rcpp::as<arma::vec>(tmp["v"])).tail(K);
    lik += Rcpp::as<double>(tmp["likelihood"]);
  }
  
  return Rcpp::List::create(Rcpp::Named("likelihood") = lik, 
                            Rcpp::Named("parameters") = theta,
                            Rcpp::Named("se") = se, 
                            Rcpp::Named("v") = v);
}
