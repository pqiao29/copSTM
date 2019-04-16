#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <map>

#include "labeled_pairs.h"
#include "make_cor.h"
#include "PoissonRegression.h"
#include "mle.h"
#include "boot_penalty.h"
#include "mle_sub.h"
// [[Rcpp::plugins(cpp11)]]
template <class T1, class T2>
std::multimap<T2, T1> swapPairs(std::map<T1, T2> m) {
  std::multimap<T2, T1> m1;
  
  for (auto&& item : m) {
    m1.emplace(item.second, item.first);
  }
  
  return m1;
};

/*
 * Overloading ModelSelection w.r.t standard error
 */
// without standard error
void ModelSelection(std::vector<int>& v, double& OldCriterion,
                    std::map< std::vector<int>, int>& CandidateModels,
                    std::map< std::vector<int>, double>& CriterionRecord,
                    int B, int maxit, double eps,
                    const int n, int d,const int t_size, const int K,
                    const arma::mat& xx, const arma::vec& y, const arma::vec& y_0,
                    const arma::vec& beta, const std::vector<double>& rho_v,
                    const int& p, const int& p_main,
                    const std::multimap<int, std::vector<int> >& labeled_pairs,
                    const std::multimap<int, std::vector<int> >& labeled_pairs0, 
                    const double add_penalty, const int& iteration,
                    bool Message_prog, const int cor_type){
  
  double lik, crt; 
  bool first_prt = true;
  // override v, OldCriterion, CriterionRecord
  for(int pp = 0; pp != p; ++pp){
    if(pp >= p_main || pp % (K + 1) != 0){
      
      v[pp] = 1 - v[pp];
      
      // get criterion for alternative v
      auto v_key_iter = CriterionRecord.find(v); 
      if(v_key_iter == CriterionRecord.end()){ // v has not been recorded
        
        arma::Col<int> v_main(std::vector<int>(v.cbegin(), v.cbegin() + p_main)), 
        v_rho(std::vector<int>(v.cbegin() + p_main, v.cend()));
        
        int p_main_sub = sum(v_main), p_sub = p_main_sub + sum(v_rho);
        
        mle_sub(lik, xx, y, beta, rho_v, v_main, v_rho, labeled_pairs, maxit, eps); // override lik
        
        // make correlation matrix with zeros in rho
        if(!isnan(lik) && !isinf(lik)){
          
          if(Message_prog && first_prt){
            Rcpp::Rcout << "iteration: " << iteration + 1<< std::endl;
            first_prt = false; 
          }
          
          arma::vec rho_zeros = arma::vec(rho_v) % v_rho;
          std::vector<double> rho_sub = arma::conv_to<std::vector<double> >::from(rho_zeros);  
          arma::mat corr = cor_mat(rho_sub, labeled_pairs0, d);
          //criterion
          if(Message_prog) Rcpp::Rcout << "Evaluating the " << pp + 1 << "th parameter: ";
          double d_star = boot_CLIC_penalty_sub(y_0, n, K, t_size, beta, corr, rho_v, v_main, v_rho, p_main_sub, p_sub, labeled_pairs, B, Message_prog);
          int nn; 
          if(cor_type == 1){
            nn = t_size*K;
          }else{
            if(cor_type == 2){
              nn = t_size*n*n;
            }else{ nn = t_size; }
          }
          crt = - 2*lik + log(nn)*d_star + 2*add_penalty*d_star*log(p_sub);
        }else{
          crt = arma::datum::inf;
        }
        // keep record
        CriterionRecord.insert(std::make_pair(v, crt));
        
      }else{
        crt = v_key_iter->second; // obtain criterion from CriterionRecord
      }
      
      // generate Gibbs sample 
      if(!isnan(crt) && !isinf(crt)){
        double s = exp(crt - OldCriterion), 
          prob = Rcpp::traits::is_infinite<REALSXP>(s) ? 1 :( v[pp] ? 1/(1 + s) : s/(1 + s));
        
        double v_pp = R::rbinom(1, prob);
        
        if(v_pp == v[pp]){ // keep updated v, update OldCriterion
          OldCriterion = crt; 
        }else{  // get back to previous v
          v[pp] = 1 - v[pp];
        }
      }else{
        v[pp] = 1 - v[pp];
      }
    }
  }
  
  // override CandidateModels
  ++CandidateModels[v];              
}


// with standard error
void ModelSelection(std::vector<int>& v, double& OldCriterion,
                    std::map< std::vector<int>, int>& CandidateModels,
                    std::map< std::vector<int>, arma::vec>& ModelStdErr,
                    std::map< std::vector<int>, double>& CriterionRecord,
                    int B, int maxit, double eps,
                    const int n, int d,const int t_size, const int K,
                    const arma::mat& xx, const arma::vec& y, const arma::vec& y_0,
                    const arma::vec& beta, const std::vector<double>& rho_v,
                    const int& p, const int& p_main,
                    const std::multimap<int, std::vector<int> >& labeled_pairs,
                    const std::multimap<int, std::vector<int> >& labeled_pairs0, 
                    const double add_penalty, const int& iteration,
                    bool Message_prog, const int cor_type){
  
  double lik, crt; 
  bool first_prt = true;
  arma::vec tmp_se, se;
  // override v, OldCriterion, CriterionRecord
  for(int pp = 0; pp != p; ++pp){
    if(pp >= p_main || pp % (K + 1) != 0){
      
      v[pp] = 1 - v[pp];
      
      // get criterion for alternative v
      auto v_key_iter = CriterionRecord.find(v); 
      if(v_key_iter == CriterionRecord.end()){ // v has not been recorded
        
        arma::Col<int> v_main(std::vector<int>(v.cbegin(), v.cbegin() + p_main)), 
        v_rho(std::vector<int>(v.cbegin() + p_main, v.cend()));
        
        int p_main_sub = sum(v_main), p_sub = p_main_sub + sum(v_rho);
        
        mle_sub(lik, xx, y, beta, rho_v, v_main, v_rho, labeled_pairs, maxit, eps); // override lik
        
        // make correlation matrix with zeros in rho
        if(!isnan(lik) && !isinf(lik)){
          
          if(Message_prog && first_prt){
            Rcpp::Rcout << "iteration: " << iteration + 1 << std::endl;
            first_prt = false; 
          }
          
          arma::vec rho_zeros = arma::vec(rho_v) % v_rho;
          std::vector<double> rho_sub = arma::conv_to<std::vector<double> >::from(rho_zeros);  
          arma::mat corr = cor_mat(rho_sub, labeled_pairs0, d);
          // bootstrap: criterion and standard error
          if(Message_prog) Rcpp::Rcout << "Evaluating the " << pp + 1 << "th parameter: ";
          double d_star = boot_CLIC_penalty_sub(y_0, n, K, t_size, beta, tmp_se, corr, rho_v, v_main, v_rho, p_main_sub, p_sub, labeled_pairs, B, Message_prog);
          int nn; 
          if(cor_type == 1){
            nn = t_size*K;
          }else{
            if(cor_type == 2){
              nn = t_size*n*n;
            }else{ nn = t_size; }
          }
          crt = - 2*lik + log(nn)*d_star + 2*add_penalty*d_star*log(p_sub);
        }else{
          crt = arma::datum::inf;
        }
        // keep record
        CriterionRecord.insert(std::make_pair(v, crt));
        
      }else{
        crt = v_key_iter->second; // obtain criterion from CriterionRecord
      }
      
      // generate Gibbs sample 
      if(!isnan(crt) && !isinf(crt)){
        double s = exp(crt - OldCriterion), 
          prob = Rcpp::traits::is_infinite<REALSXP>(s) ? 1 :( v[pp] ? 1/(1 + s) : s/(1 + s));
        
        double v_pp = R::rbinom(1, prob);
        
        if(v_pp == v[pp]){ // keep updated v, update OldCriterion and standard error
          OldCriterion = crt; 
          se = tmp_se;
        }else{  // get back to previous v
          v[pp] = 1 - v[pp];
        }
      }else{
        v[pp] = 1 - v[pp];
      }
    }
  }
  
  // override CandidateModels and ModelStdErr
  ++CandidateModels[v];        
  ModelStdErr.insert(std::make_pair(v, se));
}


// [[Rcpp::export]]
Rcpp::List copSTModelSelect_cpp(const arma::mat& x, const arma::vec& y, 
                                const int cor_type, int K, int n, 
                                int ModelCnt, int B, int maxit1, int maxit2, 
                                const double add_penalty = 0, 
                                bool Message_prog = true, bool Message_res = true,
                                bool std_err = false, double eps = 0.1){  
  
  // Input xx need to include the 1 column for intercept
  int d = K*n*n, p_main = x.n_cols; 
  // Initialize
  int tmp_maxit = maxit1;
  auto ini = PoissonRegression(x, y, maxit1); //overwrite maxit1
  if(!ini["happyending"]) throw Rcpp::exception("Unsuccessful glm fit.", false);
  arma::vec beta = ini["paramters"];
  if(Message_prog) Rcpp::Rcout << "glm initial converged in " << tmp_maxit - maxit1 << " iterations" << std::endl;
  
  //labeled_pairs
  int t_size = y.size()/d, p_rho;
  const std::multimap<int, std::vector<int> > labeled_pairs = get_pairs(K, n, p_rho, t_size, cor_type); //override p_rho
  std::vector<double> rho_v(p_rho, 0.0);
  int p = p_rho + p_main;
  
  // Full model 
  double lik = mle(x, y, beta, rho_v, labeled_pairs, maxit1, eps); // override beta0, beta, rho_v
  
  const std::multimap<int, std::vector<int> > labeled_pairs0 = get_pairs(K, n, p_rho, 1, cor_type); 
  arma::mat corr = cor_mat(rho_v, labeled_pairs0, d);
  
  arma::vec y_0 = y.head(d), se;
  
  double d_star = boot_CLIC_penalty(y_0, n, K, t_size, beta, se, corr, rho_v, p_main, p, labeled_pairs, B, Message_prog);
  double criterion = - 2*lik + log(t_size)*d_star + 2*add_penalty*d_star*log(p);
  
  // Model selection -------------------------------------------------
  std::vector<int> v(p, 1);
  
  std::map< std::vector<int>, int> CandidateModels;
  std::map< std::vector<int>, arma::vec> ModelStdErr;
  
  std::map< std::vector<int>, double> CriterionRecord;
  CriterionRecord.insert(std::make_pair(v, criterion));
  
  if(std_err){
    for(int iteration = 0; iteration != ModelCnt; ++iteration){
      ModelSelection( v, criterion, CandidateModels, ModelStdErr, CriterionRecord, B, maxit2, eps, 
                      n, d, t_size, K, x, y, y_0, beta, rho_v, p, p_main, 
                      labeled_pairs, labeled_pairs0, add_penalty, iteration, Message_prog, cor_type);
      // override v, OldCriterion, CriterionRecord, CandidateModels and ModelStdErr
    }
  }else{
    for(int iteration = 0; iteration != ModelCnt; ++iteration){
      ModelSelection( v, criterion, CandidateModels, CriterionRecord, B, maxit2, eps, 
                      n, d, t_size, K, x, y, y_0, beta, rho_v, p, p_main, 
                      labeled_pairs, labeled_pairs0, add_penalty, iteration, Message_prog, cor_type);
      // override v, OldCriterion, CriterionRecord, CandidateModels
    }
  }
  
  std::multimap<int, std::vector<int> > Generated_Models = swapPairs(CandidateModels);
  auto iter = Generated_Models.rbegin();
  std::vector<int> selected = iter->second;
  
  // print message
  if(Message_res){
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
  
  // fit selected model again ofr output
  arma::Col<int> v_main(std::vector<int>(selected.cbegin(), selected.cbegin() + p_main)), 
                 v_rho(std::vector<int>(selected.cbegin() + p_main, selected.cend()));
  auto tmp_theta = mle_sub(lik, x, y, beta, rho_v, v_main, v_rho, labeled_pairs, maxit2, eps);

  if(std_err){
    
    // extract standard error from ModelStdErr
    if(sum(v_main) + sum(v_rho) != p){
      auto iter2 = ModelStdErr.find(selected);
      se = iter2->second;
    }

    return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                              Rcpp::Named("main") = tmp_theta["beta"], 
                              Rcpp::Named("rho") = tmp_theta["rho"], 
                              Rcpp::Named("v") = v, 
                              Rcpp::Named("std_err") = se);
  }
  
  return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                            Rcpp::Named("main") = tmp_theta["beta"], 
                            Rcpp::Named("rho") = tmp_theta["rho"], 
                            Rcpp::Named("v") = v);
}



