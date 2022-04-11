#include "distribution.h"
#include "family2.h"

// [[Rcpp::export(.esvd_family)]]
SEXP esvd_family(std::string family)
{
    Distribution* distr;

    // Look up family string
    if(family == "poisson")
        distr = get_poisson();
    else if(family == "neg_binom")
        distr = get_neg_binom();
    else if(family == "neg_binom2")
        distr = get_neg_binom2();
    else
        Rcpp::stop("unimplemented family");

    return Rcpp::XPtr<Distribution>(distr, true);
}
