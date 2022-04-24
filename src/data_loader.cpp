#include "data_loader.h"

using Rcpp::NumericMatrix;
using Rcpp::S4;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using SpMatrix = Eigen::SparseMatrix<double>;
using MapMatrixXd = Eigen::Map<MatrixXd>;
using MapMatrixXi = Eigen::Map<MatrixXi>;
using MapSpMatrix = Eigen::Map<SpMatrix>;

// Create a data loader from an R matrix
// Can be dense or sparse
// [[Rcpp::export(.data_loader)]]
SEXP data_loader(SEXP mat)
{
    DataLoader* loader = nullptr;

    // Test if mat is dgCMatrix
    if (::Rf_isS4(mat))
    {
        S4 spmat_rcpp(mat);
        if (spmat_rcpp.is("dgCMatrix"))
        {
            MapSpMatrix spmat = Rcpp::as<MapSpMatrix>(spmat_rcpp);
            loader = new SparseDataLoader(spmat);
        }
    }

    // Test if mat is a dense matrix
    else if (::Rf_isMatrix(mat))
    {
        // Test the storage type
        if (::Rf_isInteger(mat))
        {
            MapMatrixXi denmat = Rcpp::as<MapMatrixXi>(mat);
            loader = new DenseDataLoader<int>(denmat);
        } else if (::Rf_isNumeric(mat)) {
            MapMatrixXd denmat = Rcpp::as<MapMatrixXd>(mat);
            loader = new DenseDataLoader<double>(denmat);
        }
    }

    // Otherwise, given an error
    else
    {
        Rcpp::stop("data_loader(): unsupported matrix type");
    }

    return Rcpp::XPtr<DataLoader>(loader, true);
}

// Print a short description of the data loader
// Used to implement R function print.esvd_data_loader
// [[Rcpp::export]]
void data_loader_description(SEXP loader_)
{
    Rcpp::XPtr<DataLoader> loader(loader_);
    loader->description();
}

// Test the data loader
// [[Rcpp::export]]
void test_data_loader(SEXP mat)
{
    Rcpp::XPtr<DataLoader> loader(data_loader(mat));
    loader->debug_info();
}

/*

 set.seed(123)
 m = 5
 n = 6
 dat = rpois(m * n, lambda = 10)
 dat[sample(m * n, m * n / 2)] = NA
 x = matrix(as.integer(dat), m)
 eSVD2:::test_data_loader(x)

 x = x + 0.1
 eSVD2:::test_data_loader(x)

 library(Matrix)
 dat[sample(m * n, m * n / 2)] = 0
 x = matrix(dat, m)
 xsp = as(x, "sparseMatrix")
 eSVD2:::test_data_loader(xsp)

 */
