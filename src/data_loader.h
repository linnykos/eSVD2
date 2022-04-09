#ifndef ESVD2_DATA_LOADER_H
#define ESVD2_DATA_LOADER_H

#include <type_traits>
#include <RcppEigen.h>

// Type of the element
enum class Flag { na, regular, zero };

// Abstract class to represent the iterator for a vector,
// typically a row or a column of a matrix
//
// The purpose of this class is to unify the interface of iterations
// over vectors, as the vector can be either dense or sparse.
// To loop over all elements in the vector, including sparse zeros,
// we can do:
//
// for(; iter; ++iter)
// {
//     double val;
//     Flag flag = iter.value(val);
// }
//
// The `flag` variable indicates whether the current element is
// NA, a regular value, or a sparse zero
//
class VecIterator
{
protected:
    std::size_t m_len;
    std::size_t m_pos;

public:
    VecIterator(std::size_t len) : m_len(len), m_pos(0) {}

    // Returns the index of the current element
    std::size_t index() const { return m_pos; }
    // Returns the flag and value for the current element
    // 0 - NA, 1 - val, 2 - sparse zero
    virtual Flag value(double& val) const = 0;
    // Move to the next element
    virtual VecIterator& operator++() { m_pos++; return *this; }
    // Whether to continue the iteration
    operator bool() const { return m_pos < m_len; }
};

// VecIterator for a dense vector of element type T
template <typename T>
class DenseVecIterator: public VecIterator
{
private:
    const T* m_data;

public:
    DenseVecIterator(const T* data, std::size_t len) :
        VecIterator(len), m_data(data)
    {}

    inline Flag value(double& val) const override
    {
        // Get the element of the original type
        const T Tval = m_data[this->m_pos];
        // Cast to double type
        val = static_cast<double>(Tval);
        // Test NA according to the original type
        const bool isna = (std::is_same<T, int>::value) ?
            Rcpp::IntegerVector::is_na(Tval) :
            Rcpp::NumericVector::is_na(val);
        return isna ? Flag::na : Flag::regular;
    }
};

// VecIterator for a sparse vector
// The element type is always double
class SparseVecIterator: public VecIterator
{
private:
    // Suppose the sparse vector is [0, 0, 1, NA, 0, 2, 3]
    // Then m_len = 7, m_pos is initialized to be 0
    const double* m_valptr;    // [1, NA, 2, 3]
    const int* m_indptr;       // [2, 3, 5, 6]
    std::size_t m_nnz;         // 4
    std::size_t m_innerpos;    // position in m_valptr and m_indptr

public:
    SparseVecIterator(const double* valptr, const int* indptr, std::size_t nnz, std::size_t len) :
        VecIterator(len), m_valptr(valptr), m_indptr(indptr), m_nnz(nnz), m_innerpos(0)
    {}

    inline Flag value(double& val) const override
    {
        if (m_nnz < 1)
            return Flag::zero;

        // Index of the next or the last nonzero value
        const std::size_t next_nonzero_ind = static_cast<std::size_t>(m_indptr[m_innerpos]);
        // If m_pos != next_nonzero_ind, then the current element is zero
        if (this->m_pos != next_nonzero_ind)
            return Flag::zero;

        // Otherwise, m_pos == next_nonzero_ind, and we retrieve the value
        val = m_valptr[m_innerpos];
        return Rcpp::NumericVector::is_na(val) ? Flag::na : Flag::regular;
    }

    inline VecIterator& operator++()
    {
        // Note that this condition also implies m_nnz > 0
        if (this->m_innerpos < m_nnz - 1)
        {
            // Index of the next nonzero value
            const std::size_t next_nonzero_ind = static_cast<std::size_t>(m_indptr[m_innerpos]);
            // If m_pos == next_nonzero_ind, increase m_innerpos
            if (this->m_pos == next_nonzero_ind)
                m_innerpos++;
        }

        m_pos++;
        return *this;
    }
};



// Abstract class for data loader
//
// The purpose of this class is to unify the interface of data loading,
// as the data matrix can be either dense or sparse. Also, for dense
// matrices, the element type can be int or double.
//
class DataLoader
{
protected:
    // Print the data using iterators, mainly used for debugging
    void print_data_iter(std::size_t nrow, std::size_t ncol)
    {
        for (std::size_t i = 0; i < nrow; i++)
        {
            VecIterator& riter = row_iter(i);
            Rcpp::Rcout << "Row iter " << i << ": ";
            for (; riter; ++riter)
            {
                double val = 0.0;
                int flag = static_cast<int>(riter.value(val));
                Rcpp::Rcout << "(" << val << ", f" << flag << ") ";
            }
            Rcpp::Rcout << std::endl;
        }
        Rcpp::Rcout << std::endl;
        for (std::size_t i = 0; i < ncol; i++)
        {
            VecIterator& citer = col_iter(i);
            Rcpp::Rcout << "Column iter " << i << ": ";
            for (; citer; ++citer)
            {
                double val = 0.0;
                int flag = static_cast<int>(citer.value(val));
                Rcpp::Rcout << "(" << val << ", f" << flag << ") ";
            }
            Rcpp::Rcout << std::endl;
        }
    }

public:
    virtual ~DataLoader() {}
    // Iterator for a row
    virtual VecIterator& row_iter(std::size_t index) = 0;
    // Iterator for a column
    virtual VecIterator& col_iter(std::size_t index) = 0;
    // Print a short description of the data loader
    virtual void description() const = 0;
    // Print verbose debugging information
    virtual void debug_info() {}
};

// Dense data loader
template <typename T>
class DenseDataLoader: public DataLoader
{
private:
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using ConstGenericMatrix = const Eigen::Ref<const Matrix>;

    ConstGenericMatrix m_mat;
    Matrix m_tmat;
    DenseVecIterator<T> m_rowiter;
    DenseVecIterator<T> m_coliter;

public:
    DenseDataLoader(ConstGenericMatrix& mat) :
        m_mat(mat), m_tmat(mat.transpose()),
        m_rowiter(nullptr, mat.cols()), m_coliter(nullptr, mat.rows())
    {}

    VecIterator& row_iter(std::size_t index) override
    {
        m_rowiter = DenseVecIterator<T>(&m_tmat.coeffRef(0, index), m_tmat.rows());
        return m_rowiter;
    }

    VecIterator& col_iter(std::size_t index) override
    {
        m_coliter = DenseVecIterator<T>(&m_mat.coeffRef(0, index), m_mat.rows());
        return m_coliter;
    }

    void description() const override
    {
        std::string element_type = std::is_same<T, int>::value ? "int" : "double";
        Rcpp::Rcout << "Dense data loader: " <<
            m_mat.rows() << " x " << m_mat.cols() <<
            ", " << element_type << " type" << std::endl;
    }

    void debug_info() override
    {
        this->description();
        Rcpp::Rcout << std::endl;
        this->print_data_iter(m_mat.rows(), m_mat.cols());
    }
};

// Sparse data loader
class SparseDataLoader: public DataLoader
{
private:
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using ConstGenericSparseMatrix = const Eigen::Ref<const SparseMatrix>;

    ConstGenericSparseMatrix m_mat;
    SparseMatrix m_tmat;
    SparseVecIterator m_rowiter;
    SparseVecIterator m_coliter;

public:
    SparseDataLoader(ConstGenericSparseMatrix& mat) :
        m_mat(mat), m_tmat(mat.transpose()),
        m_rowiter(nullptr, nullptr, 0, mat.cols()),
        m_coliter(nullptr, nullptr, 0, mat.rows())
    {}

    VecIterator& row_iter(std::size_t index) override
    {
        const int outer_start = m_tmat.outerIndexPtr()[index];
        const int nnz = m_tmat.outerIndexPtr()[index + 1] - outer_start;
        m_rowiter = SparseVecIterator(
            m_tmat.valuePtr() + outer_start, m_tmat.innerIndexPtr() + outer_start,
            nnz, m_tmat.rows()
        );
        return m_rowiter;
    }

    VecIterator& col_iter(std::size_t index) override
    {
        const int outer_start = m_mat.outerIndexPtr()[index];
        const int nnz = m_mat.outerIndexPtr()[index + 1] - outer_start;
        m_coliter = SparseVecIterator(
            m_mat.valuePtr() + outer_start, m_mat.innerIndexPtr() + outer_start,
            nnz, m_mat.rows()
        );
        return m_coliter;
    }

    void description() const override
    {
        Rcpp::Rcout << "Sparse data loader: " <<
            m_mat.rows() << " x " << m_mat.cols() <<
            ", double type" << std::endl;
    }

    void debug_info() override
    {
        this->description();
        Rcpp::Rcout << std::endl;
        this->print_data_iter(m_mat.rows(), m_mat.cols());
    }
};

#endif  // ESVD2_DATA_LOADER_H
