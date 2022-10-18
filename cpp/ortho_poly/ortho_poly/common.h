#ifndef COMMON_H
#define COMMON_H

#ifdef ORTHOPOLY_EXPORTS
#define ORTHOPOLY_API __declspec(dllexport)
#else
#define ORTHOPOLY_API __declspec(dllimport)
#endif

// Aliases
using int_t = Eigen::Index;
using Tripletd = Eigen::Triplet<double>;
using TripletListd = std::vector<Tripletd>;
using SparseMatrixXXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using MatrixXXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorXd = Eigen::VectorXd;
using VectorXi = Eigen::VectorXi;
using VectorMapd = Eigen::Map<VectorXd, Eigen::Unaligned>;
using std_vecd = std::vector<double>;
using Solver = Eigen::LeastSquaresConjugateGradient<SparseMatrixXXd>;
using QRSolver = Eigen::SparseQR<SparseMatrixXXd, Eigen::COLAMDOrdering<int>>;

inline int_t ID_1D(int_t x, int_t y, int_t width) { return (y * width + x); }


#endif // !COMMON_H
