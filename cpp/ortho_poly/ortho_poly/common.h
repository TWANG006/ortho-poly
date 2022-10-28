#ifndef COMMON_H
#define COMMON_H

#ifdef ORTHOPOLY_EXPORTS
#define ORTHOPOLY_API __declspec(dllexport)
#else
#define ORTHOPOLY_API __declspec(dllimport)
#endif

// Aliases
using int_t = Eigen::Index;
using MatrixXXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorXd = Eigen::VectorXd;
using VectorXi = Eigen::VectorXi;
using VectorMapd = Eigen::Map<VectorXd, Eigen::Unaligned>;
using MatrixMapd = Eigen::Map<MatrixXXd, Eigen::Unaligned>;
using vec_i = std::vector<int_t>;
using set_i = std::set<int_t>;
using vec_d = std::vector<double>;
using vec_v = std::vector<VectorXd>;
using vec_m = std::vector<MatrixXXd>;
using map_id = std::map<int_t, double>;
using map_iv = std::map<int_t, VectorXd>;
inline int_t ID_1D(int_t x, int_t y, int_t width) { return (y * width + x); }


#endif // !COMMON_H

