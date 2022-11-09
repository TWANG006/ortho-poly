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
using ii = std::pair<int_t, int_t>;
using map_iim = std::map<ii, MatrixXXd>;

// functions
inline int_t ID_1D(int_t x, int_t y, int_t width) { return (y * width + x); }
inline std::tuple<MatrixXXd, MatrixXXd> cart2pol(const MatrixXXd& X, const MatrixXXd& Y)
{
	return std::make_tuple(
		(X.array().square() + Y.array().square()).sqrt().matrix(),
		Y.binaryExpr(X, [](double a, double b) {return std::atan2(a, b); })
	);
}

#endif // !COMMON_H

