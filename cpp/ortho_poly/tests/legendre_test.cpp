#include "pch.h"
#include "../ortho_poly/common.h"
#include "../ortho_poly/legendre.h"
#include "../ortho_poly/matrix_io.h"

TEST(Legendre, normalize)
{
	// load data
	int rows = 0, cols = 0;
	double* X = nullptr;
	double* Y = nullptr;
	double* Z = nullptr;

	read_matrix_from_disk("../../../data/X_legendre.bin", &rows, &cols, &X);
	read_matrix_from_disk("../../../data/Y_legendre.bin", &rows, &cols, &Y);
	read_matrix_from_disk("../../../data/Z_legendre.bin", &rows, &cols, &Z);

	Eigen::Map<MatrixXXd> Xmap(X, rows, cols);
	Eigen::Map<MatrixXXd> Ymap(Y, rows, cols);
	Eigen::Map<MatrixXXd> Zmap(Z, rows, cols);

	Legendre lg;
	lg.normalize(Xmap, Ymap);

	write_matrix_to_disk("../../../data/X_lg.bin", rows, cols, lg.GetX().data());
	write_matrix_to_disk("../../../data/Y_lg.bin", rows, cols, lg.GetY().data());

	free(X); X = nullptr;
	free(Y); Y = nullptr;
	free(Z); Z = nullptr;
}

TEST(Legendre, gen_1d_p)
{
	// without normalization
	std::cout << "Calculation only" << std::endl;
	Legendre lg;
	VectorXd x(3), x_norm(3);
	x << -1, 0, 1;
	x_norm << -2, 0, 2;
	auto polys = lg(x).gen_1d_p(vec_i{ 2, 1, 0, 1, 4 });
	for (const auto& i : polys) {
		std::cout << i.transpose() << std::endl;
	}

	// with normalization
	std::cout << "Normalization then calculation" << std::endl;
	polys = lg.normalize(x_norm).gen_1d_p(vec_i{ 0, 0, 1, 1, 2, 3, 4, 5 });
	for (const auto& i : polys) {
		std::cout << i.transpose() << std::endl;
	}
}

TEST(Legendre, gen_1d_p_dp)
{
	// without normalization
	std::cout << "Calculation only" << std::endl;
	Legendre lg;
	VectorXd x(3), x_norm(3);
	x << -1, 0, 1;
	x_norm << -2, 0, 2;
	auto [p, dp] = lg(x).gen_1d_p_dp(vec_i{ 0, 1, 2 });
	for (const auto& i : p) {
		std::cout << i.transpose() << std::endl;
	}
	for (const auto& i : dp) {
		std::cout << i.transpose() << std::endl;
	}
}

TEST(Legendre, gen_2d_p)
{
	// without normalization
	std::cout << "Calculation only" << std::endl;
	Legendre lg;
	MatrixXXd X(3, 3), Y(3, 3);
	X << -1, 0, 1,
		-1, 0, 1,
		-1, 0, 1;
	Y << 1, 1, 1,
		0, 0, 0,
		-1, -1, -1;
	vec_m polys = lg(X, Y).gen_2d_p(set_i({ 1, 2, 3 }));

	// with combined polynomial
	auto [P, Ps] = lg(X, Y).gen_2d_p({
		{ 1, 1},
		{ 2, 1 },
		{ 3, 1 },
		}
	);
	std::cout << "Individual Polynomials: " << std::endl;
	for (const auto& i : polys) {
		std::cout << i << std::endl;
	}
	std::cout << "Combined Polynomial: " << std::endl;
	std::cout << P << std::endl;


	//// with normalization
	//std::cout << "Normalization then calculation" << std::endl;
	//polys = lg.normalize(x_norm).gen_1d_p(vec_i{ 0, 0, 1, 1, 2, 3, 4, 5 });
	//for (const auto& i : polys) {
	//	std::cout << i.transpose() << std::endl;
	//}
}