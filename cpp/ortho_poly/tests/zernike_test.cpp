#include "pch.h"
#include "../ortho_poly/common.h"
#include "../ortho_poly/zernike.h"
#include "../ortho_poly/matrix_io.h"

TEST(common, cart2pol)
{
	MatrixXXd X(3, 3), Y(3, 3);
	X << -1, 0, 1,
		-1, 0, 1,
		-1, 0, 1;
	Y << 1, 1, 1,
		0, 0, 0,
		-1, -1, -1;

	auto [R, Theta] = cart2pol(X, Y);

	std::cout << R << std::endl;
	std::cout << Theta << std::endl;
}

TEST(Zernike, gen_2d_p)
{
	Zernike zer;
	MatrixXXd X(3, 3), Y(3, 3);
	X << -1, 0, 1,
		-1, 0, 1,
		-1, 0, 1;
	Y << 1, 1, 1,
		0, 0, 0,
		-1, -1, -1;
	vec_m polys = zer(X, Y).gen_2d_p(set_i({ 1, 2, 3, 4 }));

	// with combined polynomial
	auto [P, Ps] = zer(X, Y).gen_2d_p({
		{ 1, 1},
		{ 2, 1 },
		{ 3, 1 },
		{ 4, 1 },
		}
	);
	std::cout << "Individual Polynomials: " << std::endl;
	for (const auto& i : polys) {
		std::cout << i << std::endl;
	}
	std::cout << "Combined Polynomial: " << std::endl;
	std::cout << P << std::endl;
}

TEST(Zernike, fitting_2d)
{
	// load data
	int rows = 0, cols = 0;
	double* X = nullptr;
	double* Y = nullptr;
	double* Z = nullptr;

	read_matrix_from_disk("../../../data/X_zernike.bin", &rows, &cols, &X);
	read_matrix_from_disk("../../../data/Y_zernike.bin", &rows, &cols, &Y);
	read_matrix_from_disk("../../../data/Z_zernike.bin", &rows, &cols, &Z);

	Eigen::Map<MatrixXXd> Xmap(X, rows, cols);
	Eigen::Map<MatrixXXd> Ymap(Y, rows, cols);
	Eigen::Map<MatrixXXd> Zmap(Z, rows, cols);

	Zernike zer;
	auto j_orders = set_i{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
	zer.normalize(Xmap, Ymap).fit(Zmap, j_orders);
	auto coeffs = zer.coeffs();
	for (const auto& jc : coeffs) {
		std::cout << jc.first << ", " << jc.second << std::endl;
	}

	auto Zfit = zer.normalize(Xmap, Ymap).fit_predict(Zmap, j_orders);

	write_matrix_to_disk("../../../data/Zfit_zernike.bin", Zfit.rows(), Zfit.cols(), Zfit.data());

	free(X); X = nullptr;
	free(Y); Y = nullptr;
	free(Z); Z = nullptr;
}