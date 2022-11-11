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