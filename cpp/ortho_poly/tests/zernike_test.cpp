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
	MatrixXXd X(3, 3), Y(3, 3);
	X << -1, 0, 1,
		-1, 0, 1,
		-1, 0, 1;
	Y << 1, 1, 1,
		0, 0, 0,
		-1, -1, -1;

	auto [R, Theta] = cart2pol(X, Y);

	Zernike zer;
	set_i j{ 1, 2, 3, 4 };
	auto Z = zer(X, Y).gen_2d_p(j);

	for (int i = 0; i < Z.size(); i++) {
		std::cout << Z[i] << std::endl << std::endl;
	}
}