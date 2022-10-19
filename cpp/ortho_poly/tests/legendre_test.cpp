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