#include "pch.h"
#include "legendre.h"

Legendre::~Legendre()
{
}

Legendre& Legendre::normalize(const MatrixXXd& X, const MatrixXXd& Y)
{
	m_X = (-1 + 2 * (X.array() - X.minCoeff()) / (X.maxCoeff() - X.minCoeff())).matrix();
	m_Y = (-1 + 2 * (Y.array() - Y.minCoeff()) / (Y.maxCoeff() - Y.minCoeff())).matrix();
	
	return *this;
}

Legendre& Legendre::normalize(const VectorXd& x)
{
	m_x = (-1 + 2 * (x.array() - x.minCoeff()) / (x.maxCoeff() - x.minCoeff())).matrix();

	return *this;
}

std::vector<VectorXd> Legendre::gen_1d_poly(vec_i& orders)
{
	// sort the orders
	std::stable_sort(orders.begin(), orders.end());

	// grab the largest order
}
