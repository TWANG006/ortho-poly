#include "pch.h"
#include "basepolynomial.h"


BasePolynomial::~BasePolynomial()
{
}

BasePolynomial& BasePolynomial::operator()(const MatrixXXd& X, const MatrixXXd& Y)
{
	m_X = X;
	m_Y = Y;
	return *this;
}

BasePolynomial& BasePolynomial::operator()(const VectorXd& x)
{
	m_x = x;
	return *this;
}

BasePolynomial& BasePolynomial::fit(const MatrixXXd& Z, const set_i& j_orders)
{
	// generate the polynomials of the j_orders
	auto Ps = gen_2d_p(j_orders);

	// build the Ax = b system


	return *this;
}

std::tuple<MatrixXXd, VectorXd> BasePolynomial::_build_Axb(const vec_m& Ps, const MatrixXXd& Z)
{
	for (auto i = 0; i < Z.rows(); i++) {
		for (auto j = 0; j < Z.cols(); j++) {
			if (isfinite(Z)) {

			}
		}
	}

	return std::tuple<MatrixXXd, VectorXd>();
}
