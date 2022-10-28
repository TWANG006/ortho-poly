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

	// build and solve the Ax = b system
	auto coeffs = _build_solve_Axb(Ps, Z);

	// update the (j, c) map member
	m_order_coeff_map.clear();
	int_t id = 0;
	for (const auto& j : j_orders) {
		m_order_coeff_map[j] = coeffs[id++];
	}

	return *this;
}

VectorXd BasePolynomial::_build_solve_Axb(const vec_m& Ps, const MatrixXXd& Z)
{
	// reserve A and b to the largest sizes
	vec_d A, b;
	A.reserve(Z.size() * Ps.size());
	b.reserve(Z.size());

	// build A and b
	for (auto i = 0; i < Z.rows(); i++) {
		for (auto j = 0; j < Z.cols(); j++) {
			if (isfinite(Z(i, j))) {
				// build one row of A
				for (const auto& P : Ps) {
					A.push_back(P(i, j));
				}
				// build one element of b
				b.push_back(Z(i, j));
			}
		}
	}

	// map the A and b
	auto m = b.size();
	auto n = A.size() / m;
	MatrixMapd Amap(A.data(), m, n);
	VectorMapd bmap(b.data(), m, n);

	// solve the Ax = b
	//return Amap.colPivHouseholderQr().solve(bmap);
	return Amap.fullPivHouseholderQr().solve(bmap);
}
