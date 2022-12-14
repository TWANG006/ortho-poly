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

BasePolynomial& BasePolynomial::normalize(const MatrixXXd& X, const MatrixXXd& Y)
{
	m_X = (-1 + 2 * (X.array() - X.minCoeff()) / (X.maxCoeff() - X.minCoeff())).matrix();
	m_Y = (-1 + 2 * (Y.array() - Y.minCoeff()) / (Y.maxCoeff() - Y.minCoeff())).matrix();

	return *this;
}

BasePolynomial& BasePolynomial::normalize(const VectorXd& x)
{
	m_x = (-1 + 2 * (x.array() - x.minCoeff()) / (x.maxCoeff() - x.minCoeff())).matrix();

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

MatrixXXd BasePolynomial::predict()
{
	auto [Zfit, Ps] = gen_2d_p(m_order_coeff_map);

	for (auto i = 0; i < Zfit.rows(); i++) {
		for (auto j = 0; j < Zfit.cols(); j++) {
			if (m_is_valid_id[ID_1D(j, i, Zfit.cols())] == false) {
				Zfit(i, j) = NAN;
			}
		}
	}

	return Zfit;
}

MatrixXXd BasePolynomial::predict(const map_id& order_coeff_map)
{
	auto [Zfit, Ps] = gen_2d_p(order_coeff_map);

	for (auto i = 0; i < Zfit.rows(); i++) {
		for (auto j = 0; j < Zfit.cols(); j++) {
			if (m_is_valid_id[ID_1D(j, i, Zfit.cols())] == false) {
				Zfit(i, j) = NAN;
			}
		}
	}

	return Zfit;
}

MatrixXXd BasePolynomial::fit_predict(const MatrixXXd& Z, const set_i& j_orders)
{
	return fit(Z, j_orders).predict();
}

VectorXd BasePolynomial::_build_solve_Axb(const vec_m& Ps, const MatrixXXd& Z)
{
	// reserve A and b to the largest sizes
	vec_d A, b;
	A.reserve(Z.size() * Ps.size());
	b.reserve(Z.size());
	m_is_valid_id.resize(Z.size());

	// build A and b
	for (auto i = 0; i < Z.rows(); i++) {
		for (auto j = 0; j < Z.cols(); j++) {
			if (isfinite(Z(i, j))) {
				m_is_valid_id[ID_1D(j, i, Z.cols())] = true;
				// build one row of A
				for (const auto& P : Ps) {
					A.push_back(P(i, j));
				}
				// build one element of b
				b.push_back(Z(i, j));
			}
			else {
				m_is_valid_id[ID_1D(j, i, Z.cols())] = false;
			}
		}
	}

	// map the A and b
	auto m = b.size();
	auto n = A.size() / m;
	MatrixMapd Amap(A.data(), m, n);
	VectorMapd bmap(b.data(), m);

	// solve the Ax = b
	return Amap.colPivHouseholderQr().solve(bmap);
}
