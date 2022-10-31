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

vec_v Legendre::gen_1d_p(vec_i& orders)
{
	// base cases
	VectorXd p0(VectorXd::Ones(m_x.size()));
	VectorXd p1(m_x);

	// map to store the (order, poly) pairs to avoid re-calculation
	map_iv order_poly_map;

	// fill in the (order, poly) map with blank vectors
	for (const auto& i : orders) {
		if (0 == i) order_poly_map[i] = p0;
		if (1 == i) order_poly_map[i] = p1;
		order_poly_map[i];
	}

	// get the highest order
	auto max_order = order_poly_map.crbegin()->first;

	// reserve the space for the output
	vec_v polys;
	polys.reserve(orders.size());

	// calculate for the highest order & fill the order_poly_map with the lower ones
	VectorXd p2 = VectorXd::Zero(m_x.size());

	// recursive generation
	for (auto i = 2; i <= max_order; i++) {
		double n = double(i - 1);
		p2 = (2 * n + 1) / (n + 1) * m_x.cwiseProduct(p1) - n / (n + 1) * p0;

		// save p2 if required
		if (order_poly_map.find(i) != order_poly_map.end()) {
			order_poly_map[i] = p2;
		}

		// switch p0 and p1
		p0 = p1;
		p1 = p2;
	}


	// build the polynomial vectors
	for (const auto& i : orders) {
		polys.push_back(order_poly_map[i]);
	}

	return polys;
}

std::tuple<vec_v, vec_v> Legendre::gen_1d_p_dp(vec_i& orders)
{
	// base cases for p
	VectorXd p0(VectorXd::Ones(m_x.size()));
	VectorXd p1(m_x);

	// base cases for dp
	VectorXd dp0(VectorXd::Zero(m_x.size()));
	VectorXd dp1(VectorXd::Ones(m_x.size()));

	// map to store the (order, poly) pairs to avoid re-calculation
	map_iv order_p_map, order_dp_map;

	// fill in the (order, poly) map with blank vectors
	for (const auto& i : orders) {
		if (0 == i) {
			order_p_map[i] = p0; order_dp_map[i] = dp0;
		}
		if (1 == i) {
			order_p_map[i] = p1; order_dp_map[i] = dp1;
		}
		order_p_map[i];
		order_dp_map[i];
	}

	// get the highest order
	auto max_order = order_p_map.crbegin()->first;

	// reserve the space for the output
	vec_v p, dp;
	p.reserve(orders.size());
	dp.reserve(orders.size());

	// calculate for the highest order & fill the order_poly_map with the lower ones
	VectorXd p2 = VectorXd::Zero(m_x.size());
	VectorXd dp2 = VectorXd::Zero(m_x.size());

	// recursive generation
	for (auto i = 2; i <= max_order; i++) {
		double n = double(i - 1);
		p2 = (2 * n + 1) / (n + 1) * m_x.cwiseProduct(p1) - n / (n + 1) * p0;
		dp2 = (2 * n + 1) / (n + 1) * (p1 + m_x.cwiseProduct(dp1)) - n / (n + 1) * dp0;

		// save p2 if required
		if (order_p_map.find(i) != order_p_map.end()) {
			order_p_map[i] = p2;
		}
		if (order_dp_map.find(i) != order_dp_map.end()) {
			order_dp_map[i] = dp2;
		}

		// switch p0 and p1
		p0 = p1;
		p1 = p2;

		dp0 = dp1;
		dp1 = dp2;
	}

	// build the polynomial vectors
	for (const auto& i : orders) {
		p.push_back(order_p_map[i]);
		dp.push_back(order_dp_map[i]);
	}

	return std::make_tuple(p, dp);
}

std::tuple<MatrixXXd, vec_m> Legendre::gen_2d_p(const map_id& jorder_coeff)
{
	// obtain the orders for x & y, respectively
	vec_i x_orders, y_orders;
	vec_d coeffs;
	x_orders.reserve(jorder_coeff.size());
	y_orders.reserve(jorder_coeff.size());
	coeffs.reserve(jorder_coeff.size());

	for (const auto& jc : jorder_coeff) {
		auto b = int_t(ceil(sqrt(double(jc.first))));
		auto a = b * b - jc.first + 1;

		auto nsm = -a / 2 * ((a % 2) == 0 ? 1 : 0) + (a - 1) / 2 * (a % 2);
		auto nam = 2 * b - abs(nsm) - 2;

		x_orders.push_back((nam + nsm) / 2);
		y_orders.push_back((nam - nsm) / 2);

		/*int_t a = int_t(floor((1 + sqrt(1 + 8.0 * (jc.first - 1))) * 0.5));
		int_t b = (jc.first - 1) - (a * (a - 1)) / 2 + 1;

		y_orders.push_back(b - 1);
		x_orders.push_back(a - b);*/

		coeffs.push_back(jc.second);
	}

	// obtain the 1d polynomials in x and y, respectively
	auto p_x = (*this)(m_X.row(0)).gen_1d_p(x_orders);
	auto p_y = (*this)(m_Y.col(0)).gen_1d_p(y_orders);

	// get the 2d polynomials: multiply each of the 
	// x & y corresponding polynomials and coefficients
	vec_m Ps(jorder_coeff.size(), MatrixXXd());
	MatrixXXd P = MatrixXXd::Zero(m_X.rows(), m_X.cols());
	Ps.reserve(jorder_coeff.size());
	for (auto i = 0; i < p_x.size(); i++) {
		Ps[i] = p_y[i] * p_x[i].transpose();
		P += coeffs[i] * Ps[i];
	}

	return std::make_tuple(P, Ps);
}

vec_m Legendre::gen_2d_p(const set_i& j_orders)
{
	// obtain the orders for x & y, respectively
	vec_i x_orders, y_orders;
	x_orders.reserve(j_orders.size());
	y_orders.reserve(j_orders.size());

	for (const auto& j : j_orders) {
		auto b = int_t(ceil(sqrt(double(j))));
		auto a = b * b - j + 1;

		auto nsm = -a / 2 * ((a % 2) == 0 ? 1 : 0) + (a - 1) / 2 * (a % 2);
		auto nam = 2 * b - abs(nsm) - 2;

		x_orders.push_back((nam + nsm) / 2);
		y_orders.push_back((nam - nsm) / 2);
		/*int_t a = int_t(floor((1 + sqrt(1 + 8.0 * (j - 1))) * 0.5));
		int_t b = (j - 1) - (a * (a - 1)) / 2 + 1;

		y_orders.push_back(b - 1);
		x_orders.push_back(a - b);*/
	}

	// obtain the 1d polynomials in x and y, respectively
	auto p_x = (*this)(m_X.row(0)).gen_1d_p(x_orders);
	auto p_y = (*this)(m_Y.col(0)).gen_1d_p(y_orders);

	// get the 2d polynomials: multiply each of the x & y corresponding polynomials
	vec_m Ps;
	Ps.reserve(j_orders.size());
	for (auto i = 0; i < p_x.size(); i++) {
		Ps.push_back(p_y[i] * p_x[i].transpose());
	}

	// output the 2d polynomials
	return Ps;
}
