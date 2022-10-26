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

	// get the highest order
	auto max_order = orders.back();

	// map to store the (order, poly) pairs to avoid re-calculation
	map_iv order_poly_map;

	// fill in the (order, poly) map with blank vectors
	for (const auto& i : orders) {
		order_poly_map[i];
	}

	// reserve the space for the output
	std::vector<VectorXd> polys(orders.size(), VectorXd::Zero(m_x.size()));

	// calculate for the highest order & fill the order_poly_map with the lower ones
	VectorXd p0(VectorXd::Ones(m_x.size()));
	VectorXd p1(m_x);

	// recursive generation
	if (max_order == 0) {
		order_poly_map[max_order] = p0;
	}
	else if (max_order == 1) {
		order_poly_map[1] = p1;
	}
	else {
		for (auto i = 2; i <= max_order; i++) {
			auto n = i - 1;
			VectorXd p2 = (2 * n + 1) / (n + 1) * m_x.cwiseProduct(p1) - n / (n + 1) * p0;

			// save p2 if required
			if (order_poly_map.find(i) != order_poly_map.end()) {
				order_poly_map[i] = p2;
			}

			// switch p0 and p1
			p0 = p1;
			p1 = p2;
		}
	}
}
