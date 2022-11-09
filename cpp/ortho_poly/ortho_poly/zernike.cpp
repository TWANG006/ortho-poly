#include "pch.h"
#include "zernike.h"

Zernike::~Zernike()
{
}

vec_m Zernike::gen_2d_p(const set_i& j_orders)
{
	// convert to polar coordinates

	// obtain the m, n orders, respectively
	for (const auto& j : j_orders) {
		auto [m_order, n_order] = _mn_from_j(j);
		
	}

	return vec_m();
}

std::tuple<MatrixXXd, vec_m> Zernike::gen_2d_p(const map_id& jorder_coeff)
{
	return std::tuple<MatrixXXd, vec_m>();
}

std::tuple<int_t, int_t> Zernike::_mn_from_j(const int& j)
{
	auto b = int_t(ceil(sqrt(double(j))));
	auto a = b * b - j + 1;

	auto m = -a / 2 * ((a % 2) == 0 ? 1 : 0) + (a - 1) / 2 * (a % 2);
	auto n = 2 * b - abs(m) - 2;

	return std::make_tuple(m, n);
}


