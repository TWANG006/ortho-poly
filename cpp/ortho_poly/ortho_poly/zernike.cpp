#include "pch.h"
#include "zernike.h"

Zernike::~Zernike()
{
}

vec_m Zernike::gen_2d_p(const set_i& j_orders)
{
	// convert to polar coordinates
	auto [R, Theta] = cart2pol(m_X, m_Y);

	// memorize the already calculated polynomials
	map_iim mn_p_map;
	mn_p_map[std::make_pair(0, 0)] = MatrixXXd::Ones(m_X.rows(), m_X.cols());
	MatrixXXd R0 = MatrixXXd::Zero(m_X.rows(), m_X.cols());

	// calculate the polynomials for j_orders
	vec_m Rs;
	Rs.reserve(j_orders.size());
	for (const auto& j : j_orders) {
		// convert j to m, n
		auto [m_order, n_order] = _mn_from_j(j);

		// calculate R_mn
		MatrixXXd R_mn;
		for (auto m_plus_n = 2; m_plus_n < m_order + n_order; m_plus_n += 2) {
			for (auto n = 1; n < n_order; n++) {
				auto m = m_plus_n - n;
				if (m < 0) break;
				else if (n < m) continue;
				else {
					if (mn_p_map.find(std::make_pair(m, n)) == mn_p_map.end()) {
						R_mn = (
							R.array() * (
								((n - 1 < abs(m - 1)) ? R0 : mn_p_map[std::make_pair(abs(m - 1), n - 1)]) +
								((n - 1 < m + 1) ? R0 : mn_p_map[std::make_pair(m + 1, n - 1)])
								).array() -
							((n - 2 < m + 1) ? R0 : mn_p_map[std::make_pair(abs(m + 1), n - 2)]).array()
							).matrix();
						mn_p_map[std::make_pair(m, n)] = R_mn;
					}
				}
			}
		}
		Rs.push_back(R_mn);
	}

	return Rs;
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


