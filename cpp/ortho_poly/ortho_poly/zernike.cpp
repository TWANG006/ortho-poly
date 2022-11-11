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
	vec_m Zs;
	Zs.reserve(j_orders.size());
	for (const auto& j : j_orders) {
		// convert j to m, n
		auto [m_order, n_order] = _mn_from_j(j);

		// calculate R_mn
		MatrixXXd R_mn = _calculate_Rmn(mn_p_map, R, R0, abs(m_order), n_order);
		if (m_order < 0) {
			Zs.push_back((R_mn.array() * (Theta.array() * (-m_order)).sin()).matrix());
		}
		else {
			Zs.push_back((R_mn.array() * (Theta.array() * (m_order)).cos()).matrix());
		}
	}

	return Zs;
}

std::tuple<MatrixXXd, vec_m> Zernike::gen_2d_p(const map_id& jorder_coeff)
{
	// convert to polar coordinates
	auto [R, Theta] = cart2pol(m_X, m_Y);

	// memorize the already calculated polynomials
	map_iim mn_p_map;
	mn_p_map[std::make_pair(0, 0)] = MatrixXXd::Ones(m_X.rows(), m_X.cols());
	MatrixXXd R0 = MatrixXXd::Zero(m_X.rows(), m_X.cols());

	// calculate the polynomials for j_orders
	vec_m Zs;
	MatrixXXd Z = MatrixXXd::Zero(m_X.rows(), m_X.cols());
	Zs.reserve(jorder_coeff.size());
	for (const auto& jc : jorder_coeff) {
		// convert j to m, n
		auto [m_order, n_order] = _mn_from_j(jc.first);

		// calculate R_mn
		MatrixXXd R_mn = _calculate_Rmn(mn_p_map, R, R0, abs(m_order), n_order);
		if (m_order < 0) {
			Zs.push_back((R_mn.array() * (Theta.array() * (-m_order)).sin()).matrix());
		}
		else {
			Zs.push_back((R_mn.array() * (Theta.array() * (m_order)).cos()).matrix());
		}

		// multiply the coefficients
		Z += jc.second * Zs.back();
	}

	return std::make_tuple(Z, Zs);
}

std::tuple<int_t, int_t> Zernike::_mn_from_j(const int& j)
{
	auto b = int_t(ceil(sqrt(double(j))));
	auto a = b * b - j + 1;

	auto m = -a / 2 * ((a % 2) == 0 ? 1 : 0) + (a - 1) / 2 * (a % 2);
	auto n = 2 * b - abs(m) - 2;

	return std::make_tuple(m, n);
}

MatrixXXd Zernike::_calculate_Rmn(map_iim& mn_p_map, const MatrixXXd& R, const MatrixXXd& R0, const int_t& m_order, const int_t& n_order)
{
	if (mn_p_map.find(std::make_pair(m_order, n_order)) != mn_p_map.end())
	{
		return mn_p_map[std::make_pair(m_order, n_order)];
	}
	else {
		MatrixXXd R_mn;
		for (auto m_plus_n = 2; m_plus_n <= m_order + n_order; m_plus_n += 2) {
			for (auto n = 1; n <= n_order; n++) {
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
							((n - 2 < m) ? R0 : mn_p_map[std::make_pair(m, n - 2)]).array()
							).matrix();
						mn_p_map[std::make_pair(m, n)] = R_mn;
					}
					else {
						R_mn = mn_p_map[std::make_pair(m, n)];
					}
				}
			}
		}
		return R_mn;
	}
}


