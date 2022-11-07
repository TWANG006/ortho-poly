#include "pch.h"
#include "zernike.h"

Zernike::~Zernike()
{
}

vec_m Zernike::gen_2d_p(const set_i& j_orders)
{
	return vec_m();
}

std::tuple<MatrixXXd, vec_m> Zernike::gen_2d_p(const map_id& jorder_coeff)
{
	return std::tuple<MatrixXXd, vec_m>();
}

map_iv Zernike::_precompute_r_powers(const VectorXd& r, vec_i m, vec_i n)
{
	// obtain the unique powers
	set_i power_set;
	for (auto i = 0; i < n.size(); i++) {
		for (auto m_abs = abs(m[i]); m_abs <= n[i]; m_abs += 2) {
			power_set.insert(m_abs);
		}
	}

	// calculate the r_powers for the unique powers
	map_iv r_power_map;
	for (const auto& power : power_set) {
		r_power_map[power] = r.array().pow(power).matrix();
	}

	return r_power_map;
}

double Zernike::_factorial(const unsigned long long& n)
{
	
}

