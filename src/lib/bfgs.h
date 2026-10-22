/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

/* Modified for FastDock Vina (selected D): column-major symmetric BFGS storage
   and two-double CPU SIMD, preserving each output accumulation order.
   See docs/PROVENANCE.md for upstream and experimental source identities. */

#ifndef VINA_BFGS_H
#define VINA_BFGS_H

#include "matrix.h"
#include <cstring>

// Only the local BFGS inverse Hessian uses this full symmetric storage.
class flmat {
	flv m_data;
	sz m_dim;
public:
	flmat() : m_dim(0) {}
	flmat(sz n, fl filler) : m_data(n * n, filler), m_dim(n) {}
	sz dim() const { return m_dim; }
	sz index(sz i, sz j) const { return i + j * m_dim; }
	sz index_permissive(sz i, sz j) const { return index(i, j); }
	const fl& operator()(sz i) const { return m_data[i]; }
	fl& operator()(sz i) { return m_data[i]; }
	const fl& operator()(sz i, sz j) const { return m_data[index(i, j)]; }
	fl& operator()(sz i, sz j) { return m_data[index(i, j)]; }
};

typedef fl bfgs_lane2 __attribute__((vector_size(16)));

template<typename Change>
void minus_mat_vec_product(const flmat& m, const Change& in, Change& out) {
	const sz n = m.dim();
	flv values(n);
	VINA_FOR(j, n) values[j] = in(j);
	sz i = 0;
	for(; i + 1 < n; i += 2) {
		bfgs_lane2 sums = {0, 0};
		VINA_FOR(j, n) {
			bfgs_lane2 coefficients;
			std::memcpy(&coefficients, &m(i + j * n), sizeof(coefficients));
			const bfgs_lane2 input = {values[j], values[j]};
			sums += coefficients * input;
		}
		fl lanes[2];
		std::memcpy(lanes, &sums, sizeof(sums));
		out(i) = -lanes[0];
		out(i + 1) = -lanes[1];
	}
	for(; i < n; ++i) {
		fl sum = 0;
		VINA_FOR(j, n)
			sum += m(i + j * n) * values[j];
		out(i) = -sum;
	}
}

template<typename Change>
inline fl scalar_product(const Change& a, const Change& b, sz n) {
	fl tmp = 0;
	VINA_FOR(i, n)
		tmp += a(i) * b(i);
	return tmp;
}

template<typename Change>
inline bool bfgs_update(flmat& h, const Change& p, const Change& y, const fl alpha) {
	const fl yp  = scalar_product(y, p, h.dim());
	if(alpha * yp < epsilon_fl) return false; // FIXME?
	Change minus_hy(y); minus_mat_vec_product(h, y, minus_hy);
	const fl yhy = - scalar_product(y, minus_hy, h.dim());
	const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
	const sz n = p.num_floats();
	flv p_values(n), hy_values(n);
	VINA_FOR(i, n) {
		p_values[i] = p(i);
		hy_values[i] = minus_hy(i);
	}
	VINA_FOR(j, n)
		for(sz i = 0; i <= j; ++i) {
			h(i, j) +=   alpha * r * (hy_values[i] * p_values[j]
	                                + hy_values[j] * p_values[i]) +
			           + alpha * alpha * (r*r * yhy  + r) * p_values[i] * p_values[j]; // s * s == alpha * alpha * p * p
			h(j, i) = h(i, j);
		}
	return true;
}

template<typename F, typename Conf, typename Change>
fl line_search(F& f, sz n, const Conf& x, const Change& g, const fl f0, const Change& p, Conf& x_new, Change& g_new, fl& f1, int& evalcount) { // returns alpha
	const fl c0 = 0.0001;
	const unsigned max_trials = 10;
	const fl multiplier = 0.5;
	fl alpha = 1;

	const fl pg = scalar_product(p, g, n);

	VINA_U_FOR(trial, max_trials) {
		x_new = x; x_new.increment(p, alpha);
		f1 = f(x_new, g_new);
		evalcount++;
		if(f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
			break;
		alpha *= multiplier;
	}
	return alpha;
}

inline void set_diagonal(flmat& m, fl x) {
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

template<typename Change>
void subtract_change(Change& b, const Change& a, sz n) { // b -= a
	VINA_FOR(i, n)
		b(i) -= a(i);
}

template<typename F, typename Conf, typename Change>
fl bfgs(F& f, Conf& x, Change& g, const unsigned max_steps, const fl average_required_improvement, const sz over,
		int& evalcount) { // x is I/O, final value is returned
	sz n = g.num_floats();
	flmat h(n, 0);
	set_diagonal(h, 1);

	Change g_new(g);
	Conf x_new(x);
	fl f0 = f(x, g);
	evalcount++;

	fl f_orig = f0;
	Change g_orig(g);
	Conf x_orig(x);

	Change p(g);

	flv f_values; f_values.reserve(max_steps+1);
	f_values.push_back(f0);

	VINA_U_FOR(step, max_steps) {
		minus_mat_vec_product(h, g, p);
		fl f1 = 0;
		const fl alpha = line_search(f, n, x, g, f0, p, x_new, g_new, f1, evalcount);
		Change y(g_new); subtract_change(y, g, n);

		f_values.push_back(f1);
		f0 = f1;
		x = x_new;
		if(!(std::sqrt(scalar_product(g, g, n)) >= 1e-5)) break; // breaks for nans too // FIXME !!?? 
		g = g_new; // ?

		if(step == 0) {
			const fl yy = scalar_product(y, y, n);
			if(std::abs(yy) > epsilon_fl)
				set_diagonal(h, alpha * scalar_product(y, p, n) / yy);
		}

		bool h_updated = bfgs_update(h, p, y, alpha);
	}
	if(!(f0 <= f_orig)) { // succeeds for nans too
		f0 = f_orig;
		x = x_orig;
		g = g_orig;
	}
	return f0;
}

#endif
