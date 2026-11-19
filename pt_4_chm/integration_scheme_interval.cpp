#include "integration_scheme_interval.h"

namespace com_methods {
	integration_scheme_interval::integration_scheme_interval(integration_scheme_types type) {
		switch (type) {
			case gauss1: {
				weights = { 2 };
				points = { Point(0, 0, 0) };
				break;
			}
			case gauss2: {
				weights = { 1, 1 };
				points = { Point(-1.0 / sqrt(3.0), 0, 0), Point(1.0 / sqrt(3.0), 0, 0) };
				break;
			}
			case gauss3: {
				weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
				points = { Point(-sqrt(3.0 / 5.0), 0, 0), Point(0.0, 0.0, 0.0), Point(sqrt(3.0 / 5.0), 0, 0) };
				break;
			}
			case gauss4: {
				weights = { (18.0 - sqrt(30.0)) / 36.0,
					(18.0 + sqrt(30.0)) / 36.0,
					(18.0 + sqrt(30.0)) / 36.0,
					(18.0 - sqrt(30.0)) / 36.0 };

				points = { Point(-sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0), 0, 0),
					Point(-sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), 0, 0),
					Point(sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), 0, 0),
					Point(sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0), 0, 0) };

				break;
			}
			case gauss5: {
				weights = { (322.0 - 13.0 * sqrt(70.0)) / 900.0,
					(322.0 + 13.0 * sqrt(70.0)) / 900.0,
					128.0 / 225.0,
					(322.0 + 13.0 * sqrt(70.0)) / 900.0,
					(322.0 - 13.0 * sqrt(70.0)) / 900.0 };

				points = { Point(-sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0, 0.0, 0.0),
					Point(-sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0, 0.0, 0.0),
					Point(0.0, 0.0, 0.0),
					Point(sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0, 0.0, 0.0),
					Point(sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0, 0.0, 0.0) };

				break;
			}
			case trap: {
				weights = { 1, 1 };
				points = { Point(-1.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
				break;
			}
			case parab: {
				weights = {1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0 };
				points = { Point(-1.0, 0.0, 0.0), Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
				break;
			}
		}
	}

	double integration_scheme_interval::calculate_integral(const Point& begin, const Point& end, int num_segments, const std::function<double(const Point& p)>& func) const {
		double result = 0.0; // квадратурная сумма
		double x0; // начальная точка сегмента
		double h = (end.x() - begin.x()) / double(num_segments); // длина отрезка (сегмента)
		for (int i = 0; i < num_segments; i++) {
			x0 = begin.x() + i * h;
			for (int integ_point = 0; integ_point < points.size(); integ_point++) {
				auto p = Point(x0 + (1.0 + points[integ_point].x()) * h / 2.0, 0.0, 0.0);
				result += weights[integ_point] * func(p);
			}
		}
		return result * (h / 2.0);
	}

	double integration_scheme_interval::runge_estimate(const double& I_h, const double& I_h2, const int& k) const {
		return (I_h2 - I_h) / (pow(2.0, k) - 1.0);
	}

	double integration_scheme_interval::richardson_refinement(const double& I_h, const double& I_h2, const int& k) const {
		return I_h2 + (I_h2 - I_h) / (pow(2.0, k) - 1.0);
	}

	double integration_scheme_interval::estimate_order(const double& I_h, const double& I_h2, const double& I_true) const {
		double denominator = I_true - I_h2;
		if (fabs(denominator) < 1e-15) {
			return 0.0;
		}
		return log2(fabs(1.0 + (I_h2 - I_h) / denominator));
	}
}