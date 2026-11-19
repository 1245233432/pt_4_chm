#pragma once
#include "integration_scheme.h"
#include <functional>

namespace com_methods {
	class integration_scheme_interval : public integration_scheme {
	public:
		integration_scheme_interval(integration_scheme_types type); // конструктор
		double calculate_integral(const Point& begin, 
			const Point& end, 
			int num_segments, 
			const std::function<double(const Point& p)>& func) const; // метод для вычисления интеграла
																	  // begin, ends - начало, конец отрезка
																	  // num_segments - число сегментов
																	  // func - подынтегральная f(x)
		double runge_estimate(const double& I_h, const double& I_h2, const int& k) const; // метод для оценки погрешности по правилу Рунге
		double richardson_refinement(const double& I_h, const double& I_h2, const int& k) const; // уточнение по Ричардсону
		bool is_runge_applicable(const double& I_2h, const double& I_h, const double& I_h2, const int& k) const; // проверка на применимость уточнения по Ричардсону
		double estimate_order(const double& I_h, const double& I_h2, const double& I_true) const; // оценка порядка аппроксимации метода относительно шага сеточного разбиения
	};
}