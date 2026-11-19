#pragma once
#include "point.h"
#include <vector>

namespace com_methods {
	class integration_scheme {
	protected:
		std::vector<Point> points; // узлы интегрирования
		std::vector<double> weights; // веса квадратурной формулы (альфы)

	public:
		enum integration_scheme_types { // типы квадратурных формул
			gauss1 = 1, // он же метод прямоугольников
			gauss2,
			gauss3,
			gauss4,
			gauss5,
			trap,
			parab
		};
		int get_order(integration_scheme_types type) const { // получение порядка аппроксимации для правила Рунге
			switch (points.size()) {
				case trap: return 2;
				case parab: return 4;
				case gauss1: return 2;
				case gauss2: return 4;
				case gauss3: return 6;
				case gauss4: return 8;
				case gauss5: return 10;
				default: return 2;
			}
		}
	};
}