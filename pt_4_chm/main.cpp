#include <iostream>
#include <iomanip>
#include <functional>
#include "point.h"
#include "integration_scheme_interval.h"

int main() {
	setlocale(LC_ALL, "rus");

	std::function<double(const com_methods::Point& p)> phi = [](const com_methods::Point& p) {
		return exp(p.x()) * cos(p.x());
	};

	const double I_true = 0.8309116713772953; // аналитическое значение интеграла от fi

	com_methods::Point begin(0.07, 0.0, 0.0);
	com_methods::Point end(0.69, 0.0, 0.0);

	const int N = 1; // число сегментов

	com_methods::integration_scheme::integration_scheme_types schemes[] = { // выбранные схемы интегрирования
		com_methods::integration_scheme::parab,
		com_methods::integration_scheme::gauss3 };

	const char* names[] = { "Parab", "Gauss-3" };

	for (int s = 0; s < 2; s++) {
		com_methods::integration_scheme_interval quadrature(schemes[s]);
		int k = quadrature.get_order(schemes[s]);

		double I_h = quadrature.calculate_integral(begin, end, N, phi);
		double I_h2 = quadrature.calculate_integral(begin, end, N * 2, phi);
		double I_h4 = quadrature.calculate_integral(begin, end, N * 4, phi);

		double err_h = I_true - I_h;
		double err_h2 = I_true - I_h2;
		double runge = quadrature.runge_estimate(I_h, I_h2, k);
		double I_R = quadrature.richardson_refinement(I_h, I_h2, k);
		double err_R = I_true - I_R;
		double appro_estimation_h2 = quadrature.estimate_order(I_h, I_h2, I_true);
		double appro_estimation_h4 = quadrature.estimate_order(I_h2, I_h4, I_true);

		std::cout << std::scientific << std::setprecision(10);
		std::cout << "Метод " << names[s] << std::endl;
		std::cout << "h = " << (0.69 - 0.07) / N << std::endl;
		std::cout << "I_h = " << I_h << std::endl;
		std::cout << "I_h2 = " << I_h2 << std::endl;
		std::cout << "I_h4 = " << I_h4 << std::endl;
		std::cout << "|I_true - I_h| = " << err_h << std::endl;
		std::cout << "Погрешность по Рунге: " << runge << std::endl;
		std::cout << "Уточнение по Ричардсону I_R = "<< I_R << std::endl;
		std::cout << "|I_true - I_R| = " << err_R << std::endl;
		std::cout << "Оценка порядка аппроксимации k при шаге h/2: " << appro_estimation_h2 << std::endl;
		std::cout << "Оценка порядка аппроксимации k при шаге h/4: " << appro_estimation_h4 << std::endl;
		std::cout << std::endl;

	}

	com_methods::integration_scheme_interval quadrature(com_methods::integration_scheme::gauss3);
	double I_h4 = quadrature.calculate_integral(begin, end, N * 4, phi);
	double I_h8 = quadrature.calculate_integral(begin, end, N * 8, phi);
	double appro_estimation_h8 = quadrature.estimate_order(I_h4, I_h8, I_true);
	std::cout << "I_h8 = " << I_h8 << std::endl;
	std::cout << "|I_true - I_h8| = " << I_true - I_h8 << std::endl;
	std::cout << "Оценка порядка аппроксимации k при шаге h/8: " << appro_estimation_h8 << std::endl << std::endl;

	double I_h16 = quadrature.calculate_integral(begin, end, N * 16, phi);
	double appro_estimation_h16 = quadrature.estimate_order(I_h8, I_h16, I_true);
	std::cout << "I_h16 = " << I_h16 << std::endl;
	std::cout << "|I_true - I_h16| = " << I_true - I_h16 << std::endl;
	std::cout << "Оценка порядка аппроксимации k при шаге h/16: " << appro_estimation_h16 << std::endl << std::endl;


	double I_h32 = quadrature.calculate_integral(begin, end, N * 32, phi);
	double appro_estimation_h32 = quadrature.estimate_order(I_h16, I_h32, I_true);
	std::cout << "I_h32 = " << I_h32 << std::endl;
	std::cout << "|I_true - I_h32| = " << I_true - I_h32 << std::endl;
	std::cout << "Оценка порядка аппроксимации k при шаге h/32: " << appro_estimation_h32 << std::endl << std::endl;

	return 0;
}