#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>


// Определение системы уравнений
std::vector<double> system(double x, const std::vector<double>& y,
	std::function<double(double)> p,
	std::function<double(double)> q,
	std::function<double(double)> f) {
	double y1 = y[0];
	double y2 = y[1];
	double dy1dx = y2;
	double dy2dx = f(x) - p(x) * y2 - q(x) * y1;
	return{ dy1dx, dy2dx };
}

std::vector<std::vector<double>> rungeKutta(double a, double b, const std::vector<double>& y0,
	std::function<double(double)> p,
	std::function<double(double)> q,
	std::function<double(double)> f,
	double h = 0.01) {
	int n = static_cast<int>((b - a) / h);
	std::vector<std::vector<double>> result(n + 1, std::vector<double>(2));
	double x = a;
	result[0] = y0;

	for (int i = 0; i < n; ++i) {
		std::vector<double> k1 = system(x, result[i], p, q, f);
		std::vector<double> y_temp = { result[i][0] + h * k1[0] / 2, result[i][1] + h * k1[1] / 2 };
		std::vector<double> k2 = system(x + h / 2, y_temp, p, q, f);
		y_temp = { result[i][0] + h * k2[0] / 2, result[i][1] + h * k2[1] / 2 };
		std::vector<double> k3 = system(x + h / 2, y_temp, p, q, f);
		y_temp = { result[i][0] + h * k3[0], result[i][1] + h * k3[1] };
		std::vector<double> k4 = system(x + h, y_temp, p, q, f);

		result[i + 1][0] = result[i][0] + h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
		result[i + 1][1] = result[i][1] + h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
		x += h;
	}

	return result;
}


// Определение метода стрельбы
std::vector<std::vector<double>> shootingMethod(double a, double b, double i1, double j1, double g1,
	double i2, double j2, double g2,
	std::function<double(double)> p,
	std::function<double(double)> q,
	std::function<double(double)> f,
	double initial_guess, double h = 1) {
	auto boundary_conditions = [&](double s) {
		double y_a = (g1-s*i1)/j1;  // Начальное значение y(a)
		std::vector<double> y0 = { y_a, s };
		auto sol = rungeKutta(a, b, y0, p, q, f, h);
		double y_b = sol.back()[0];
		double y_b_prime = sol.back()[1];
		return i2 * y_b_prime + j2 * y_b - g2;
	};

	// Метод Ньютона для поиска корня
	double s_opt = initial_guess;
	double tol = 1e-8;  // уменьшенная погрешность
	int max_iter = 100;

	for (int iter = 0; iter < max_iter; ++iter) {
		double F_s = boundary_conditions(s_opt);
		double delta = 1e-6;
		double F_s_prime = (boundary_conditions(s_opt + delta) - boundary_conditions(s_opt - delta)) / (2 * delta);

		double s_new = s_opt - F_s / F_s_prime;
		if (std::abs(s_new - s_opt) < tol) {
			s_opt = s_new;
			break;
		}
		s_opt = s_new;
	}

	// Проверка начальных условий и корректировка s_opt, если необходимо
	double y_a = (g1 - s_opt*i1) / j1;  // Начальное значение y(a)
	auto sol = rungeKutta(a, b, { y_a, s_opt }, p, q, f, h);

	return sol;
}

std::vector<std::vector<double>> rungekutta(double a, double b, const std::vector<double>& y0,
	std::function<double(double)> p,
	std::function<double(double)> q,
	std::function<double(double)> f,
	double h = 0.01) {
	int n = static_cast<int>((b - a) / h);
	std::vector<std::vector<double>> result(n + 1, std::vector<double>(2));
	double x = a;
	result[0] = y0;

	for (int i = 0; i < n; ++i) {
		std::vector<double> k1 = system(x, result[i], p, q, f);
		std::vector<double> y_temp = { result[i][0] + h/2 * k1[0]/2, result[i][1] + h/2 * k1[1]/2 };
		std::vector<double> k2 = system(x + h, y_temp, p, q, f);

		result[i + 1][0] = result[i][0] + h * (k1[0] + k2[0]);
		result[i + 1][1] = result[i][1] + h * (k1[1] + k2[1]);
		x += h;
	}

	return result;
}

int main() {
	// Пример функции
	auto p = [](double x) { return 1.0; };
	auto q = [](double x) { return 1.0; };
	auto f = [](double x) { return 2+2*x+x*x; };

	// Граничные условия
	double a = 0.0;
	double b = 4.0;
	double i1 = 1.0, j1 = 2.0, g1 = 0.0;
	double i2 = 3.0, j2 = 4.0, g2 = 88.0;

	double initial_guess = g1/(i1*j1)*i1;
	double h = 0.1;

	// Решение задачи
	auto solution = shootingMethod(a, b, i1, j1, g1, i2, j2, g2, p, q, f, initial_guess,h);

	// Вывод результата
	std::cout << std::fixed << std::setprecision(6);
	for (const auto& vec : solution) {
		std::cout<<"x = "<<a << ", y(x) = " << vec[0] << ", y'(x) = " << vec[1] << std::endl;
		a += h;
	}
	system("pause");
	return 0;
}