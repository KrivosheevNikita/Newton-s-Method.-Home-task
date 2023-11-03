#include "equation.h"
#include <iostream>

double f(double x)
{
	//функция f
	return x + log10(x) - 0.5;
}

double fd(double x)
{
	//производная функции f
	return (log(10) * x + 1) / (log(10) * x);
}

std::pair<double, double> localization()
{
	//локализация корня
	double b = 100, a = 0;
	int N = 10;
	double h = (b - a) / N;
	for (int k = 0; k < N; ++k)
		if (f(a + k * h) * f(a + (k + 1) * h) < 0)
			return std::make_pair(a + k * h, a + (k + 1) * h);
}

double newton()
{
	//метод ньютона для решения неленийного уравнения
	auto interval = localization();
	double x = interval.second, tmp;
	int count = 0;

	do {
		++count;
		tmp = x;
		x = x - f(x) / fd(x);

		if (x > interval.second || x < interval.first)
		{ //проверка на "вылет" приближения  из интервала
			x = (interval.first + interval.second) / 2;
			std::cout << "Метод половинного деления сработал на " << count << " шаге" << std::endl;
		}

		if (f(x) < 0) interval = std::make_pair(x, interval.second); //пересчитывание интервала
		else interval = std::make_pair(interval.first, x);
	} while (abs(tmp - x) > 0.0001);

	std::cout << "Решение: " << x << std::endl;
	std::cout << "Количество итераций: " << count << std::endl;
	return x;
}
