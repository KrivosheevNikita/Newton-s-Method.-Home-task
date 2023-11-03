#include "equation.h"
#include "systemOfEquation.h"
#include <iostream>

int main() {
	setlocale(LC_ALL, "Russian");
	std::cout << "Решение нелинейного уравнения методом Ньютона:" << std::endl;
	newton();
	std::cout << std::endl << std::endl << "Решение системы уравнений методом Ньютона:" << std::endl << std::endl;
	std::cout << "Начальное приближение находится графическим методом: " << std::endl;
	system_equations({ -0.1, -0.8 }); // найдено графически
	std::cout << std::endl;
	std::cout << "Начальное приближение находится с помощью замечания 1.4: " << std::endl;
	system_equations(initial_approximation());
	return 0;
}