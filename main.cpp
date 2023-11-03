#include "equation.h"
#include "systemOfEquation.h"
#include <iostream>

int main() {
	setlocale(LC_ALL, "Russian");
	std::cout << "������� ����������� ��������� ������� �������:" << std::endl;
	newton();
	std::cout << std::endl << std::endl << "������� ������� ��������� ������� �������:" << std::endl << std::endl;
	std::cout << "��������� ����������� ��������� ����������� �������: " << std::endl;
	system_equations({ -0.1, -0.8 }); // ������� ����������
	std::cout << std::endl;
	std::cout << "��������� ����������� ��������� � ������� ��������� 1.4: " << std::endl;
	system_equations(initial_approximation());
	return 0;
}