#include "equation.h"
#include <iostream>

double f(double x)
{
	//������� f
	return x + log10(x) - 0.5;
}

double fd(double x)
{
	//����������� ������� f
	return (log(10) * x + 1) / (log(10) * x);
}

std::pair<double, double> localization()
{
	//����������� �����
	double b = 100, a = 0;
	int N = 10;
	double h = (b - a) / N;
	for (int k = 0; k < N; ++k)
		if (f(a + k * h) * f(a + (k + 1) * h) < 0)
			return std::make_pair(a + k * h, a + (k + 1) * h);
}

double newton()
{
	//����� ������� ��� ������� ����������� ���������
	auto interval = localization();
	double x = interval.second, tmp;
	int count = 0;

	do {
		++count;
		tmp = x;
		x = x - f(x) / fd(x);

		if (x > interval.second || x < interval.first)
		{ //�������� �� "�����" �����������  �� ���������
			x = (interval.first + interval.second) / 2;
			std::cout << "����� ����������� ������� �������� �� " << count << " ����" << std::endl;
		}

		if (f(x) < 0) interval = std::make_pair(x, interval.second); //�������������� ���������
		else interval = std::make_pair(interval.first, x);
	} while (abs(tmp - x) > 0.0001);

	std::cout << "�������: " << x << std::endl;
	std::cout << "���������� ��������: " << count << std::endl;
	return x;
}
