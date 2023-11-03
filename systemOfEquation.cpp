#include "systemOfEquation.h"
#include <vector>
#include <iostream>
std::vector<double> find_x(std::vector<std::vector <double>> U, std::vector <double> y)
{
	//нахождение вектора x в QR методе
	int n = y.size();
	std::vector<double> x(n, 0);
	for (int i = n - 1; i > -1; --i)
	{
		x[i] = y[i];
		for (int k = i + 1; k < n; ++k)
		{
			x[i] -= U[i][k] * x[k];

		}
		x[i] /= U[i][i];
	}
	return x;
}

double norma_vec(std::vector<double> vec)
{
	//вычисление нормы вектора 
	double norm = 0;
	int n = vec.size();
	for (int i = 0; i != n; ++i) {
		norm += pow(vec[i], 2);
	}
	if (norm == 0) return 0;
	norm = sqrt(norm);
	return norm;
}

std::vector<double> sum(std::vector<double> a, std::vector<double> b)
{
	//сложение двух векторов
	int n = a.size();
	for (int i = 0; i != n; ++i)
		a[i] += b[i];
	return a;
}

std::vector<double> mult_num(double a, std::vector<double>& z)
{
	//умножение вектора на число 
	int n = z.size();
	for (int i = 0; i < n; ++i)
		z[i] *= a;
	return z;
}

std::vector<std::vector<double>> mult_matr(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
	//умножение двух матриц
	int n = A.size();
	std::vector<double> v(n, 0);
	std::vector<std::vector<double>> R(n, v);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				R[i][j] += A[i][k] * B[k][j];
	return R;
}

std::vector<double> mult(std::vector<std::vector<double>> a, std::vector<double> b)
{
	//умножение матрицы на вектор
	int n = b.size();
	std::vector<double> vec(n, 0);
	for (int i = 0; i != n; ++i) {
		for (int j = 0; j != n; ++j) {
			vec[i] += a[i][j] * b[j];
		}
	}
	return vec;
}

std::vector<double> QR(std::vector<std::vector<double>> R, std::vector<double> b)
{
	//решение СЛАУ методом QR
	int n = b.size();
	double a;
	std::vector<double> x(n, 0), y, v, z;
	std::vector<std::vector<double>> Q(n, x), Qn, Rn, temp;
	for (int i = 0; i < n; ++i)
		Q[i][i] = 1;
	Qn = Q;
	for (int k = 0; k < n - 1; ++k) {
		y.clear();
		for (int i = k; i < n; ++i)
			y.push_back(R[i][k]);
		a = norma_vec(y);
		z.clear();
		z.push_back(1);
		for (int i = 1; i < n - k; ++i)
			z.push_back(0);
		std::vector<double> w;
		mult_num(-a, z);
		w = sum(y, z);
		w = mult_num(1 / norma_vec(w), w);
		Rn.clear();
		temp.clear();
		for (int i = 0; i < n - k; ++i)
		{
			Rn.push_back(v);
			temp.push_back(v);
		}
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				if (i >= k && j >= k)
				{
					Q[i][j] = (-2 * w[i - k] * w[j - k]);
					if (i == j) Q[i][i] += 1;
					temp[i - k].push_back(Q[i][j]);
					Rn[i - k].push_back(R[i][j]);
				}
				else
				{
					if (i == j) Q[i][j] = 1;
					else Q[i][j] = 0;
				}
			}
		}
		Qn = mult_matr(Qn, Q);
		Rn = mult_matr(temp, Rn);
		for (int i = k; i < n; ++i)
		{
			for (int j = k; j < n; ++j)
			{
				R[i][j] = Rn[i - k][j - k];
			}
		}

	}
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			Q[i][j] = Qn[j][i];
	y = mult(Q, b);
	x = find_x(R, y);
	return x;
}

std::vector<double> F(std::vector<double> x)
{
	//функция F(x, y)
	return { cos(x[1] + 0.5) + x[0] - 0.8,
			 sin(x[0]) - 2 * x[1] - 1.6 };
}

std::vector<std::vector<double>> dF(std::vector<double> x)
{
	//матрица Якоби для F(x, y)
	return { { 1, -sin(x[1] + 0.5)},
			 { cos(x[0]), -2 } };
}

std::vector<double> Fi(std::vector<double> x, double h)
{
	//функция Ф(h, x, y)
	return { (h / 10) * cos(x[1] + 0.5) + x[0] - 0.8, (h / 10) * sin(x[0]) - 2 * x[1] - 1.6 };
}
std::vector<std::vector<double>> dFi(std::vector<double> x, double h)
{
	//матрица Якоби для Ф(h, x, y)
	return { { 1, -(h / 10) * sin(x[1] + 0.5)},
			 { (h / 10) * cos(x[0]), -2 } };
}

std::vector<double> initial_approximation()
{
	//нахождение начального приближения используя замечание 1.4
	std::vector<double> x = { 0.8, -0.8 }, fi;
	for (int i = 1; i < 10; ++i) {
		fi = Fi(x, i);
		x = sum(QR(dFi(x, i), mult_num(-1, fi)), x);
	}
	return x;
}

std::vector<double> system_equations(std::vector<double> x)
{
	//метод Ньютона для системы уравнений
	std::vector<double> tmp, d, f;
	int count = 0;
	do {
		++count;
		tmp = x;
		f = F(x);
		d = QR(dF(x), mult_num(-1, f)); // вычисление delta X(k)
		x = sum(x, d); // вычисление X(k+1)
	} while (norma_vec(sum(x, mult_num(-1, tmp))) > 0.0001);
	std::cout << "Решение: ";
	for (auto& c : x)
		std::cout << c << " ";
	std::cout << std::endl << "Количество итераций: " << count << std::endl;
	return x;
}
