#ifndef SYSTEMOFEQUATION_H
#define SYSTEMOFEQUATION_H
#include <vector>
std::vector<double> find_x(std::vector<std::vector <double>> U, std::vector <double> y);

double norma_vec(std::vector<double> vec);

std::vector<double> sum(std::vector<double> a, std::vector<double> b);

std::vector<double> mult_num(double a, std::vector<double>& z);

std::vector<std::vector<double>> mult_matr(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);

std::vector<double> mult(std::vector<std::vector<double>> a, std::vector<double> b);

std::vector<double> QR(std::vector<std::vector<double>> R, std::vector<double> b);

std::vector<double> F(std::vector<double> x);

std::vector<std::vector<double>> dF(std::vector<double> x);

std::vector<double> Fi(std::vector<double> x, double h);

std::vector<std::vector<double>> dFi(std::vector<double> x, double h);

std::vector<double> initial_approximation();

std::vector<double> system_equations(std::vector<double> x);


#endif