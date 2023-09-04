#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>


#include <bits/stdc++.h>


#define rhoc 1000
#define sigma 0.07
#define c2 2.047

#define epsilon 1
#define dk 0.006

#define M_MC 500
#define N_MC 500

#define M_NC 50
#define N_NC 50

#define M_GL 5
#define N_GL 5



double func(double xi, double fbv)
{
	double cf = std::pow(fbv, 2.0/3.0) + std::pow(1.0 - fbv, 2.0/3.0) - 1.0;
	double chi = 12.0 * cf * sigma / c2 / rhoc / std::pow(epsilon, 2.0/3.0) / std::pow(dk, 5.0/3.0) / std::pow(xi, 11.0/3.0);
	return std::pow(1.0 + xi, 2.0) * std::pow(xi, -11.0/3.0) * std::exp(-chi);

}

std::vector<double> randomGenerator(const double min, const double max, const int n) // min, max > 0
{
	std::vector<double> arr(n);
	for (int i = 0; i < n; i++)
	{
		arr[i] = min + static_cast<double> (rand()) / (static_cast<double> (RAND_MAX / (max - min)));
	}
	return arr;
}

std::vector<double> gausslegenre(const std::vector<double>& f)
{
	std::vector<double> A1(M_GL,0.0);
	std::vector<double> A2(N_GL,0.0);

	std::vector<double> P1(M_GL,0.0);
	std::vector<double> P2(N_GL,0.0);
	
//-------------------------------------------------------------------------//

	A1[0] = A2[0] = 0.2369268851;
	A1[1] = A2[1] = 0.2369268851;
	A1[2] = A2[2] = 0.4786286705;
	A1[3] = A2[3] = 0.4786286705;
	A1[4] = A2[4] = 0.5688888889;

	P1[0] = P2[0] = -0.9061798459;
	P1[1] = P2[1] = 0.9061798459;
	P1[2] = P2[2] = -0.5384693101;
	P1[3] = P2[3] = 0.5384693101;
	P1[4] = P2[4] = 0.0000000000;

//-------------------------------------------------------------------------//

	double a = 0.0;
	double b = 1.0;
	double c = 0.0;
	double d = 1.0;
	int N_f = f.size();

	std::vector<double> numerator(N_f,0.0);
	std::vector<double> result(N_f, 0.0);

	for (int i = 0; i < N_f; i++)
	{
		double Fbv = f[i];
		for (int j = 0; j < M_GL; j++)
		{
			double Xi = (a + b) / 2.0 + (b - a) * P1[j] / 2.0;
			numerator[i] += A1[j] * func(Xi,Fbv);
		}
		numerator[i] = numerator[i] * (b - a) / 2.0;
	}

	double denominator = 0.;
	for (int i = 0; i < M_GL; i++)
	{
		double Xi = (a + b) / 2.0 + (b - a) * P1[i] / 2.0;
		for (int j = 0; j < N_GL; j++)
		{
			double Fj = (c + d) / 2.0 + (d - c) * P2[j] / 2.0;
			denominator += A1[i] * A2[j] * func(Xi,Fj);
		}
	}
	denominator = denominator * (b - a) * (d - c) / 4.0;

	for (int k = 0; k < N_f; k++)
	{
		result[k] = 2 * numerator[k] / denominator;
	}
	return result;
}

std::vector<double> newtoncotes(const std::vector<double>& f)
{
	double a = 0.001;
	double b = 1.0 - 0.001;
	double c = 0.001;
	double d = 1.0 - 0.001;
	int N_f = f.size();
	
	double h = (b - a) / M_NC;
	double t = (d - c) / N_NC;

	std::vector<double> X(M_NC + 1, 0.0);
	std::vector<double> Y(N_NC + 1, 0.0);

	for (int k = 0; k < M_NC + 1; k++)
	{
		X[k] = a + h * k;
		Y[k] = c + t * k;
	}

	std::vector<std::vector<double>> lamd(M_NC + 1);
	for (auto& eachRow : lamd)
	{
		eachRow.resize(N_NC);
	}

	for (int i = 0; i < M_NC + 1; i++)
	{
		for (int j = 0; j < N_NC + 1; j++)
		{
			if ((i == 0 && j == 0) || (i == 0 && j == N_NC) || (i == M_NC && j == 0) || (i == M_NC && j == N_NC))
				lamd[i][j] = 1.0;
			else if ((i >= 1 && i <= M_NC - 1) && (j == 0 || j == N_NC))
				lamd[i][j] = 2.0;
			else if ((j >= 1 && j <= N_NC - 1) && (i == 0 || i == M_NC))
				lamd[i][j] = 2.0;
			else
				lamd[i][j] = 4.0;
		}
	}

	std::vector<double> numerator(N_f,0.0);
	std::vector<double> result(N_f, 0.0);

	for (int i = 0; i < N_f; i++)
	{
		double Fbv = f[i];
		numerator[i] = func(a, Fbv) + func(b, Fbv);
		for (int j = 0; j < M_NC + 1; j++)
		{
			double Xi = X[j];
			numerator[i] += 2 * func(Xi,Fbv);
		}
		numerator[i] = numerator[i] * h / 2.0;
	}

	double denominator = 0.;
	for (int i = 0; i < M_NC + 1; i++)
	{
		double Xi = X[i];
		for (int j = 0; j < N_NC + 1; j++)
		{
			double Fj = Y[j];
			denominator += lamd[i][j] * func(Xi,Fj);
		}
	}
	denominator = denominator * h * t / 4.0;

	for (int k = 0; k < N_f; k++)
	{
		result[k] = 2 * numerator[k] / denominator;
	}
	return result;
}


std::vector<double> montecarlo(const std::vector<double>& f)
{
	srand(time(0));

	double a = 0.0;
	double b = 1.0;
	double c = 0.0;
	double d = 1.0;
	int N_f = f.size();

	std::vector<double> X = randomGenerator(a,b,M_MC);
	std::vector<double> Y = randomGenerator(c,d,N_MC);

	std::vector<double> numerator(N_f,0.0);
	std::vector<double> result(N_f, 0.0);

	for (int i = 0; i < N_f; i++)
	{
		double Fbv = f[i];
		for (int j = 0; j < M_MC; j++)
		{
			double Xi = X[j];
			numerator[i] += func(Xi, Fbv);
		}

		numerator[i] = numerator[i] * (b - a) / M_MC;
	}


	
	double denominator = 0.;
	for (int i = 0; i < M_MC; i++)
	{
		for (int j = 0; j < N_MC; j++)
		{
			double Xi = X[i];
			double Fj = Y[j];
			denominator += func(Xi,Fj);
		}
	}
	denominator = denominator * (b - a) * (d - c) / M_MC / N_MC;
	
	for (int k = 0; k < N_f; k++)
	{
		result[k] = 2 * numerator[k] / denominator;
	}
	return result;
}


int main()
{

	double f_min = 0.1;
	double f_max = 0.9;
	double delta_f = 0.05;

	int N_f = (f_max - f_min) / delta_f + 1;

	std::vector<double> fbv(N_f, 0.0);
	for (int k = 0; k < N_f; k++)
	{
		fbv[k] = f_min + k * delta_f;
	}

//	std::vector<double> x_random = randomGenerator(0.0, 1.0, x_num);

	std::clock_t start_mc, end_mc;
	std::clock_t start_gl, end_gl;
	std::clock_t start_nc, end_nc;
//------------------------------------------------------------//
	start_mc = clock();

	std::vector<double> dsd_mc = montecarlo(fbv);

	std::cout << "Results for Monte Carlo:" << std::endl;

	for (int k = 0; k < N_f; k++)
	{
		std::cout << dsd_mc[k] << std::endl;
	}
	
	end_mc = clock();
	double time_mc = double(end_mc - start_mc) / double(CLOCKS_PER_SEC);
	std::cout << "Time elapsed by monte carlo is: " << std::fixed << time_mc << std::setprecision(5);
	std::cout << " sec " << std::endl;
//------------------------------------------------------------//
	start_gl = clock();

	std::vector<double> dsd_gl = gausslegenre(fbv);

	std::cout << "Results for Gauss Legendre:" << std::endl;

	for (int k = 0; k < N_f; k++)
	{
		std::cout << dsd_gl[k] << std::endl;
	}

	end_gl = clock();
	double time_gl = double(end_gl - start_gl) / double(CLOCKS_PER_SEC);
	std::cout << "Time elapsed by gauss legendre is: " << std::fixed << time_gl << std::setprecision(5);
	std::cout << " sec " << std::endl;
//------------------------------------------------------------//
	start_nc = clock();

	std::vector<double> dsd_nc = newtoncotes(fbv);

	std::cout << "Results for Newton Cotes:" << std::endl;

	for (int k = 0; k < N_f; k++)
	{
		std::cout << dsd_nc[k] << std::endl;
	}

	end_nc = clock();
	double time_nc = double(end_nc - start_nc) / double(CLOCKS_PER_SEC);
	std::cout << "Time elapsed by newton cotes is: " << std::fixed << time_nc << std::setprecision(5);
	std::cout << " sec " << std::endl;
}



