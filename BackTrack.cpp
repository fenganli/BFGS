#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include "backTrackingSearch.h"

using namespace boost::numeric::ublas;
using namespace std;

matrix<double> inverseMat(matrix<double>* inputMat);
double rosen(std::vector<double> & point);
std::vector<double> rosenDeriv(std::vector<double> & point);
std::vector<double> rosenPk (std::vector<double> & point);
std::vector<double> rosenPk2 (std::vector<double> & point);
std::vector<double> rosenPk3 (std::vector<double> & point);
double get_element(matrix<double> const &m,int i,int j);
int main(int argc, char * argv[])
{
	std::vector<double> inputVec;
	inputVec.push_back(atof(argv[2]));
	inputVec.push_back(atof(argv[3]));
	BackTrackingSearch * searchLength;
	if (atof(argv[1]) == 1)
		searchLength = new BackTrackingSearch(rosen, rosenDeriv, rosenPk, inputVec);
	else if (atof(argv[1]) == 2)
		searchLength = new BackTrackingSearch(rosen, rosenDeriv, rosenPk2, inputVec);
	else if (atof(argv[1]) == 3)
		searchLength = new BackTrackingSearch(rosen, rosenDeriv, rosenPk3, inputVec);
	std::vector<double> minimizer = searchLength -> search();
	cout << minimizer[0] << "      " << minimizer[1] << endl;
	return 0;
}
double rosen(std::vector<double> & point)
{
	return 100 * pow(point[1] - point[0] * point[0], 2) + pow(1 - point[0], 2);
}
std::vector<double> rosenDeriv(std::vector<double> & point)
{
	std::vector<double> returnVec;
	returnVec.push_back(-400 * point[0] * point[1] + 400 * pow(point[0], 3) + 2 * point[0] -2);
	returnVec.push_back(200 * (point[1] - point[0] * point[0]));
	return returnVec;
}
std::vector<double> rosenPk(std::vector<double> & point)
{
	std::vector<double> vec = rosenDeriv(point);
	vec[0] = 0 - vec[0];
	vec[1] = 0 - vec[1];
	return vec;
}
std::vector<double> rosenPk2(std::vector<double> &point)
{
	double a11 = -400 * point[1] + 1200 * point[0] * point[0] + 2;
	double a12 = -400 * point[0];
	double a21 = -400 * point[0];
	double a22 = 200;
	std::vector<double> deriv = rosenDeriv(point);
	std::vector<double> returnVec;
	double first = -(deriv[0] * a22 + deriv[1] * (-a12)) / (a11 * a22 - a12 * a21);
	double second = -(deriv[0] * (-a21) + deriv[1] * a11) / (a11 * a22 - a12 * a21);
	returnVec.push_back(first);
	returnVec.push_back(second);
	return returnVec;
}
std::vector<double> rosenPk3 (std::vector<double> & point)
{
	using namespace boost::numeric::ublas;
	static matrix<double> * old_x = NULL;
	static matrix<double> * old_diff = NULL;
	static matrix<double> * old_heissan = NULL;
	if (old_x == NULL)
	{
		int size = point.size();
		//initialize these old values
		old_x = new matrix<double> (1, size);
		for (int i = 0; i < size; i++)
			(*old_x)(0, i) = point[i];
		old_diff = new matrix<double> (1, size);
		(*old_diff)(0, 0) = (-400 * point[0] * point[1] + 400 * pow(point[0], 3) + 2 * point[0] - 2);
		(*old_diff)(0, 1) = (200 * (point[1] - point[0] * point[0]));
		old_heissan = new matrix<double> (size, size);
		(*old_heissan)(0, 0) = -400 * point[1] + 1200 * point[0] * point[0] + 2;
		(*old_heissan)(0, 1) = -400 * point[0];
		(*old_heissan)(1, 0) = -400 * point[0];
		(*old_heissan)(1, 1) = 200;
		matrix<double> invMat = inverseMat(old_heissan);
		matrix<double> pk = prod(invMat, trans(*old_diff));
		std::vector<double> returnVec;
		returnVec.push_back(0 - pk(0, 0));
		returnVec.push_back(0 - pk(1, 0));
		return returnVec;
	}
	else
	{
		int size = point.size();
		matrix<double> sk(1, size);
		matrix<double> yk(1, size);
		for (int i = 0; i < size; i++)
			sk(0, i) = point[i] - (*old_x)(0, i);
		matrix<double> new_diff(1, size);
		new_diff(0, 0) = (-400 * point[0] * point[1] + 400 * pow(point[0], 3) + 2 * point[0] - 2);
		new_diff(0, 1) = (200 * (point[1] - point[0] * point[0]));
		for (int i = 0; i < size; i++)
			yk(0, i) = new_diff(0, i) - (*old_diff)(0, i);
		std::cout << sk(0,0) * yk(0,0) + sk(0, 1) * yk(0, 1) <<std::endl;
		matrix<double> new_heissan(size, size);
		matrix<double> matrix1 = (*old_heissan);
		matrix1 = prod(matrix1, trans(sk));
		matrix1 = prod(matrix1, sk);
		matrix1 = prod(matrix1, *old_heissan);
		matrix<double> matrix_1 = sk;
		matrix_1 = prod(matrix_1, *old_heissan);
		matrix_1 = prod(matrix_1, trans(sk));
		matrix1 /= matrix_1(0, 0);
		matrix<double> matrix2 = prod(trans(yk), yk);
		matrix2 /= (prod(yk, trans(sk)))(0,0);
		new_heissan = (*old_heissan) - matrix1 + matrix2;
		matrix<double> invMat = inverseMat(&new_heissan);
		matrix<double> pk = prod(invMat, trans(new_diff));
		std::vector<double> returnVec;
		returnVec.push_back(0 - pk(0, 0));
		returnVec.push_back(0 - pk(1, 0));
		(*old_x)(0, 0) = point[0];
		(*old_x)(0, 1) = point[1];
		(*old_diff)(0, 0) = new_diff(0,0);
		(*old_diff)(0, 1) = new_diff(0,1);
		(*old_heissan)(0, 0) = new_heissan(0,0);
		(*old_heissan)(0, 1) = new_heissan(0,1);
		(*old_heissan)(1, 0) = new_heissan(1,0);
		(*old_heissan)(1, 1) = new_heissan(1,1);
		return returnVec;
	}
}
matrix<double> inverseMat(matrix<double>* inputMat)
{
	using namespace boost::numeric::ublas;
	double denomitor =(*inputMat)(0, 0) * (*inputMat)(1, 1) -(*inputMat)(0, 1) * (*inputMat)(1, 0);
	matrix<double> returnMat(2,2);
	returnMat(0, 0) = (*inputMat)(1, 1) / denomitor;
	returnMat(0, 1) = - (*inputMat)(0, 1) / denomitor;
	returnMat(1, 0) = - (*inputMat)(1, 0) / denomitor;
	returnMat(1, 1) = (*inputMat)(0, 0) / denomitor;
	return returnMat;
}
double get_element(matrix<double> const &m,int i,int j) {
    return m(i,j);
}
