#ifndef BACK_TRACKING_SEARCH_H
#define BACK_TRACKING_SEARCH_H
#define EPS 1e-9
#include <iostream>
#include <vector>
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace std;
class BackTrackingSearch
{
	//shrink rate
	double rate;
	//initial search point
	vector<double> searchPoint;
	//c value
	//0 < c1 < c2 < 1
	double c1;
	double c2;
	//alpha value
	//initial alpha
	double maxAlpha;
	double minAlpha;
	double lo;
	double hi;
	//store the former values

	//function pointer to origin function
	double (*func)(vector<double> &);
	//derivative of origin function
	vector<double> (*deriv)(vector<double> &);
	//search direction function
	vector<double> (*direc)(vector<double> &);
	public:
	BackTrackingSearch();
	//constructor
	BackTrackingSearch(double (*_func)(vector<double> &), vector<double> (*_deriv)(vector<double> &), vector<double> (* _direc)(vector<double> &), vector<double>& _searchPoint);	

	void setParameter(double _rate, double _c1, double _c2, double _maxAlpha, double _minAlpha);
	//search method
	vector<double> search();
	//search range
	bool searchRange(vector<double> pk, double & alpha);
	//accept the step length
	double zoomIn(vector<double> pk);
	//destructor
	~BackTrackingSearch();
};
#endif
