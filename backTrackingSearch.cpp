#include "backTrackingSearch.h"
BackTrackingSearch::BackTrackingSearch()
{
	//empty constructor
}
//people really care about function, derivative function, and initial points
BackTrackingSearch::BackTrackingSearch(double (*_func)(vector<double> &), vector<double> (*_deriv)(vector<double> &), vector<double> (* _direc)(vector<double> &), vector<double>& _searchPoint)
{
		func = _func;
		deriv = _deriv;
		direc = _direc;
		searchPoint = _searchPoint;
		//set up some init parameters
		rate = 2;
		c1 = 1e-2;
		c2 = 0.1;
		maxAlpha = 1000;
		minAlpha = 1e-3;
}

void BackTrackingSearch::setParameter(double _rate, double _c1, double _c2, double _maxAlpha, double _minAlpha)
{
		rate = _rate;
		c1 = _c1;
		c2 = _c2;
		maxAlpha = _maxAlpha;
		minAlpha = _minAlpha;
}
	//search function returns a local minimizer
vector<double> BackTrackingSearch::search()
{
		int step = 0;
		while (true) //the while loop for searching the min
		{
			lo = 0;
			hi = minAlpha;
			double alpha = 0;
			vector<double> pk = direc(searchPoint);//the search direction
			while (true)
			{
				if (searchRange(pk, alpha))
					break;
			}//keep searching until return true
			cout << step << ": searchRange " << hi << endl;
			if (fabs(alpha) < EPS)
				alpha = zoomIn(pk);
			cout << step++ << ": zoomIn " << alpha << endl;
			double originValue = func(searchPoint);
			searchPoint[0] += pk[0] * alpha;
			searchPoint[1] += pk[1] * alpha;
			cout << step << ": search points are: " << searchPoint[0] << "  " << searchPoint[1] << endl;
			if (fabs(func(searchPoint) - originValue) < EPS)
				return searchPoint;
		}
}
bool BackTrackingSearch::searchRange(vector<double> pk, double & alpha)
{
	//base case, terminate search direcly
	if (hi >= maxAlpha)
	{
		alpha = hi;
		return true;
	}
	//get the sufficient value
	double originValue = func(searchPoint);
	vector<double> derivative = deriv(searchPoint);
	int size = searchPoint.size();
	for (int i = 0; i < size; i++)
		originValue += derivative[i] * pk[i] * c1 * hi;
	vector<double> loPoint = searchPoint;
	vector<double> hiPoint = searchPoint;
	for (int i = 0; i < size; i++)
	{
		loPoint[i] += pk[i] * lo;
		hiPoint[i] += pk[i] * hi;
	}
	if (func(hiPoint) >= originValue || func(hiPoint) >= func(loPoint))
		return true; //zoom in
	double diZero = derivative[0] * pk[0] + derivative[1] * pk[1];
	vector<double> deriHi = deriv(hiPoint);
	double diHi = deriHi[0] * pk[0] + deriHi[1] * pk[1];
	if (fabs(diHi) <= fabs(diZero) * c2) //this is the point that we want
	{
		alpha = hi;
		return true; 
	}
	if (diHi > 0)
		return true; //zoom in
	lo = hi;
	hi = hi * rate;
	return false; //keep searching
}
double BackTrackingSearch::zoomIn(vector<double> pk)
{
	vector<double> derivative = deriv(searchPoint);
	while (true)
	{
		double mid = (lo + hi) / 2;
		vector<double> midPoint = searchPoint;
		vector<double> lowPoint = searchPoint;
		midPoint[0] += mid * pk[0];
		midPoint[1] += mid * pk[1];
		lowPoint[0] += lo * pk[0];
		lowPoint[1] += lo * pk[1];
		double originValue = func(searchPoint);
		originValue += derivative[0] * pk[0] * c1 * mid;
		originValue += derivative[1] * pk[1] * c1 * mid;
		if (func(midPoint) >= func(lowPoint) || func(midPoint) >= originValue)
		{
			hi = mid;
			continue;
		}
		vector<double> derivMid = deriv(midPoint);
		double deMid = derivMid[0] * pk[0] + derivMid[1] * pk[1];
		double deZero = derivative[0] * pk[0] + derivative[1] * pk[1];
		if (fabs(deMid) <= fabs(deZero) * c2)
			return mid;
		if (derivMid[0] * pk[0] * (hi - lo) + derivMid[1] * pk[1] * (hi - lo) >= 0)
			hi = lo;
		lo = mid;
	}
}
BackTrackingSearch::~BackTrackingSearch()
{
	//nothing to do here
}
