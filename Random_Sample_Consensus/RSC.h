#ifndef RSC_H
#define RSC_H

#include <iostream>
#include <vector>
#include<fstream>
#include "Site_type.h"
using namespace std;

class RSC
{
public:
	double Total_Distance(vector<Site>Site_Array , vector<double>Function_Parameter);
	pair<vector<double>, double> Fittest_Parameter(int Times_K, vector<Site>Site_Array, double Average);
	int Count_inliers(vector<Site> Site_Array, double Average, vector<double> Function_Parameter);
	vector<double> Random_Parameter(vector<Site>Site_Array);
private:

};


#endif // !RSC_H