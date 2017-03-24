#include "RSC.h"
#include "Fit_Equation.h"
#include <math.h>
#include <time.h>

#include <map>
#include<fstream>
using namespace std;

#define Site_Number 5 
int const Parameter_Number = 3;


double RSC::Total_Distance(vector<Site>Site_Array, vector<double>Function_Parameter) {
	double Y_p = 0.0;
	long double Sum_Distance = 0.0;
	for (int Index_Site = 110; Index_Site < 290/*Site_Array.size()*/; Index_Site++) {
		for (int Index_Pow_X = Function_Parameter.size() - 1; Index_Pow_X >= 0; Index_Pow_X--)
			Y_p += Function_Parameter[Function_Parameter.size() - Index_Pow_X - 1] * pow(Site_Array[Index_Site].get_x(), Index_Pow_X);
		if (Y_p > 0)
			Sum_Distance += pow((Y_p - Site_Array[Index_Site].get_y()), 2);
		Y_p = 0.0;
	}
	return Sum_Distance;
}

pair<vector<double>, double> RSC::Fittest_Parameter(int Times_K, vector<Site> Site_Array, double Average) {
	vector<double> First_Parameter = Random_Parameter(Site_Array);
	int First_Distance = Count_inliers(Site_Array, Average, First_Parameter);
	pair<vector<double>, int> Fittest(First_Parameter, First_Distance);
	int MAX_Count = First_Distance;

	for (int Index_Times = 1; Index_Times < Times_K; Index_Times++) {
		vector<double> One_Parameter = Random_Parameter(Site_Array);
		int Inliers_Number = Count_inliers(Site_Array, Average, One_Parameter);
		if (MAX_Count < Inliers_Number) {
			Fittest.first = One_Parameter;
			Fittest.second = Inliers_Number;
		}
	}
	return Fittest;
}


vector<double> RSC::Random_Parameter(vector<Site>Site_Array) {
	vector<Site> Fitting_Site;
	srand((unsigned)time(NULL));
	map<int, double> Random_number;
	do {
		int Ran = rand() % 40;
		if (Random_number.find(Ran) == Random_number.end())
			Random_number.insert(pair<int, double>(Ran, (double)Ran));
	} while (Random_number.size() <= Site_Number - 1);
	for (auto i = Random_number.begin(); i != Random_number.end(); ++i) 
		Fitting_Site.push_back(Site_Array[130 + (i->first)]);
	Fit_Equation RANSAC;
	return RANSAC.Handle(Fitting_Site, Parameter_Number);
}

int RSC::Count_inliers(vector<Site> Site_Array, double Average, vector<double> Function_Parameter) {
	int Count = 0;
	double Y_p = 0.0;
	long double Sum_Distance = 0.0;
	for (int Index_Site = 110; Index_Site < 290/*Site_Array.size()*/; Index_Site++) {
		for (int Index_Pow_X = Function_Parameter.size() - 1; Index_Pow_X >= 0; Index_Pow_X--)
			Y_p += Function_Parameter[Function_Parameter.size() - Index_Pow_X - 1] * pow(Site_Array[Index_Site].get_x(), Index_Pow_X);
		if (Y_p > 0 && Y_p <= Average)
			Count++;
	}
	return Count;
}