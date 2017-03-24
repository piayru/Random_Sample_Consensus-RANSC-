#ifndef FIT_EQUATION_H
#define FIT_EQUATION_H

#include "Site_type.h"
#include <iostream>
#include <vector>
#include<fstream>
using namespace std;


class Fit_Equation
{
public:
	vector<double> Handle(double site_array[] , int Parameter_Number);
	vector<Site> Set_Site_Array(double site_array[]);
	vector<vector<double>> Equation_Array(vector<Site> Site_Array, int Parameter_Number);
	vector<vector<double>> Traspose_Array(vector<vector<double>> Target_Array);
	vector<vector<double>> Dot_Array(vector<vector<double>> Array_One, vector<vector<double>> Array_Two);
	vector<vector<double>> Inverse_Array(vector<vector<double>> Target_Array);
	vector<vector<double>> Minor(vector<vector<double>> Target_Array, int Index_Matrix, int Matrix_Number);
	float Det(vector<vector<double>> Target_Array, int Matrix_Number);
	vector<vector<double>> Transpose(vector<vector<double>> Target_Array, int Matrix_Number, float Det);
	vector<vector<double>> Cofactor(vector<vector<double>> Target_Array, int Matrix_Number, float Determinte);
	vector<vector<double>> Inverse(vector<vector<double>> Target_Array, int Matrix_Number, float Det);
	vector<double> Handle(vector<Site> site_array, int Parameter_Number);
private:

};



#endif // !FIT_EQUATION_H