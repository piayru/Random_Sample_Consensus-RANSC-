#include <iostream>
#include <vector>

using namespace std;

class Site
{
public:
	void set_x(double x) { this->x = x; };
	void set_y(double y) { this->y = y; };
	double get_x() { return this->x; };
	double get_y() { return this->y; };
private:
	double x;
	double y;
};

class RANSAC
{
public:
	vector<Site> Set_Site_Array(double site_array[]);
	vector<vector<double>> Equation_Array(vector<Site> Site_Array, int Parameter_Number);
	vector<vector<double>> Traspose_Array(vector<vector<double>> Target_Array);
	vector<vector<double>> Dot_Array(vector<vector<double>> Array_One, vector<vector<double>> Array_Two);
	vector<vector<double>> Inverse_Array(vector<vector<double>> Target_Array);


	vector<vector<double>> minor(vector<vector<double>> a, int i, int n);
	float det(vector<vector<double>> a, int n);
	vector<vector<double>> RANSAC::transpose(vector<vector<double>> c, int n, float det);
	vector<vector<double>> RANSAC::cofactor(vector<vector<double>> a, int n, float determinte);
	vector<vector<double>> RANSAC::inverse(vector<vector<double>> a, int n, float det);
private:

};

