#include "RANSAC.h"

int const Parameter_Number = 3;
//Site *Site_Array;
double site[] = { -9 , 66 , -5 , 18 , 13 , 198 , 17 , 326 };
int main() {
	RANSAC RANSAC;
	
	vector<vector<double>> X_Array = RANSAC.Equation_Array(RANSAC.Set_Site_Array(site) , Parameter_Number);
	vector<vector<double>> Transpose_X = RANSAC.Traspose_Array(X_Array);
	vector<vector<double>> a = RANSAC.Dot_Array(Transpose_X, X_Array);
	vector<vector<double>> b = RANSAC.Inverse_Array(a);
	vector<vector<double>> q = RANSAC.Dot_Array(b, Transpose_X);
	vector<vector<double>> y;
	for (int i = 0; i < 4; i++) {
		vector<double> T{ site[2 * i + 1] };
		y.push_back(T);
	}
	vector<vector<double>> ans = RANSAC.Dot_Array(q, y);
	int i = 0;
	cout << i << endl;
	return 0;
}