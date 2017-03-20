#include "RANSAC.h"
#include <math.h>

//int array to Site array.
vector<Site> RANSAC::Set_Site_Array(double site_array[]) {
	int Array_Number = sizeof(site_array);
	vector<Site> Site_Array;
	for (int Index_Number = 0; Index_Number < Array_Number / 2; Index_Number++) {
		Site Temp_Site;
		Temp_Site.set_x(site_array[2 * Index_Number]);
		Temp_Site.set_y(site_array[2 * Index_Number + 1]);
		Site_Array.push_back(Temp_Site);
		//int a = Site_Array[Index_Number].get_x();
		//cout << Site_Array[Index_Number].get_x() << Site_Array[Index_Number].get_y() << endl;
	}
	return  Site_Array;
}


vector<vector<double>> RANSAC::Equation_Array(vector<Site> Site_Array, int Parameter_Number) {
	vector<vector<double>> Temp_Equation_Array;
	int index = 0;
	for (int Index_Site_Number = 0; Index_Site_Number < Site_Array.size(); Index_Site_Number++) {
		vector<double> Temp_row;
		for (int Index_Parameter_Number = Parameter_Number - 1; Index_Parameter_Number >= 0; Index_Parameter_Number--) {
			Temp_row.push_back(pow(Site_Array[Index_Site_Number].get_x(), Index_Parameter_Number));
			index++;
		}
		Temp_Equation_Array.push_back(Temp_row);
	}

	return Temp_Equation_Array;
}

vector<vector<double>> RANSAC::Traspose_Array(vector<vector<double>> Target_Array) {
	vector<vector<double>> Temp_array;
	for (int Index_Target_Number = 0; Index_Target_Number < Target_Array[0].size(); Index_Target_Number++) {
		vector<double> Temp_row;
		for (int Index_Target_Array = 0; Index_Target_Array < Target_Array.size(); Index_Target_Array++) {
			Temp_row.push_back(Target_Array[Index_Target_Array][Index_Target_Number]);
		}
		Temp_array.push_back(Temp_row);
	}
	return Temp_array;
}

vector<vector<double>> RANSAC::Dot_Array(vector<vector<double>> Array_One, vector<vector<double>> Array_Two) {
	vector<vector<double>> One_Dot_Two;
	for (int Index_Array_One = 0; Index_Array_One < Array_One.size(); Index_Array_One++) {
		vector<double> Temp_Dot;
		for (int Index_Array_Two = 0; Index_Array_Two < Array_Two[0].size(); Index_Array_Two++) {
			Temp_Dot.push_back(0.0);
			for (int Index_Dot = 0; Index_Dot < Array_One[Index_Array_One].size(); Index_Dot++) {
				Temp_Dot[Index_Array_Two] += Array_One[Index_Array_One][Index_Dot] * Array_Two[Index_Dot][Index_Array_Two];
			}
		}
		One_Dot_Two.push_back(Temp_Dot);
	}
	return One_Dot_Two;
}

vector<vector<double>> RANSAC::Inverse_Array(vector<vector<double>> a) {
	int i, j, n;
	float deter;
	vector<vector<double>> ans;
	n = 3;	// read function
	deter = (float)det(a, n);	// read function
	ans = inverse(a, n, deter);	// read function
	return ans;
}



//---------------------------------------------------
//	calculate minor of matrix OR build new matrix : k-had = minor
vector<vector<double>> RANSAC::minor(vector<vector<double>> a, int i, int n) {
	vector<vector<double>> b;
	int j, l, h = 0, k = 0;
	vector<double> c;
	for (l = 1; l < n; l++)
		for (j = 0; j < n; j++) {
			if (j == i)
				continue;
			c.push_back(a[l][j]);
			k++;
			if (k == (n - 1)) {
				b.push_back(c);
				c.clear();
				h++;
				k = 0;
			}
		}
	return b;
}// end function

 //---------------------------------------------------
 //	calculate determinte of matrix
float RANSAC::det(vector<vector<double>> a, int n) {
	int i;
	vector<vector<double>> b;
	float sum = 0;
	if (n == 1)
		return a[0][0];
	else if (n == 2)
		return (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
	else
		for (i = 0; i < n; i++) {
			b = minor(a, i, n);	// read function
			sum = (float)(sum + a[0][i] * pow(-1, i)*det(b, (n - 1)));	// read function	// sum = determinte matrix
		}
	return sum;
}// end function

 //---------------------------------------------------
 //	calculate transpose of matrix
vector<vector<double>> RANSAC::transpose(vector<vector<double>> c, int n, float det) {
	vector<vector<double>> d;
	int i, j;
	float b[500][500];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			b[i][j] = c[j][i];

	for (i = 0; i < n; i++) {
		vector<double > e;
		for (j = 0; j < n; j++)
			e.push_back(b[i][j] / det);	// array d[][] = inverse matrix
		d.push_back(e);
	}


	return d;
}// end function

 //---------------------------------------------------
 //	calculate cofactor of matrix
vector<vector<double>> RANSAC::cofactor(vector<vector<double>> a, int n, float determinte) {
	vector<vector<double>> b, c;
	int l, h, m, k, i, j;
	vector<double> e;
	for (h = 0; h < n; h++) {
		vector<double> f;
		for (l = 0; l < n; l++) {
			m = 0;
			k = 0;
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (i != h && j != l) {
						e.push_back(a[i][j]);
						if (k < (n - 2))
							k++;
						else {
							b.push_back(e);
							e.clear();
							k = 0;
							m++;
						}
					}
				}
			}
			f.push_back((float)pow(-1, (h + l))*det(b, (n - 1)));	// c = cofactor Matrix
			b.clear();
		}
		c.push_back(f);
	}
	return transpose(c, n, determinte);	// read function
}// end function

 //---------------------------------------------------
 //	calculate inverse of matrix
vector<vector<double>> RANSAC::inverse(vector<vector<double>> a, int n, float det) {
	vector<vector<double>> d;
	if (det == 0)
		cout << "\nInverse of Entered Matrix is not possible\n";
	else if (n == 1)
		d[0][0] = 1;
	else
		d = cofactor(a,  n, det);	// read function
	return d;
}// end function