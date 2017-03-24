#include "Fit_Equation.h"
#include <math.h>


//int array to Site array.
vector<Site> Fit_Equation::Set_Site_Array(double site_array[]) {
	int Array_Number = sizeof(site_array);
	vector<Site> Site_Array;
	for (int Index_Number = 0; Index_Number < Array_Number / 2; Index_Number++) {
		Site Temp_Site;
		Temp_Site.set_x(site_array[2 * Index_Number]);
		Temp_Site.set_y(site_array[2 * Index_Number + 1]);
		Site_Array.push_back(Temp_Site);
	}
	return  Site_Array;
}

vector<vector<double>> Fit_Equation::Equation_Array(vector<Site> Site_Array, int Parameter_Number) {
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

vector<vector<double>> Fit_Equation::Traspose_Array(vector<vector<double>> Target_Array) {
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

vector<vector<double>> Fit_Equation::Dot_Array(vector<vector<double>> Array_One, vector<vector<double>> Array_Two) {
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

vector<vector<double>> Fit_Equation::Inverse_Array(vector<vector<double>> Target_Array) {
	return Inverse(Target_Array, Target_Array.size(), Det(Target_Array, Target_Array.size()));
}

//	calculate minor of matrix OR build new matrix : k-had = minor
vector<vector<double>> Fit_Equation::Minor(vector<vector<double>> Target_Array, int Index_Matrix, int Matrix_Number) {
	vector<vector<double>> Minor_Array;
	int h = 0, k = 0;
	vector<double> Minor_row;
	for (int l = 1; l < Matrix_Number; l++)
		for (int j = 0; j < Matrix_Number; j++) {
			if (j == Index_Matrix)
				continue;
			Minor_row.push_back(Target_Array[l][j]);
			k++;

			if (k == (Matrix_Number - 1)) {
				Minor_Array.push_back(Minor_row);
				Minor_row.clear();
				h++;
				k = 0;
			}
		}
	return Minor_Array;
}

//	calculate determinte of matrix
float Fit_Equation::Det(vector<vector<double>> Target_Array, int Matrix_Number) {
	vector<vector<double>> Minor_Array;
	float Sum = 0;
	if (Matrix_Number == 1)
		return Target_Array[0][0];
	else if (Matrix_Number == 2)
		return (Target_Array[0][0] * Target_Array[1][1] - Target_Array[0][1] * Target_Array[1][0]);
	else
		for (int Index_Matrix = 0; Index_Matrix < Matrix_Number; Index_Matrix++) {
			Minor_Array = Minor(Target_Array, Index_Matrix, Matrix_Number);
			Sum = (float)(Sum + Target_Array[0][Index_Matrix] * pow(-1, Index_Matrix)*Det(Minor_Array, (Matrix_Number - 1)));	// Sum = determinte matrix
		}
	return Sum;
}

//	calculate transpose of matrix
vector<vector<double>> Fit_Equation::Transpose(vector<vector<double>> Target_Array, int Matrix_Number, float Det) {
	vector<vector<double>> Transpose_Array;
	/*float b[500][500];
	for (int Index_row = 0; Index_row < Matrix_Number; Index_row++)
		for (int Index_column = 0; Index_column < Matrix_Number; Index_column++)
			b[Index_row][Index_column] = Target_Array[Index_column][Index_row];*/

	for (int Index_row = 0; Index_row < Matrix_Number; Index_row++) {
		vector<double > Transpose_row;
		for (int Index_column = 0; Index_column < Matrix_Number; Index_column++)
			Transpose_row.push_back(/*b*/Target_Array[Index_row][Index_column] / Det);	// array d[][] = inverse matrix
		Transpose_Array.push_back(Transpose_row);
	}
	return Transpose_Array;
}

//	calculate cofactor of matrix
vector<vector<double>> Fit_Equation::Cofactor(vector<vector<double>> Target_Array, int Matrix_Number, float Determinte) {
	vector<vector<double>> b, c;
	int m, k;
	vector<double> e;
	for (int Index_row = 0; Index_row < Matrix_Number; Index_row++) {
		vector<double> f;
		for (int Index_column = 0; Index_column < Matrix_Number; Index_column++) {
			m = 0;
			k = 0;
			for (int i = 0; i < Matrix_Number; i++)
				for (int j = 0; j < Matrix_Number; j++)
					if (i != Index_row && j != Index_column) {
						e.push_back(Target_Array[i][j]);
						if (k < (Matrix_Number - 2))
							k++;
						else {
							b.push_back(e);
							e.clear();
							k = 0;
							m++;
						}
					}
			f.push_back((float)pow(-1, (Index_row + Index_column))*Det(b, (Matrix_Number - 1)));	// c = cofactor Matrix
			b.clear();
		}
		c.push_back(f);
	}
	return Transpose(c, Matrix_Number, Determinte);
}

//	calculate inverse of matrix
vector<vector<double>> Fit_Equation::Inverse(vector<vector<double>> Target_Array, int Matrix_Number, float Det) {
	vector<vector<double>> Inverse_Array;
	if (Det == 0)
		cout << "\nInverse of Entered Matrix is not possible\n";
	else if (Matrix_Number == 1)
		Inverse_Array[0][0] = 1;
	else
		Inverse_Array = Cofactor(Target_Array, Matrix_Number, Det);
	return Inverse_Array;
}

vector<double> Fit_Equation::Handle(double site_array[], int Parameter_Number) {
	//X
	vector<vector<double>> X_Array = Equation_Array(Set_Site_Array(site_array), Parameter_Number);
	//X^T
	vector<vector<double>> Transpose_X = Traspose_Array(X_Array);
	//Y
	vector<vector<double>> Y_Array;
	for (int Index_Y = 0; Index_Y < sizeof(site_array)/2; Index_Y++) {
		vector<double> One_Y{ site_array[2 * Index_Y + 1] };
		Y_Array.push_back(One_Y);
	}
	
	//P=(X^T¡DX)^-1¡DX^T¡DY
	vector<vector<double>> Parameter = Dot_Array(Dot_Array(Inverse_Array(Dot_Array(Transpose_X, X_Array)), Transpose_X), Y_Array);
	vector<double > Parameter_Array;
	for (int i = 0; i < Parameter.size(); i++) {
		Parameter_Array.push_back(Parameter[i][0]);
	}
	return Parameter_Array;
}

vector<double> Fit_Equation::Handle(vector<Site> site_array, int Parameter_Number) {
	//X
	vector<vector<double>> X_Array = Equation_Array(site_array, Parameter_Number);
	//X^T
	vector<vector<double>> Transpose_X = Traspose_Array(X_Array);
	//Y
	vector<vector<double>> Y_Array;
	for (int Index_Site = 0; Index_Site < site_array.size(); Index_Site++) {
		vector<double> One_Y;
		One_Y.push_back(site_array[Index_Site].get_y());
		Y_Array.push_back(One_Y);
	}

	//P=(X^T¡DX)^-1¡DX^T¡DY
	vector<vector<double>> Parameter = Dot_Array(Dot_Array(Inverse_Array(Dot_Array(Transpose_X, X_Array)), Transpose_X), Y_Array);
	vector<double > Parameter_Array;
	for (int i = 0; i < Parameter.size(); i++) {
		Parameter_Array.push_back(Parameter[i][0]);
	}
	return Parameter_Array;
}


