#include "Fit_Equation.h"
#include <time.h>
#include <opencv\highgui.h>
#include <opencv2\opencv.hpp>
#include "RSC.h"

int const Parameter_Number = 3;
using namespace cv;
double site[12];
double x[] = { 1,6,-5,-7,3,-10 };
#define One_Line_SIZE 100
#define Site_Number 5 
char Y_Site[One_Line_SIZE];

int main() {
	for (int i = 0; i < 6; i++) {
		site[i * 2] = x[i];
		site[i * 2 + 1] =  + (-4 * x[i] * x[i]) + 7 * x[i] + 13;
	}
	vector<Site> Site_Array;
	fstream fin;
	fin.open("D:\\Visual Studio 2015 專案\\其他\\Random_Sample_Consensus\\07_000028.txt",ios::in);
	Mat Image = imread("D:\\Visual Studio 2015 專案\\其他\\Random_Sample_Consensus\\ti.png");
	int Index_X = 0;
	while (fin.getline(Y_Site, sizeof(Y_Site), '\n')) {	
		Site Temp_Site;
		Temp_Site.set_x(Index_X-110);
		Temp_Site.set_y(atof(Y_Site)/50);
		Site_Array.push_back(Temp_Site);
		Index_X++;
		//cout << "X = " <<Temp_Site.get_x() << "Y = " << Temp_Site.get_y() << endl;
	}

	vector<Site> Fitting_Site;
	srand((unsigned)time(NULL));
	map<int,double> Random_number;
	do {
		int Ran = rand() % 40;
		if (Random_number.find(Ran) == Random_number.end())
			Random_number.insert(pair<int, double>(Ran, (double)Ran));
	} while (Random_number.size() <= Site_Number -1 );

	for (auto i = Random_number.begin(); i != Random_number.end(); ++i) {
		//Fitting_Site.push_back(Site_Array[130+(rand() % 40)]);
		//Fitting_Site.push_back(Site_Array[i]);
		Fitting_Site.push_back(Site_Array[130+(i -> first)]);
	}
	Fit_Equation RANSAC;
	vector<double> Function_Parameter = RANSAC.Handle(Fitting_Site, Parameter_Number);


	for (int i = 0; i < Function_Parameter.size(); i++)
		cout << Function_Parameter[i] << " , ";
	for (int site_X = 0; site_X < 450; site_X++) {
		int Y = 0;
		for (int Index_Pow_X = Function_Parameter.size()-1; Index_Pow_X >=0 ; Index_Pow_X--) {
			Y+= Function_Parameter[Function_Parameter.size()- Index_Pow_X -1]*pow(site_X, Index_Pow_X);
		}
		int X = 70 + site_X  +110;
		Y = 620 - Y;
		if( Y > 20 && Y < 620 )
			if( X > 70 && X < 520)
				Image.at<Vec3b>(Point(X, Y)) = { 0,0,255 };
	}
	imshow("asd", Image);
	imwrite("D:\\Visual Studio 2015 專案\\其他\\Random_Sample_Consensus\\123.png", Image);


	RSC test_RSC;
	double Average_Distance = (test_RSC.Total_Distance(Site_Array, Function_Parameter))/(290-110);
	pair<vector<double>, double> KK = test_RSC.Fittest_Parameter(300, Site_Array, Average_Distance);

	cout << KK.second << endl;
	waitKey(0);
	return 0;
}