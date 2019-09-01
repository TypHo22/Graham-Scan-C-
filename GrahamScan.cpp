#include <iostream>
#include <graham.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace std;

int main()
{
	srand(time(NULL));
	//vector which holds the random values
	vector<vector<double>> xyPair;
	xyPair.push_back(vector<double>());
	xyPair.push_back(vector<double>());
	  for(unsigned int a=0 ; a<50 ; a++)
	  {
		  xyPair[0].push_back(rand() );
		  xyPair[1].push_back(rand() );
	  }
	

	//vector which holds the convex Hull
	vector<vector<double>> convexHull;
	convexHull.push_back(vector<double>());
	convexHull.push_back(vector<double>());
	time_t tstart, tend;
	tstart = time(0);

	grahamScan<double>(xyPair, convexHull);

	tend = time(0);
	cout << "It took " << difftime(tend, tstart) << " second(s)." << endl;
	cout << "Amount of convex Hull points:" << convexHull[0].size() - 1  << endl;

	for (unsigned int a = 0; a < convexHull[0].size(); a++)
	{
		cout << "x: " << convexHull[0][a] << " y: " << convexHull[1][a] << endl;
	}
	// do File writing for further Processing. Can be removed if not needed
	std::ofstream scatterPointsFile;
	scatterPointsFile.open("scatterPoints.csv");
	for (unsigned int a = 0; a < xyPair[0].size(); a++)
	{
		scatterPointsFile << xyPair[0][a] << ";" << xyPair[1][a] << "\n";
	}
	scatterPointsFile.close();

	std::ofstream hullPointsFile;
	hullPointsFile.open("hull.csv");
	for (unsigned int a = 0; a < convexHull[0].size(); a++)
	{
		hullPointsFile << convexHull[0][a] << ";" << convexHull[1][a] << "\n";
	}
	hullPointsFile.close();

}
