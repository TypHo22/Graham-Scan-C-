#include <QCoreApplication>
#include <grahamscan.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <ctime>
using namespace std;

int main(int argc, char *argv[])
{
    srand(time(NULL));
    //vector which holds the random values
    vector<vector<double>> xyPair;
    xyPair.push_back(vector<double>());
    xyPair.push_back(vector<double>());
  /*  for(uint a=0;a<2000;a++)
    {
        xyPair[0].push_back(rand());
        xyPair[1].push_back(rand());
    }
*/
    xyPair[0].push_back(0.5);
    xyPair[0].push_back(1);
    xyPair[0].push_back(2);
    xyPair[0].push_back(3);
    xyPair[0].push_back(3);

    xyPair[1].push_back(2);
    xyPair[1].push_back(1);
    xyPair[1].push_back(4);
    xyPair[1].push_back(9);
    xyPair[1].push_back(2);
    //vector which holds the convex Hull
    vector<vector<double>> convexHull;
    convexHull.push_back(vector<double>());
    convexHull.push_back(vector<double>());
    time_t tstart, tend;
    tstart = time(0);

    grahamScan<double>(xyPair,convexHull);

   tend = time(0);
   cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    QCoreApplication a(argc, argv);
    return a.exec();
}
