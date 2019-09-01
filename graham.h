#pragma once
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>

const double M_PI = 3.14159265358979323846;

/**
* @brief   Provides the connecting of two ways
*    Algorithm for calculating a convex hull around a point cloud Graham's scan is a method of finding the convex hull of a finite set of points in the plane with time complexity O(n log n).
*    It is named after Ronald Graham, who published the original algorithm in 1972.
*    The algorithm finds all vertices of the convex hull ordered along its boundary. It uses a stack to detect and remove concavities in the boundary efficiently.
*
*    Further Informations: https://en.wikipedia.org/wiki/Graham_scan https://www.youtube.com/watch?v=VP9ylElm1yY
*    If there is any Troubleshooting with this content I have done this Algorithm before: https://github.com/TypHo22/Graham-Scan (written in MATLAB)
*
*    Input grahamscan: std::vector<std::vector<T>> &xyPointcloud,std::vector<std::vector<T>> &convexHull
*    T means any numeric Value.
*    My convention:
*    xyPointcloud[0] and convexHull[0] holds all the x-Values
*    xyPointcloud[1] and convexHull[1] holds all the y-Values,
*
* @author  Andreas Bernatzky
*
* @date    23 April 2019
*/

struct Point
{
	double X;
	double Y;
};
/*
polyCheck: One corner pi in a polygon course p0p1, p1p2, ...,
pn-1p0 means convex, if for the left angle ? between the edges
pi-1pi and pipi+1 is 0° less than or equal to ? < 180° (i-1 and i+1 are calculated modulo n).
Otherwise the corner is called concave.
concave = clockwise rotation
convex = counterclockwise rotation
I have chosen to determine the cross product,
whether the polygon course is concave or convex.  If the sum of the cross product does have a negative value
then the polygon course is concave. With a positive value for polycrumming the result is
a convex polygon course
PolyVal < 0 = concave = clockwise rotation
PolyVal > 0 = convex = left rotation
You could also check for convex or concav via the distances P0P1 and P0P2
Instead of turning to the left or to the right, it is also possible to use a
*/
double polyCheck(Point& P1, Point& P2, Point& P3);

//fill the points for a polycourse
template <typename T>
static void fillPoint(Point& mP, T& x, T& y)
{
	mP.X = x;
	mP.Y = y;
}

//remove Element from ArcValues because it makes the poly course concave
void removeElements(std::map<double, unsigned int>& ArcValues, int& k);

// get the indice for the next P1,P2,P3 xycoordinates
//Arcvalues has as key the angle values between Z and the respective point.
//The vector position is stored in arcvalues for the value.
unsigned int getIndice(std::map<double, unsigned int>& ArcValues, unsigned int k, int offSet);

template <typename T>
static void calcArc(std::map<double, unsigned int>& ArcValues, std::vector<std::vector<T>>& xyPointcloud, Point Z)
{
	const double toDeg = 180 / M_PI;
	double dX, dY, ArcRad, ArcDeg;
	for (unsigned int a = 0; a < xyPointcloud[0].size(); a++)
	{
		dX = xyPointcloud[0][a] - Z.X;
		dY = xyPointcloud[1][a] - Z.Y;
		ArcRad = atan2(dY, dX);
		if (ArcRad < 0)
		{
			ArcDeg = 360 - ((abs(ArcRad) / M_PI) * 180);
		}
		else
		{
			ArcDeg = ArcRad * toDeg;
		}
		//check if the key already exists!
	
		std::map<double, unsigned int>::iterator it = ArcValues.find(ArcDeg);//check if there exists points with the same key (arc value)

		if (it == ArcValues.end())// the key is unique you do not have to check for the distances
		{
			ArcValues.insert(std::pair<double, unsigned int>(ArcDeg, a));
		}
		else //the mapkey is not unique! we have to check the distances 
		{
			double actDistance = sqrt(dX * dX + dY * dY); //distance of the actual point (this point must be checked if its distance is greater than the actual point with the same 
												          //distance.
			double dXcomp = xyPointcloud[0][it->second] - Z.X;//get the x value and calculate dX from Z Point
			double dYcomp = xyPointcloud[1][it->second] - Z.Y;//get the y value and calculate dY from Z Point
			double compDistance = sqrt(dXcomp * dXcomp + dYcomp * dYcomp); //distance of the point which gets compared to the actualDistance point

			if (actDistance >= compDistance) //the actual Point is more likely to be a hullpoint because its distance is much greater
			{
				it->second = a;//will overwrite the old map value with the indice of the new point		
			}
		}
		
		
	}
	
}


template <typename T>
static void grahamScan(std::vector<std::vector<T>>& xyPointcloud, std::vector<std::vector<T>>& convexHull)
{

	auto xMax = max_element(std::begin(xyPointcloud[0]), std::end(xyPointcloud[0]));
	auto xMin = min_element(std::begin(xyPointcloud[0]), std::end(xyPointcloud[0]));
	auto yMax = max_element(std::begin(xyPointcloud[1]), std::end(xyPointcloud[1]));
	auto yMin = min_element(std::begin(xyPointcloud[1]), std::end(xyPointcloud[1]));
	// calculate the center point of the pointcloud
	Point Z;
	Z.X = *xMin + (*xMax - *xMin) / 2; //x-Coordinate of the center point
	Z.Y = *yMin + (*yMax - *yMin) / 2; //y-Coordinate of the center point

	//Arcvalues has as key the angle values between Z and the respective point.
	//The vector position is stored in arcvalues for the value.
	std::map<double, unsigned int> ArcValues; //holds the Arc between every point of the xyCloud and the center point (Z)
	calcArc(ArcValues, xyPointcloud, Z); //calculate the arcs

// map values are sorted ascending for the keyValue (important for the algorithm to work)
//Calculate convexHull
	bool noMoreRemove = false; //condition for staying in while-Loop. the algorithm runs until there isn't anything more to remove
	int k = 0; //inkrement of while-loop
	unsigned int cas = 0; //case for map endings
	unsigned int mapSize = ArcValues.size(); // actual size of map
	std::map<double, unsigned int>::iterator it; //iterator object vor ArcValues
	unsigned int indi; //takes the indice of Arcvalue (second)
	Point P1, P2, P3; //for polygon course. PolygonCourse consists of P1P2P3
	double polyTest; // result of polyCheck if polyTest > 0 means convex polygon course, polyTest < 0 means concav polygon course
	while (noMoreRemove == false)
	{
		if (k < mapSize - 2)//standard case far away from the vector end
		{
			indi = getIndice(ArcValues, k, 0);
			fillPoint(P1, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			indi = getIndice(ArcValues, k, 1);
			fillPoint(P2, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			indi = getIndice(ArcValues, k, 2);
			fillPoint(P3, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			cas = 0;
		}
		else if (k == mapSize - 2)//k is the penultimate element must now be k+1 = last vector element and k+2 == 1st vector element
		{
			indi = getIndice(ArcValues, k, 0);
			fillPoint(P1, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			indi = getIndice(ArcValues, k, 1);
			fillPoint(P2, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			indi = getIndice(ArcValues, 0, 0);
			fillPoint(P3, xyPointcloud[0][1], xyPointcloud[1][1]);
			cas = 0;
		}
		else if (k == mapSize)//k is the last element must now be k+1 == 1st vector element and k+2 == 2nd vector element
		{
			indi = getIndice(ArcValues, k, - 1);
			fillPoint(P1, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			indi = getIndice(ArcValues, 0, 0);
			fillPoint(P2, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			indi = getIndice(ArcValues, 1, 0);
			fillPoint(P3, xyPointcloud[0][indi], xyPointcloud[1][indi]);
			cas = 1;
		}
		//check if current poly course is concave or convex
		polyTest = polyCheck(P1, P2, P3);

		if (polyTest <= 0 && cas == 0)
		{
			removeElements(ArcValues, k);

			k = k - 3;
			if (k < -1)
			{
				k = -1;
			}
		}
		else if (polyTest < 0 && cas == 1)
		{

			it = ArcValues.begin();
			ArcValues.erase(it);//remove first element

			k = k - 3;
			if (k < -1)
			{
				k = -1;
			}
		}
		else if (k == mapSize && polyTest > 0)
		{
			noMoreRemove = true;
			break;
		}
		mapSize = ArcValues.size();
		k++;
	}
	// create the convexPolygon
	it = ArcValues.begin();
	for (unsigned int a = 0; a < ArcValues.size(); a++)
	{
		convexHull[0].push_back(xyPointcloud[0][it->second]);
		convexHull[1].push_back(xyPointcloud[1][it->second]);
		it++;
	}
}
