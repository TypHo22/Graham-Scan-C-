#include <graham.h>

double polyCheck(Point& P1, Point& P2, Point& P3)
{
	double p1p2x, p1p2y, p1p2z; //mathematical x,y,z vector
	double p1p3x, p1p3y, p1p3z; //mathematical x,y,z vector
	p1p2x = P2.X - P1.X;
	p1p2y = P2.Y - P1.Y;
	p1p2z = 0; //there will be no z coordinate
	p1p3x = P3.X - P1.X;
	p1p3y = P3.Y - P1.Y;
	p1p3z = 0; //there will be no z coordinate
	double crossX, crossY, crossZ;
	crossX = p1p2y * p1p3z - p1p2z * p1p3y;
	crossY = p1p2z * p1p3x - p1p2x * p1p3z;
	crossZ = p1p2x * p1p3y - p1p2y * p1p3x;
	double PolyVal = crossX + crossY + crossZ;
	return  PolyVal;
}

void removeElements(std::map<double, unsigned int>& ArcValues, int& k)
{
	std::map<double, unsigned int>::iterator it;
	it = ArcValues.begin();

	int b = 0;
	while (b < k + 1)
	{
		it++;
		b++;
	}
	ArcValues.erase(it);
}

unsigned int getIndice(std::map<double, unsigned int>& ArcValues, unsigned int k, int offSet)
{
	std::map<double, unsigned int>::iterator it;
	it = ArcValues.begin();
	for (unsigned int a = 0; a < k + offSet; a++)
	{
		it++;
	}

	unsigned int indi = it->second;
	return indi;
}