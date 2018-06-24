#include <iostream>
#include "Coord.h"
using namespace std;

Coord::Coord() {
	x = 0;
	y = 0;
	z = 0;
}
Coord::Coord(int i, int j, int k) {
	x = i;
	y = j;
	z = k;
}
Coord::Coord(double i, double j, double k) {
	x = i;
	y = j;
	z = k;
}
Coord Coord::operator+(const Coord& c) const {
	Coord cNew;
	cNew.x = c.x + x;
	cNew.y = c.y + y;
	cNew.z = c.z + z;
	return cNew;
}
Coord Coord::operator/(const Coord& c) const {
	Coord cNew;
	cNew.x = x / c.x;
	cNew.y = y / c.y;
	cNew.z = z / c.z;
	return cNew;
}
Coord Coord::operator*(const Coord& c) const {
	Coord cNew;
	cNew.x = c.x * x;
	cNew.y = c.y * y;
	cNew.z = c.z * z;
	return cNew;
}
Coord Coord::operator-(const Coord& c) const {
	Coord cNew;
	cNew.x = x - c.x;
	cNew.y = y - c.y;
	cNew.z = z - c.z;
	return cNew;
}
Coord Coord::operator=(const Coord& c) {
	x = c.x;
	y = c.y;
	z = c.z;
	return *this;
}
Coord Coord::operator/(const double& i) const {
	Coord cNew;
	cNew.x = x / i;
	cNew.y = y / i;
	cNew.z = z / i;
	return cNew;
}
Coord Coord::operator*(const double& i) const {
	Coord cNew;
	cNew.x = x * i;
	cNew.y = y * i;
	cNew.z = z * i;
	return cNew;
}
Coord Coord::operator+(const double& i) const {
	Coord cNew;
	cNew.x = x + i;
	cNew.y = y + i;
	cNew.z = z + i;
	return cNew;
}
Coord Coord::operator-(const double& i) const {
	Coord cNew;
	cNew.x = x - i;
	cNew.y = y - i;
	cNew.z = z - i;
	return cNew;
}

