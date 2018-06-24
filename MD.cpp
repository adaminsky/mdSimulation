#include <iostream>
#include <cstdlib>
#include <cmath>
#include "Coord.h"
#include "gnuplot_i.hpp"
using namespace std;

double tmp = 240;
double etot;
const double DT = 1e-2;
const double DENSITY_FACTOR = 1.709975947;
const double X_SIZE = 1/(10.25985568);
const double PI = 3.14159;


Coord velocities[216];
Coord x[216];
Coord xprev[216];
Coord f[216];
double en;

// Set the initial velocities of each atom so that the velocity CM is 0 and 
// determine the initial kinetic energy of the system.
void init() {
	srand(time(0));

	Coord sumV;
	Coord sumV2;

	// Assigning the particles to a lattice.
	for (int i = 0; i < 216; i++) {
		x[i] = Coord((i/36)*DENSITY_FACTOR, ((i%36)/6)*DENSITY_FACTOR, (i%6)*DENSITY_FACTOR);
		velocities[i] = Coord(((double)rand() / RAND_MAX) - .5, ((double)rand() / RAND_MAX) - .5, ((double)rand() / RAND_MAX) - .5);
		sumV = sumV + velocities[i];
		sumV2 = sumV2 + (velocities[i]*velocities[i]);
	}	
	sumV = sumV / 216;
	sumV2 = sumV2 / 216;

	// Velocity scale factor.
	double fsx = sqrt(3*tmp/sumV2.x);
	double fsy = sqrt(3*tmp/sumV2.y);
	double fsz = sqrt(3*tmp/sumV2.z);
	Coord fs = Coord(fsx, fsy, fsz);

	for (int i = 0; i < 216; i++) {
		// Subtracting the CM velocity from each velocity to make the
		// CM velocity zero.
		velocities[i] = (velocities[i] - sumV);

		// Previous locations are calculated based on the current velocity.
		xprev[i] = x[i] - velocities[i]*DT;
	}
}

void force() {
	en = 0;
	// Setting forces to zero.
	for (int i = 0; i < 216; i++) {
		f[i].x = 0;
		f[i].y = 0;
		f[i].z = 0;
	}
	Coord xrP;
	double r2;
	// Looping over all particle pairs.
	for (int i = 0; i < 215; i++) {
		for (int j = i+1; j < 216; j++) {
			xrP = (x[i] - x[j]);

			// periodic boundary conditions:
			xrP.x -= static_cast<int>(xrP.x * X_SIZE + 0.5) / X_SIZE;
			xrP.y -= static_cast<int>(xrP.y * X_SIZE + 0.5) / X_SIZE;
			xrP.z -= static_cast<int>(xrP.z * X_SIZE + 0.5) / X_SIZE;

			r2 = (xrP.x*xrP.x)+(xrP.y*xrP.y)+(xrP.z*xrP.z);
			
			// Checking for cutoff.
			if (r2 < 26.31615964) {
				double r2i = 1/r2;
				double r6i = r2i*r2i*r2i;
				double ff = 48*r2i*r6i*(r6i-0.5);
				f[i] = f[i] + xrP*ff;
				f[j] = f[j] - xrP*ff;
				en = en + 4*r6i*(r6i-1)-(-2.19466695e-4);
			}
			
		}	
	}
}
void integrate(double t) {
	Coord sumv;
	Coord sumv2;
	for (int i = 0; i < 216; i++) {
		Coord xx, vi;
		xx = x[i]*2 - xprev[i] + f[i]*(DT)*(DT);
		vi = (xx - xprev[i]) / (2*(DT));
		sumv = sumv + vi;
		sumv2 = sumv2 + (vi*vi);
		xprev[i] = x[i];
		x[i] = xx;
	}
	tmp = sqrt((sumv2.x*sumv2.x)+(sumv2.y*sumv2.y)+(sumv2.z*sumv2.z)) / (3*216);
	etot = (en + 0.5*sqrt((sumv2.x*sumv2.x)+(sumv2.y*sumv2.y)+(sumv2.z*sumv2.z))) / 216;
}

void radDistro() {
	int ngr = 0;
	double delg = 10.2598 / (2*100);
	double g[100] = {0};
	ngr++;
	for (int i = 0; i < 215; i++) {
		for (int j = i + 1; j < 216; j++) {
			Coord xr;
			xr = x[i] - x[j];
			xr.x -= static_cast<int>(xr.x * X_SIZE + 0.5) / X_SIZE;
			xr.y -= static_cast<int>(xr.y * X_SIZE + 0.5) / X_SIZE;
			xr.z -= static_cast<int>(xr.z * X_SIZE + 0.5) / X_SIZE;
			double r = sqrt((xr.x*xr.x)+(xr.y*xr.y)+(xr.z*xr.z));
			if (r < (1/X_SIZE)/2) {
				int ig = static_cast<int>(r/delg);
				g[ig] = g[ig] + 2;
			}
		}
	}
	for (int i = 0; i < 100; i++) {
		double r = delg*(i+.5);
		double vb = ((i+1)*(i+1)*(i+1) - i*i*i)*delg*delg*delg;
		double nid = (4.0/3)*PI*vb*.2;
		g[i] = g[i] / (216*nid);
	}
	for (int i = 0; i < 100; i++) {
		cout << g[i] << endl;
	}
}

int main() {
	Gnuplot g("points");
	vector<double> x1, y1, z1;
	init();
//	for (int i = 0; i < 216; i++) {
//		cout << x[i].x << ", " << x[i].y << ", " << x[i].z << "   ";
//		if ((i+1) % 6 == 0)
//			cout << endl;
//	}
//	for (int i = 0; i < 216; i++) {
//		cout << velocities[i].x << ", " << velocities[i].y << ", " << velocities[i].z << "   ";
//		if ((i+1) % 6 == 0)
//			cout << endl;
//	}
	for (int i = 0; i < 216; i++) {
		x1.push_back(x[i].x);
		y1.push_back(x[i].y);
		z1.push_back(x[i].z);
	}
	g.reset_plot();
	g.unset_grid();
	g.plot_xyz(x1,y1,z1, "Test");
	//cin.get();
	int index = 0;
	Coord avV;
	for (double t = 0; t < 1e1; t+=DT) {
		force();
		integrate(t);
		//cout << etot << endl;
		//cout << tmp << endl;
		x1.clear();
		y1.clear();
		z1.clear();
		if (index >= 61)
			g.remove_oldest_tmpfile();
		for (int i = 0; i < 216; i++) {
			x1.push_back(x[i].x);
			y1.push_back(x[i].y);
			z1.push_back(x[i].z);
		}
		g.reset_plot();
		g.unset_grid();
		g.plot_xyz(x1,y1,z1, "Test");
		usleep(10000);
		index++;
	}
	for (int i = 0; i < 216; i++) {
		cout << x[i].x << ", " << x[i].y << ", " << x[i].z << "   ";
		if ((i+1) % 6 == 0)
			cout << endl;
	}

	radDistro();

}
