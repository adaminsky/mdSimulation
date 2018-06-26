#include <iostream>
#include <cstdlib>
#include <cmath>
#include "Coord.h"
#include "gnuplot_i.hpp"
using namespace std;

double tmp = 2;
double etot;
// Time step.
const double DT = 1e-3;
// Distance between two particles in the original configuration (measured in A).
double DENSITY_FACTOR; // = 1.077217345;//1.709975947;
// Total number of particles in the system.
const int NPART = 216;
double SIDE_LENGTH; // = 6*DENSITY_FACTOR;
// Reciprocal of the side length of the original cube layout.
double X_SIZE; // = 1/SIDE_LENGTH;
const double PI = 3.14159;
double ECUT;// = 6547.0464;//15.97439995;
double ngr = 0;
double g[100] = {0};


Coord velocities[NPART];
Coord x[NPART];
Coord xprev[NPART];
Coord f[NPART];

double en;

// Set the initial velocities of each atom so that the velocity CM is 0 and 
// determine the initial kinetic energy of the system.
void init() {
	srand(time(0));

	Coord sumV;
	Coord sumV2;

	// Assigning the particles to a lattice.
	for (int i = 0; i < NPART; i++) {
		x[i] = Coord((i/36)*DENSITY_FACTOR, ((i%36)/6)*DENSITY_FACTOR, (i%6)*DENSITY_FACTOR);
		velocities[i] = Coord(((double)rand() / RAND_MAX) - .5, ((double)rand() / RAND_MAX) - .5, ((double)rand() / RAND_MAX) - .5);
		sumV = sumV + velocities[i];
		sumV2 = sumV2 + (velocities[i]*velocities[i]);
	}	
	sumV = sumV / NPART;
	sumV2 = sumV2 / NPART;

	// Velocity scale factor.
	double fsx = sqrt(3*tmp/sumV2.x);
	double fsy = sqrt(3*tmp/sumV2.y);
	double fsz = sqrt(3*tmp/sumV2.z);
	Coord fs = Coord(fsx, fsy, fsz);

	for (int i = 0; i < NPART; i++) {
		// Subtracting the CM velocity from each velocity to make the
		// CM velocity zero.
		velocities[i] = (velocities[i] - sumV)*fs;

		// Previous locations are calculated based on the current velocity.
		xprev[i] = x[i] - velocities[i]*DT;
	}
}

void force() {
	en = 0;
	// Setting forces to zero.
	for (int i = 0; i < NPART; i++) {
		f[i].x = 0;
		f[i].y = 0;
		f[i].z = 0;
	}
	Coord xrP;
	double r2;
	// Looping over all particle pairs.
	for (int i = 0; i < NPART - 1; i++) {
		for (int j = i+1; j < NPART; j++) {
			xrP = (x[i] - x[j]);

			// periodic boundary conditions:
		//	xrP.x -= static_cast<int>(xrP.x * X_SIZE + 0.5) / X_SIZE;
		//	xrP.y -= static_cast<int>(xrP.y * X_SIZE + 0.5) / X_SIZE;
		//	xrP.z -= static_cast<int>(xrP.z * X_SIZE + 0.5) / X_SIZE;
			xrP.x -= round(xrP.x*X_SIZE) / X_SIZE;
			xrP.y -= round(xrP.y*X_SIZE) / X_SIZE;
			xrP.z -= round(xrP.z*X_SIZE) / X_SIZE;

			r2 = (xrP.x*xrP.x)+(xrP.y*xrP.y)+(xrP.z*xrP.z);
			
			// Checking for cutoff.
			if (r2 < (SIDE_LENGTH/2)*(SIDE_LENGTH/2)) {
				double r2i = 1./r2;
				double r6i = r2i*r2i*r2i;
				double ff = 48*r2i*r6i*(r6i-0.5);
				f[i] = f[i] + (xrP*ff);
				f[j] = f[j] - (xrP*ff);
				en = en + 4*r6i*(r6i-1)-(ECUT);
			}
			
		}	
	}
}
void integrate(double t) {
	Coord sumv;
	Coord sumv2;
	for (int i = 0; i < NPART; i++) {
		Coord xx, vi;
		xx = x[i]*2;
		xx = xx - xprev[i];
		xx = xx + (f[i]*DT*DT);

		vi = (xx - xprev[i]) / (2*(DT));
		sumv = sumv + vi;
		sumv2 = sumv2 + (vi*vi);
		xprev[i] = x[i];
		x[i] = xx;
		//	x[i].x -= floor(x[i].x / SIDE_LENGTH)*SIDE_LENGTH;
		//	x[i].y -= floor(x[i].y / SIDE_LENGTH)*SIDE_LENGTH;
		//	x[i].z -= floor(x[i].z / SIDE_LENGTH)*SIDE_LENGTH;
	}
	tmp = sqrt((sumv2.x*sumv2.x)+(sumv2.y*sumv2.y)+(sumv2.z*sumv2.z)) / (3*NPART);
	etot = (en + 0.5*sqrt((sumv2.x*sumv2.x)+(sumv2.y*sumv2.y)+(sumv2.z*sumv2.z))) / NPART;
}

// Pass the number of buckets, at most 100.
void radDistro(int bucketn, double d, int i) {
	double delg = (SIDE_LENGTH) / (2*bucketn);
	ngr+=1;
	for (int i = 0; i < NPART - 1; i++) {
		for (int j = i + 1; j < NPART; j++) {
			Coord xr;
			xr = x[i] - x[j];
		//	xr.x -= static_cast<int>((xr.x / (SIDE_LENGTH)) + 0.5) * (SIDE_LENGTH);
		//	xr.y -= static_cast<int>((xr.y / (SIDE_LENGTH)) + 0.5) * (SIDE_LENGTH);
		//	xr.z -= static_cast<int>((xr.z / (SIDE_LENGTH)) + 0.5) * (SIDE_LENGTH);
			xr.x -= round(xr.x / SIDE_LENGTH) * SIDE_LENGTH;
			xr.y -= round(xr.y / SIDE_LENGTH) * SIDE_LENGTH;
			xr.z -= round(xr.z / SIDE_LENGTH) * SIDE_LENGTH;
			double r = sqrt((xr.x*xr.x)+(xr.y*xr.y)+(xr.z*xr.z));
			if (r < (SIDE_LENGTH)/2) {
				int ig = static_cast<int>(r/delg);
				g[ig] +=  2;
			}
		}
	}
	if (i == 19) {
	for (int l = 0; l < bucketn; l++) {
		double r = delg*(l+0.5);
		double vb = ((l+1)*(l+1)*(l+1) - (l*l*l))*(delg*delg*delg);
		double nid = (4.0/3.0)*PI*vb*d;
		g[l] = (double) g[l] / (ngr*NPART*nid);
	}

	for (int i = 0; i < bucketn; i++) {
		cout << g[i] << endl;
	}
//	double max = 0;
//	for (int i = 0; i < bucketn; i++) {
//		if (g[i] > max)
//			max = g[i];
//	}
//	for (int i = 0; i < bucketn; i++)
//		g[i] /= max;
	Gnuplot r("point");
	vector<double> x1,y1;
	for (int j = 0; j < bucketn; j++) {
		y1.push_back(g[j]);
		x1.push_back((j+0.5)*delg);
	}
	r.set_smooth();
	r.reset_plot();
	r.unset_grid();
	r.plot_xy(x1,y1, "Test");
	cin.get();
	}
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "Usage: ./MD <density>\n";
		return 0;
	}


	double d = stod(argv[1]);
	DENSITY_FACTOR = pow((1/d), 0.333333333);
	SIDE_LENGTH = 6*DENSITY_FACTOR;
	X_SIZE = 1/SIDE_LENGTH;
	ECUT = 4*((1/pow((SIDE_LENGTH/2), 12.0)) - (1/pow(SIDE_LENGTH/2, 6.0)));



	if (argc == 3 && strncmp(argv[2], "-g", 2) == 0) {
		// Setting up the initial condition for the simulation including
		// the temperature, density, and initial velocities.
		init();
		vector<double> x1, y1, z1;
		Gnuplot g("points");
		for (int i = 0; i < NPART; i++) {
			x1.push_back(x[i].x);
			y1.push_back(x[i].y);
			z1.push_back(x[i].z);
		}
		int index = 0;
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
			for (int i = 0; i < NPART; i++) {
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
	}
	else {
		for (int i = 0; i < 20; i++) {
			for (int i = 0; i < NPART; i++) {
				f[i].x = 0;
				f[i].y = 0;
				f[i].z = 0;
				x[i].x = 0;
				x[i].y = 0;
				x[i].z = 0;
				velocities[i].x = 0;
				velocities[i].y = 0;
				velocities[i].z = 0;
			}
			tmp = 2;
		// Setting up the initial condition for the simulation including
		// the temperature, density, and initial velocities.
		init();
		for (double t = 0; t < 1e0; t+=DT) {
			force();
			integrate(t);
			//cout << etot << endl;
			//cout << tmp << endl;
		}
		radDistro(100, d, i);
		}
	}
	cout << tmp << endl << endl;

	// Create the radial distribution function and then plot it.
	//radDistro(50, d, 0);
	cin.get();

}
