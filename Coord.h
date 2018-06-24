class Coord {
	public:
		Coord();
		Coord(int i, int j, int k);
		Coord(double i, double j, double k);
		Coord operator+(const Coord& c) const;
		Coord operator+(const double& i) const;
		Coord operator/(const Coord& c) const;
		Coord operator/(const double& i) const;
		Coord operator*(const Coord& c) const;
		Coord operator*(const double& i) const;
		Coord operator-(const Coord& c) const;
		Coord operator-(const double& i) const;
		Coord operator=(const Coord& c);
		double x;
		double y;
		double z;
};
