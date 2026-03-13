#ifndef FLOAT256_H
#define FLOAT256_H

// float256: quad-double arithmetic (4 non-overlapping doubles).
// Value represented = x[0] + x[1] + x[2] + x[3], where |x[i+1]| <= ulp(x[i])/2.
// Gives ~212 bits of mantissa — enough for exact predicates on double coords.

class float256
{
public:
	// data members: x[0] is most significant, x[3] is least
	double x[4];

	// constructors
	float256();
	float256(double a);
	float256(int a);
	float256(double x0, double x1, double x2, double x3);

	// conversion operators
	explicit operator double() const;
	explicit operator int() const;

	// member functions
	float256 Sqrt() const;

	// operator redefinition

	// Addition
	float256 operator+(const float256& z) const;
	float256 operator+(double z) const;
	float256 operator+(int z) const;
	friend float256 operator+(int x, const float256& y);
	friend float256 operator+(double x, const float256& y);
	float256& operator+=(const float256& z);
	float256& operator+=(double z);
	float256& operator+=(int z);

	// Subtraction
	float256 operator-() const;
	float256 operator-(const float256& z) const;
	float256 operator-(double z) const;
	float256 operator-(int z) const;
	friend float256 operator-(int x, const float256& y);
	friend float256 operator-(double x, const float256& y);
	float256& operator-=(const float256& z);
	float256& operator-=(double z);
	float256& operator-=(int z);

	// Multiplication
	float256 operator*(const float256& z) const;
	float256 operator*(double z) const;
	float256 operator*(int z) const;
	friend float256 operator*(int x, const float256& y);
	friend float256 operator*(double x, const float256& y);
	float256& operator*=(const float256& z);
	float256& operator*=(double z);
	float256& operator*=(int z);

	// Division
	float256 operator/(const float256& z) const;
	float256 operator/(double z) const;
	float256 operator/(int z) const;
	friend float256 operator/(int x, const float256& y);
	friend float256 operator/(double x, const float256& y);
	float256& operator/=(const float256& z);
	float256& operator/=(double z);
	float256& operator/=(int z);

	// Assignment
	float256& operator=(double z);
	float256& operator=(int z);

	// Comparison
	bool operator==(const float256& z) const;
	bool operator==(double z) const;
	bool operator==(int z) const;

	bool operator!=(const float256& z) const;
	bool operator!=(double z) const;
	bool operator!=(int z) const;

	bool operator<(const float256& z) const;
	bool operator<(double z) const;
	bool operator<(int z) const;

	bool operator<=(const float256& z) const;
	bool operator<=(double z) const;
	bool operator<=(int z) const;

	bool operator>(const float256& z) const;
	bool operator>(double z) const;
	bool operator>(int z) const;

	bool operator>=(const float256& z) const;
	bool operator>=(double z) const;
	bool operator>=(int z) const;

	// Error-free arithmetic primitives (public so free functions can use them)
	static void two_sum(double a, double b, double& s, double& e);
	static void two_prod(double a, double b, double& p, double& e);

private:
	static void split(double a, double& hi, double& lo);

	// Renormalize a 4-component expansion so components are non-overlapping
	void renormalize();

	// Internal: add a double to this expansion
	static float256 add_d(const float256& a, double b);
	// Internal: multiply expansion by a double
	static float256 mul_d(const float256& a, double b);
};

#endif // FLOAT256_H
