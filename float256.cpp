#include "float256.h"
#include <cmath>

// ============================================================
// Error-free arithmetic primitives
// ============================================================

// Two-Sum (Knuth): s + e = a + b exactly. Works for any a, b.
inline void float256::two_sum(double a, double b, double& s, double& e) {
	s = a + b;
	double v = s - a;
	e = (a - (s - v)) + (b - v);
}

// Quick-Two-Sum: s + e = a + b exactly. REQUIRES |a| >= |b|.
static inline void quick_two_sum(double a, double b, double& s, double& e) {
	s = a + b;
	e = b - (s - a);
}

// Dekker splitting: a = hi + lo, each with 26 significant bits.
inline void float256::split(double a, double& hi, double& lo) {
	static const double SPLIT = 134217729.0; // 2^27 + 1
	double t = SPLIT * a;
	hi = t - (t - a);
	lo = a - hi;
}

// Two-Product (Dekker): p + e = a * b exactly.
inline void float256::two_prod(double a, double b, double& p, double& e) {
#ifdef FP_FAST_FMA
	p = a * b;
	e = fma(a, b, -p);
#else
	double a_hi, a_lo, b_hi, b_lo;
	p = a * b;
	split(a, a_hi, a_lo);
	split(b, b_hi, b_lo);
	e = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
#endif
}

// Three-sum: replaces (a,b,c) with non-overlapping (a,b,c) where a+b+c is preserved.
static inline void three_sum(double& a, double& b, double& c) {
	double t1, t2, t3;
	float256::two_sum(a, b, t1, t2);
	float256::two_sum(c, t1, a, t3);
	float256::two_sum(t2, t3, b, c);
}

// Three-sum2: like three_sum but only 2 outputs (c is just added to b).
static inline void three_sum2(double& a, double& b, double c) {
	double t1, t2, t3;
	float256::two_sum(a, b, t1, t2);
	float256::two_sum(c, t1, a, t3);
	b = t2 + t3;
}

// ============================================================
// Renormalize (4 components)
// Based on the QD library by Hida/Li/Bailey.
// Ensures components are non-overlapping with decreasing magnitude.
// ============================================================
void float256::renormalize() {
	double s0, s1, s2 = 0.0, s3 = 0.0;
	double e;

	// First pass: bottom-up cascade using quick_two_sum
	// After this pass, x[0] absorbs all magnitude.
	quick_two_sum(x[2], x[3], s0, x[3]);
	quick_two_sum(x[1], s0, s0, x[2]);
	quick_two_sum(x[0], s0, x[0], x[1]);

	s0 = x[0];
	s1 = x[1];

	// Second pass: top-down to ensure non-overlapping
	quick_two_sum(x[0], x[1], s0, s1);
	if (s1 != 0.0) {
		quick_two_sum(s1, x[2], s1, s2);
		if (s2 != 0.0)
			quick_two_sum(s2, x[3], s2, s3);
		else
			quick_two_sum(s1, x[3], s1, s2);
	} else {
		quick_two_sum(s0, x[2], s0, s1);
		if (s1 != 0.0)
			quick_two_sum(s1, x[3], s1, s2);
		else
			quick_two_sum(s0, x[3], s0, s1);
	}
	x[0] = s0;
	x[1] = s1;
	x[2] = s2;
	x[3] = s3;
}

// Renormalize 5 components into 4.
static void renorm5(double& c0, double& c1, double& c2, double& c3, double c4) {
	double s0, s1, s2 = 0.0, s3 = 0.0;
	double e;

	// Bottom-up cascade
	quick_two_sum(c3, c4, c3, c4);
	quick_two_sum(c2, c3, c2, c3);
	quick_two_sum(c1, c2, c1, c2);
	quick_two_sum(c0, c1, c0, c1);

	s0 = c0;
	s1 = c1;

	// Top-down pass
	quick_two_sum(c0, c1, s0, s1);
	if (s1 != 0.0) {
		quick_two_sum(s1, c2, s1, s2);
		if (s2 != 0.0)
			quick_two_sum(s2, c3, s2, s3);
		else
			quick_two_sum(s1, c3, s1, s2);
	} else {
		quick_two_sum(s0, c2, s0, s1);
		if (s1 != 0.0)
			quick_two_sum(s1, c3, s1, s2);
		else
			quick_two_sum(s0, c3, s0, s1);
	}
	c0 = s0;
	c1 = s1;
	c2 = s2;
	c3 = s3;
}

// ============================================================
// Constructors
// ============================================================

float256::float256() {
	x[0] = x[1] = x[2] = x[3] = 0.0;
}

float256::float256(double a) {
	x[0] = a;
	x[1] = x[2] = x[3] = 0.0;
}

float256::float256(int a) {
	x[0] = static_cast<double>(a);
	x[1] = x[2] = x[3] = 0.0;
}

float256::float256(double x0, double x1, double x2, double x3) {
	x[0] = x0; x[1] = x1; x[2] = x2; x[3] = x3;
}

// ============================================================
// Conversion operators
// ============================================================

float256::operator double() const {
	// Sum from least to most significant for best rounding
	return ((x[3] + x[2]) + x[1]) + x[0];
}

float256::operator int() const {
	return static_cast<int>(static_cast<double>(*this));
}

// ============================================================
// Addition: QD + QD (sloppy add from QD library)
// ============================================================

float256 float256::operator+(const float256& b) const {
	double s0, s1, s2, s3;
	double t0, t1, t2, t3;

	two_sum(x[0], b.x[0], s0, t0);
	two_sum(x[1], b.x[1], s1, t1);
	two_sum(x[2], b.x[2], s2, t2);
	two_sum(x[3], b.x[3], s3, t3);

	two_sum(s1, t0, s1, t0);
	three_sum(s2, t0, t1);
	three_sum2(s3, t0, t2);
	t0 = t0 + t1 + t3;

	renorm5(s0, s1, s2, s3, t0);
	return float256(s0, s1, s2, s3);
}

// QD + double
float256 float256::add_d(const float256& a, double b) {
	double s0, s1, s2, s3;
	double e;

	two_sum(a.x[0], b, s0, e);
	two_sum(a.x[1], e, s1, e);
	two_sum(a.x[2], e, s2, e);
	s3 = a.x[3] + e;

	renorm5(s0, s1, s2, s3, 0.0);
	return float256(s0, s1, s2, s3);
}

float256 float256::operator+(double z) const {
	return add_d(*this, z);
}

float256 float256::operator+(int z) const {
	return add_d(*this, static_cast<double>(z));
}

float256 operator+(int x, const float256& y) {
	return float256::add_d(y, static_cast<double>(x));
}

float256 operator+(double x, const float256& y) {
	return float256::add_d(y, x);
}

float256& float256::operator+=(const float256& z) { *this = *this + z; return *this; }
float256& float256::operator+=(double z) { *this = add_d(*this, z); return *this; }
float256& float256::operator+=(int z) { *this = add_d(*this, (double)z); return *this; }

// ============================================================
// Negation and Subtraction
// ============================================================

float256 float256::operator-() const {
	return float256(-x[0], -x[1], -x[2], -x[3]);
}

float256 float256::operator-(const float256& z) const { return *this + (-z); }
float256 float256::operator-(double z) const { return add_d(*this, -z); }
float256 float256::operator-(int z) const { return add_d(*this, (double)(-z)); }

float256 operator-(int x, const float256& y) { return float256((double)x) + (-y); }
float256 operator-(double x, const float256& y) { return float256(x) + (-y); }

float256& float256::operator-=(const float256& z) { *this = *this - z; return *this; }
float256& float256::operator-=(double z) { *this = add_d(*this, -z); return *this; }
float256& float256::operator-=(int z) { *this = add_d(*this, (double)(-z)); return *this; }

// ============================================================
// Multiplication: QD * double
// ============================================================

float256 float256::mul_d(const float256& a, double b) {
	double p0, p1, p2, p3;
	double e0, e1, e2, e3;

	two_prod(a.x[0], b, p0, e0);
	two_prod(a.x[1], b, p1, e1);
	two_prod(a.x[2], b, p2, e2);
	p3 = a.x[3] * b;

	// Cascade errors downward through the components
	two_sum(p1, e0, p1, e0);
	two_sum(p2, e0, p2, e0);
	two_sum(p2, e1, p2, e1);
	p3 += e0 + e1 + e2;

	renorm5(p0, p1, p2, p3, 0.0);
	return float256(p0, p1, p2, p3);
}

// ============================================================
// Multiplication: QD * QD (sloppy mul from QD library)
// ============================================================

float256 float256::operator*(const float256& b) const {
	double p0, p1, p2, p3, p4, p5;
	double q0, q1, q2, q3, q4, q5;
	double s0, s1, s2;
	double t0, t1;

	two_prod(x[0], b.x[0], p0, q0);
	two_prod(x[0], b.x[1], p1, q1);
	two_prod(x[1], b.x[0], p2, q2);
	two_prod(x[0], b.x[2], p3, q3);
	two_prod(x[1], b.x[1], p4, q4);
	two_prod(x[2], b.x[0], p5, q5);

	// Start accumulation
	three_sum(p1, p2, q0);

	// Six-three sum of p2, q1, q2, p3, p4, p5
	three_sum(p2, q1, q2);
	three_sum(p3, p4, p5);

	// (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5)
	two_sum(p2, p3, s0, t0);
	two_sum(q1, p4, s1, t1);
	s2 = q2 + p5;
	two_sum(s1, t0, s1, t0);
	s2 += (t0 + t1);

	// O(eps^3) terms: just accumulate without error tracking
	s1 += x[0] * b.x[3] + x[1] * b.x[2] + x[2] * b.x[1] + x[3] * b.x[0]
	    + q0 + q3 + q4 + q5;

	renorm5(p0, p1, s0, s1, s2);
	return float256(p0, p1, s0, s1);
}

float256 float256::operator*(double z) const { return mul_d(*this, z); }
float256 float256::operator*(int z) const { return mul_d(*this, (double)z); }

float256 operator*(int x, const float256& y) { return float256::mul_d(y, (double)x); }
float256 operator*(double x, const float256& y) { return float256::mul_d(y, x); }

float256& float256::operator*=(const float256& z) { *this = *this * z; return *this; }
float256& float256::operator*=(double z) { *this = mul_d(*this, z); return *this; }
float256& float256::operator*=(int z) { *this = mul_d(*this, (double)z); return *this; }

// ============================================================
// Division: iterative refinement (from QD library)
// ============================================================

float256 float256::operator/(const float256& b) const {
	double q0, q1, q2, q3;
	float256 r;

	q0 = x[0] / b.x[0];
	r = *this - mul_d(b, q0);

	q1 = r.x[0] / b.x[0];
	r -= mul_d(b, q1);

	q2 = r.x[0] / b.x[0];
	r -= mul_d(b, q2);

	q3 = r.x[0] / b.x[0];

	renorm5(q0, q1, q2, q3, 0.0);
	return float256(q0, q1, q2, q3);
}

float256 float256::operator/(double z) const { return *this / float256(z); }
float256 float256::operator/(int z) const { return *this / float256((double)z); }

float256 operator/(int x, const float256& y) { return float256((double)x) / y; }
float256 operator/(double x, const float256& y) { return float256(x) / y; }

float256& float256::operator/=(const float256& z) { *this = *this / z; return *this; }
float256& float256::operator/=(double z) { *this = *this / float256(z); return *this; }
float256& float256::operator/=(int z) { *this = *this / float256((double)z); return *this; }

// ============================================================
// Assignment
// ============================================================

float256& float256::operator=(double z) {
	x[0] = z; x[1] = x[2] = x[3] = 0.0;
	return *this;
}

float256& float256::operator=(int z) {
	x[0] = (double)z; x[1] = x[2] = x[3] = 0.0;
	return *this;
}

// ============================================================
// Comparison
// ============================================================

bool float256::operator==(const float256& z) const {
	return x[0] == z.x[0] && x[1] == z.x[1] && x[2] == z.x[2] && x[3] == z.x[3];
}
bool float256::operator==(double z) const {
	return x[0] == z && x[1] == 0.0 && x[2] == 0.0 && x[3] == 0.0;
}
bool float256::operator==(int z) const { return *this == (double)z; }

bool float256::operator!=(const float256& z) const { return !(*this == z); }
bool float256::operator!=(double z) const { return !(*this == z); }
bool float256::operator!=(int z) const { return !(*this == z); }

bool float256::operator<(const float256& z) const {
	for (int i = 0; i < 4; i++) {
		if (x[i] < z.x[i]) return true;
		if (x[i] > z.x[i]) return false;
	}
	return false;
}
bool float256::operator<(double z) const { return *this < float256(z); }
bool float256::operator<(int z) const { return *this < float256(z); }

bool float256::operator<=(const float256& z) const { return !(z < *this); }
bool float256::operator<=(double z) const { return *this <= float256(z); }
bool float256::operator<=(int z) const { return *this <= float256(z); }

bool float256::operator>(const float256& z) const { return z < *this; }
bool float256::operator>(double z) const { return *this > float256(z); }
bool float256::operator>(int z) const { return *this > float256(z); }

bool float256::operator>=(const float256& z) const { return !(*this < z); }
bool float256::operator>=(double z) const { return *this >= float256(z); }
bool float256::operator>=(int z) const { return *this >= float256(z); }

// ============================================================
// Sqrt: Newton-Raphson
// ============================================================

float256 float256::Sqrt() const {
	if (x[0] == 0.0) return float256(0.0);
	if (x[0] < 0.0) return float256(0.0);

	float256 r(std::sqrt(x[0]));
	// 3 Newton iterations: r = (r + this/r) * 0.5
	for (int i = 0; i < 3; i++) {
		r = (r + *this / r) * 0.5;
	}
	return r;
}
