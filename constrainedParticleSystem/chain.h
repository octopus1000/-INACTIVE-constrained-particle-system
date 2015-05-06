#include "openGL-headers.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <GL\glut.h>
#include <cmath>


#ifndef _CHAIN_H_
#define _CHAIN_H_

#define MAX_PARTICLE_NUM 101
#define pi 3.141592653589793238462643383279 

struct Vector{
	double x, y, z;
};

typedef Vector Particle;
typedef Vector Point;

struct Chain{
	Point p[MAX_PARTICLE_NUM];
	Vector v[MAX_PARTICLE_NUM];
	int n;//number of particles
	double m;//mass for single particle
	double l;//length for rigid edge between 2 particles

	double dt;//timestamp
	double b;//damping parameter
	double a;
	char *integrator;
};



void showParticles(Chain &chain);
void initChain(Chain &chain);
void vertPos(Chain &);
void horiPos(Chain &);
void state(Chain &chain);

inline void set(Vector &p, double x, double y, double z);
inline Vector sub(Vector &a, Vector &b);
inline double dotProduct(Vector &a, Vector &b);
inline double length(Vector &v);

inline Vector scal(Vector &a, double s){
	Vector c;
	set(c, a.x *s, a.y * s, a.z*s);
	return c;
}
inline Vector summ(Vector &a, Vector &b){
	Vector c;
	set(c, a.x + b.x, a.y + b.y, a.z + b.z);
	return c;
}

#endif

//for Particle