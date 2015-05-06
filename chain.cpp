#include "chain.h"
#include <assert.h> 

//compute function

extern Vector inputf;
extern int controlMode;
extern int cring;

inline void set(Vector &p, double x, double y, double z){
	p.x = x;
	p.y = y;
	p.z = z;
}
inline Vector sub(Vector &a, Vector &b){
	Vector c;
	set(c, a.x - b.x, a.y - b.y, a.z - b.z);
	return c;
}
//inline Vector summ(Vector &a, Vector &b){
//	Vector c;
//	set(c, a.x + b.x, a.y + b.y, a.z + b.z);
//	return c;
//}
//inline Vector scal(Vector &a, double &s){
//	Vector c;
//	set(c, a.x *s, a.y * s, a.z*s);
//	return c;
//}

inline double dotProduct(Vector &a, Vector &b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline double length(Vector &v){
	return sqrt(dotProduct(v, v));
}

void showParticles(Chain &chain){
	GLUquadric *quad = gluNewQuadric();
	double r = chain.l / 5.;

	gluQuadricDrawStyle(quad, GLU_FILL);
	gluQuadricNormals(quad, GLU_SMOOTH);
	//draw spheres

	//draw Ring
	glPushMatrix();
	glTranslated(0, -0.5,0);
	glutWireTorus(0.01, 0.5, 10, 20);
	glPopMatrix();

	for (int i = 0; i < chain.n; i++){
		glPushMatrix();
		glTranslated(chain.p[i].x, chain.p[i].y,chain.p[i].z);


		gluSphere(quad, /*i == 0 || i == chain.n - 1 ? r * 2 :*/ r, 10, 10);
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,'p');
		
			//draw cylinder
		if (i < chain.n - 1){
			Vector v = sub(chain.p[i + 1], chain.p[i]);
			v.z = v.z ? v.z : 0.0001;
			double len = length(v);
			double ax = 57.2957795*acos(v.z / len);
			ax = v.z < 0. ? -ax : ax;
			
			glRotated(ax, -v.y * v.z, v.x * v.z, 0.);
			gluCylinder(quad, chain.l / 20., chain.l / 20., len, 10, 10);
		}
		glPopMatrix();
	}

	gluDeleteQuadric(quad);
}

void initChain(Chain &chain){

	Particle *p = chain.p;

	
	//set params
	assert(chain.n > 1);
	chain.m = 1. / chain.n;
	chain.l = (1. + 0.1/ chain.n) / (chain.n - 1);

	//set init positions for particles
	set(p[0], 0, 0, 0);
	set(chain.v[0], 0, 0, 0);

	//vertPos(chain);
	horiPos(chain);

	set(inputf, 0, 0, 0);
}

void horiPos(Chain &chain){
	Vector *p = chain.p;
	for (int i = 1; i < chain.n; i++){
		//line up toward y axis
		set(p[i], p[i - 1].x + chain.l, 0, 0);
		set(chain.v[i], 0, 0, 0);
	}
}

void vertPos(Chain &chain){
	Vector *p = chain.p;
	int i = 1;
	//for (int i = 1; i < chain.n; i++){
	//	//line up toward y axis
	//	//set(p[i], p[i - 1].x + chain.l,0, 0);
	//	set(p[i], 0, p[i - 1].y - chain.l, 0);
	//	set(chain.v[i], 0, 0, 0);
	//}

	while (i < ceil(chain.n / 2.)){
		set(p[i], 0, p[i - 1].y - chain.l, 0);
		set(chain.v[i], 0, 0, 0);
		i++;
	}

	while (i < chain.n){
		set(p[i], p[i - 1].x + chain.l, -0.5, 0);
		set(chain.v[i], 0, 0, 0);
		i++;
	}

}


void state(Chain &chain){
	printf("integrator:%s\n", chain.integrator);
	printf("dt:%f\n", chain.dt);
	printf("b:%f, a:%f\n", chain.b, chain.a);
	printf("particle number:%d\n", chain.n);
	printf("force on point:");
	if (controlMode){
		printf("all points\n");
	}
	else{
		printf("last point\n");
	}

	if (cring){
		printf("cring:last point on ring\n");
	}
	else{
		printf("cring:last point free\n");
	}

	printf("\n");

}