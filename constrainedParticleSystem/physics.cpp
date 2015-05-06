/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/


#include "chain.h"

#include "physics.h"
#include <cassert>
#include <cmath>

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
extern Vector inputf;
extern int cring;
extern int controlMode;

double time = 0;

int print_matrix(FILE *f, const gsl_matrix *m)
{
	int status, n = 0;

	for (size_t i = 0; i < m->size1; i++) {
		for (size_t j = 0; j < m->size2; j++) {
			if ((status = fprintf(f, "%g,", gsl_matrix_get(m, i, j))) < 0)
				return -1;
			n += status;
		}

		if ((status = fprintf(f, "\n")) < 0)
			return -1;
		n += status;
	}

	return n;
}


void computeC1(Chain &chain, gsl_matrix *m){

	Vector *p = chain.p;
	int i;

	//second point
	gsl_matrix_set(m, 0, 0, 2. * p[1].x);
	gsl_matrix_set(m, 0, 1, 2. * p[1].y);

	//rigid edge
	for (i = 0; i < chain.n - 1; i++){
		double x0, x1, y0, y1;
		x0 = 2 * p[i].x - 2 * p[i + 1].x;
		y0 = 2 * p[i].y - 2 * p[i + 1].y;
		x1 = -x0;
		y1 = -y0;

		gsl_matrix_set(m, i, i * 2, x0); //set x[i-1]
		gsl_matrix_set(m, i, i * 2 +1, y0); //set y[i-1]
		gsl_matrix_set(m, i, i * 2 + 2, x1); //set x[i]
		gsl_matrix_set(m, i, i * 2 + 3, y1); //set y[i]
	}

	gsl_matrix_set(m, i++, 0, 1); //pin the first particle on origin x0 = 0 y0 = 0
	gsl_matrix_set(m, i++, 1, 1);

	//attach the last on ring
	if (cring){
		gsl_matrix_set(m, i, m->size2 - 2, 2. * p[chain.n - 1].x); //pin the first particle on origin
		gsl_matrix_set(m, i, m->size2 - 1, 2. * p[chain.n - 1].y + 1 );
	}

}

void computeC2(Chain &chain, gsl_matrix *m){

	Vector *v = chain.v;
	int i;
	
	//rigid connector

	for (i = 0; i < chain.n - 1; i++){
		double xy[4];
		xy[0] = 2 * v[i].x - 2 * v[i + 1].x;
		xy[1] = 2 * v[i].y - 2 * v[i + 1].y;
		xy[2] = -xy[0];
		xy[3] = -xy[1];

		gsl_matrix_set(m, i, i * 2, xy[0]); //set x[i-1]
		gsl_matrix_set(m, i, i * 2 + 1, xy[1]); //set y[i-1]
		gsl_matrix_set(m, i, i * 2 + 2, xy[2]); //set x[i]
		gsl_matrix_set(m, i, i * 2 + 3, xy[3]); //set y[i]
	}
	i += 2;
	//attach the last on ring
	if (cring){
		gsl_matrix_set(m, i, m->size2 - 2, 2. * v[chain.n - 1].x); //pin the first particle on origin
		gsl_matrix_set(m, i, m->size2 - 1, 2. * v[chain.n - 1].y);
	}

}

void computeC(Chain &chain, gsl_vector *c){
	Vector *p = chain.p;
	int i;
	//rigid connector
	for (i = 0; i < chain.n - 1; i++){
		double len = length(sub(p[i], p[i + 1]));
		//gsl_vector_set(c, i, pow(p[i].x - p[i + 1].x, 2) + pow(p[i].y - p[i + 1].y, 2)- chain.l*chain.l);
		gsl_vector_set(c, i, len * len -  chain.l*chain.l);
	}

	i += 2; //pin the first particle on origin 

	//attach the last on ring
	if (cring){
		gsl_vector_set(c, i, pow(p[chain.n - 1].x, 2) + pow(p[chain.n - 1].y + .5, 2) - .25); //pin the first particle on origin
	}

}

void computeAcceleration(Chain &chain, Vector a[])
{
	/* for you to implement ... */

	size_t c_num, n; //constrain equation num, particle num
	gsl_matrix *A = NULL, *C1 = NULL, *C2 = NULL;
	gsl_matrix_view AM, AC1, AC1T;
	gsl_vector *b = NULL, *x = NULL;
	gsl_vector_view f, pC2, ac;

	n = chain.n;
	c_num = n + 2;

	//AX = b
	A = gsl_matrix_calloc(2 * n + c_num, 2 * n + c_num);//set to zeros
	AM = gsl_matrix_submatrix(A, 0, 0, n * 2, n * 2);

	AC1 = gsl_matrix_submatrix(A, n * 2, 0, c_num, n * 2);
	AC1T = gsl_matrix_submatrix(A, 0, n * 2, n * 2, c_num);

	//right part of equation
	b = gsl_vector_calloc(2*n + c_num);
	f = gsl_vector_subvector(b, 0,n * 2);
	pC2 = gsl_vector_subvector(b, n * 2, c_num);


	C1 = gsl_matrix_calloc(c_num , n * 2);//set to zeros
	C2 = gsl_matrix_calloc(c_num, n * 2);//set to zeros

	x = gsl_vector_alloc(c_num + n * 2);

	/*get C1 C2*/
	/*not implement*/
	
	computeC1(chain,C1);
	computeC2(chain, C2);
	

	//printf("C1:\n");
	//print_matrix(stdout, C1);
	//printf("C2:\n");
	//print_matrix(stdout, C2);



	//set M
	gsl_matrix_set_identity(&AM.matrix);
	gsl_matrix_scale(&AM.matrix, chain.m);

	//set dC/dq
	gsl_matrix_memcpy(&AC1.matrix, C1);
	//set transpos dC/dq
	gsl_matrix_transpose_memcpy(&AC1T.matrix, C1);

	//set f
	assert(n * 2 == f.vector.size);
	for (int i = 0; i < n; i++){
		//damping -a*v
		gsl_vector_set(&f.vector, 2 * i,(i == n - 1 ||  controlMode ? 2 * i, inputf.x : 0) - chain.a * chain.v[i].x);
		gsl_vector_set(&f.vector, 2 * i + 1, (i == n - 1 || controlMode ? inputf.y : 0) + -1 * chain.m - chain.a * chain.v[i].y);
	}

	//set -dC'/dq * q' (not stable)
	gsl_vector *bv = gsl_vector_alloc(n * 2);
	for (int i = 0; i < n; i++){
		gsl_vector_set(bv, 2 * i, chain.v[i].x);
		gsl_vector_set(bv, 2 * i + 1, chain.v[i].y);
	}
	gsl_blas_dgemv(CblasNoTrans, -1., C2, bv, 0, &pC2.vector);

	//stablized
	
	gsl_vector *C = gsl_vector_calloc(c_num);//clear to zero
	computeC(chain, C);


	//double constraintsError = 0;
	//for (int i = 0; i<c_num; i++) {
	//	constraintsError += gsl_vector_get(C, i);

	//}

	//printf("%f\t%f\n", time, constraintsError);
	//time += chain.dt;

	//gsl_vector_fprintf(stdout, C, "%g");
	//printf("\n");

	gsl_blas_dgemv(CblasNoTrans, -2.*chain.b, C1, bv, 1, &pC2.vector);
	gsl_vector_scale(C, -chain.b * chain.b);
	gsl_vector_add(&pC2.vector, C);
	gsl_vector_free(C);

	gsl_vector_free(bv);

	//print_matrix(stdout, A);
	//gsl_vector_fprintf(stdout, b, "%g");

	gsl_matrix *V = gsl_matrix_alloc(A->size1,A->size1);
	gsl_vector *S = gsl_vector_alloc(A->size2), *work = gsl_vector_alloc(A->size2);

	//print_matrix(stdout, A);
	gsl_linalg_SV_decomp(A, V,S,work);//A is replaced by U
	//check singularity

	double eps = 1e-6;
	double maxs = 0.;

	//get the largest singular value
	for (int i = 0; i < S->size; i++){
		maxs = max(gsl_vector_get(S, i), maxs);
	}

	for (int i = 0; i<S->size; i++) {
		if (gsl_vector_get(S, i) < eps * maxs) {
			gsl_vector_set(S, i, 0.0);
		}
	}

	gsl_linalg_SV_solve(A,V,S,b,x);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);

	ac = gsl_vector_subvector(x, 0, n * 2);



	for (int i = 0; i < ac.vector.size; i++){
		if (fabs(gsl_vector_get(&ac.vector, i)) < eps)
			gsl_vector_set(&ac.vector, i, 0);
	}
	for (int i = 0; i < chain.n; i++){
		set(a[i], gsl_vector_get(&ac.vector, 2 * i), gsl_vector_get(&ac.vector, 2 * i + 1), 0.);
	}

	gsl_matrix_free(C2);
	gsl_matrix_free(C1);
	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_vector_free(x);
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(Chain &chain)
{
  Vector *a;

  a = new Vector[chain.n];
  
  computeAcceleration(chain, a);

  for (int i = 0; i < chain.n; i++){
	//printf("a%d:%f, %f\n", i,a[i].x,a[i].y);
	  chain.p[i].x += chain.dt * chain.v[i].x;
	  chain.p[i].y += chain.dt * chain.v[i].y;
	  chain.v[i].x += chain.dt * a[i].x;
	  chain.v[i].y += chain.dt * a[i].y;
  }
 // printf("\n");
  delete[] a;
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */

void RK4(Chain &chain)
{

  Vector *a, *F1p, *F2p, *F3p, *F4p, *F1v, *F2v, *F3v, *F4v;
  int i;

  Chain buffer  = chain;


  a = new Vector[chain.n];
  F1p = new Vector[chain.n];
  F2p = new Vector[chain.n];
  F3p = new Vector[chain.n];
  F4p = new Vector[chain.n];

  F1v = new Vector[chain.n];
  F2v = new Vector[chain.n];
  F3v = new Vector[chain.n];
  F4v = new Vector[chain.n];

  computeAcceleration(chain, a);
  

  for ( i = 0; i < chain.n; i++){
	  F1p[i] = scal(chain.v[i], chain.dt);
	  F1v[i] = scal(a[i], chain.dt);
	  buffer.p[i] = scal(F1p[i], .5);
	  buffer.v[i] = scal(F1v[i], .5);
	  buffer.p[i] = summ(chain.p[i] , buffer.p[i]);
	  buffer.v[i] = summ(chain.v[i], buffer.v[i]);
  }

  computeAcceleration(buffer, a);

  for ( i = 0; i < chain.n; i++){
	  F2p[i] = scal(buffer.v[i], chain.dt);
	  F2v[i] = scal(a[i], chain.dt);
	  buffer.p[i] = scal(F2p[i], .5);
	  buffer.v[i] = scal(F2v[i], .5);
	  buffer.p[i] = summ(chain.p[i], buffer.p[i]);
	  buffer.v[i] = summ(chain.v[i], buffer.v[i]);
  }


  computeAcceleration(buffer, a);

  for ( i = 0; i < chain.n; i++){
	  F3p[i] = scal(buffer.v[i], chain.dt);
	  F3v[i] = scal(a[i], chain.dt);
	  buffer.p[i] = scal(F3p[i], .5);
	  buffer.v[i] = scal(F3v[i], .5);
	  buffer.p[i] = summ(chain.p[i], buffer.p[i]);
	  buffer.v[i] = summ(chain.v[i], buffer.v[i]);
  }
         
  computeAcceleration(buffer, a);



  for ( i = 0; i < chain.n; i++){
	  F4p[i] = scal(chain.v[i], chain.dt);
	  F4v[i] = scal(a[i], chain.dt);
	  buffer.p[i] = scal(F2p[i], 2);
	  buffer.v[i] = scal(F3p[i], 2);

	  buffer.p[i] = summ(buffer.v[i], buffer.p[i]);
	  buffer.p[i] = summ(buffer.p[i],F1p[i]);
	  buffer.p[i] = summ(buffer.p[i], F4p[i]);
	  buffer.p[i] = scal(buffer.p[i], 1. / 6);
	  chain.p[i] = summ(buffer.p[i] , chain.p[i]);


	  buffer.p[i] = scal(F2v[i], 2.);
	  buffer.v[i] = scal(F3v[i], 2.);
	  buffer.p[i] = summ(buffer.p[i], buffer.v[i]);
	  buffer.p[i] = summ(buffer.p[i], F1v[i]);
	  buffer.p[i] = summ(buffer.p[i], F4v[i]);
	  buffer.p[i] = scal(buffer.p[i], 1. / 6);
	  chain.v[i] = summ(buffer.p[i], chain.v[i]);

	  

  }

  delete [] a;
  delete[] F1p;
  delete[] F2p;
  delete[] F3p;
  delete[] F4p;
  delete[] F1v;
  delete[] F2v;
  delete[] F3v;
  delete[] F4v;

  return;  
}