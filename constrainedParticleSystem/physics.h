/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include "chain.h"

void computeAcceleration(Chain &chain, Vector a[]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(Chain &chain);
void RK4(Chain &chain);

#endif

