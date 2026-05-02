#ifndef EXTRACTBOUNDARYDATA_H
#define EXTRACTBOUNDARYDATA_H

class CartesianGridAndControlPoints;

void extractBoundaryDataForContinousFunction(CartesianGridAndControlPoints* G, double** u, double* bdry_u, double* bdry_un);

#endif