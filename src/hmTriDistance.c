#include "hmTriDistance.h"
#include "hmVec3.h"
#include "hmUtility.h"
#include "hmConstants.h"
#include "hmVectorSizeT.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include<OpenCL/opencl.h>

void hmTriDistanceInitialize( hmTriDistance* distance )
{
   distance->surface         = NULL;

   hmDenseMatrixInitialize( &distance->isSource,      0, 0 );
   hmDenseMatrixInitialize( &distance->distance,      0, 0 );
   hmDenseMatrixInitialize( &distance->heatNeumann,   0, 0 );
   hmDenseMatrixInitialize( &distance->heatDirichlet, 0, 0 );
   hmDenseMatrixInitialize( &distance->potential,     0, 0 );

   hmSparseMatrixInitialize( &distance->laplacian,         0, 0, 0 );
   hmSparseMatrixInitialize( &distance->heatFlowNeumann,   0, 0, 0 );
   hmSparseMatrixInitialize( &distance->heatFlowDirichlet, 0, 0, 0 );

   hmCholeskyFactorInitialize( &distance->laplacianFactor );
   hmCholeskyFactorInitialize( &distance->heatFlowNeumannFactor );
   hmCholeskyFactorInitialize( &distance->heatFlowDirichletFactor );

   /* by default, use pure Neumann boundary conditions */
   distance->boundaryConditions = 0.;

   hmVectorDoubleInitialize( &distance->startTimesProcessor );
   hmVectorDoubleInitialize( &distance->startTimesWallClock );
}

void hmTriDistanceDestroy( hmTriDistance* distance )
{
   hmDenseMatrixDestroy( &distance->isSource );
   hmDenseMatrixDestroy( &distance->distance );
   hmDenseMatrixDestroy( &distance->heatNeumann );
   hmDenseMatrixDestroy( &distance->heatDirichlet );
   hmDenseMatrixDestroy( &distance->potential );

   hmSparseMatrixDestroy( &distance->laplacian );
   hmSparseMatrixDestroy( &distance->heatFlowNeumann );
   hmSparseMatrixDestroy( &distance->heatFlowDirichlet );

   hmCholeskyFactorDestroy( &distance->laplacianFactor );
   hmCholeskyFactorDestroy( &distance->heatFlowNeumannFactor );
   hmCholeskyFactorDestroy( &distance->heatFlowDirichletFactor );

   hmVectorDoubleDestroy( &distance->startTimesProcessor );
   hmVectorDoubleDestroy( &distance->startTimesWallClock );
}

void hmTriDistanceEstimateTime( hmTriDistance* distance )
{
   size_t nFaces = distance->surface->nFaces;
   size_t* facesBegin = distance->surface->faces;
   size_t* facesEnd = facesBegin + 3*nFaces;
   size_t* f;
   double* vertices = distance->surface->vertices;
   double *p0, *p1, *p2;
   hmVec3 e01, e12, e20;
   double meanEdgeLength = 0.;
   double nEdges = 0.;

   /* iterate over faces */
   for( f = facesBegin; f != facesEnd; f += 3 )
   {
      /* get vertex coordinates p0, p1, p2 */
      p0 = &vertices[ f[0]*3 ];
      p1 = &vertices[ f[1]*3 ];
      p2 = &vertices[ f[2]*3 ];

      /* add edge lengths to mean */
      hmVec3Sub( e01, p1, p0 );
      hmVec3Sub( e12, p2, p1 );
      hmVec3Sub( e20, p0, p2 );
      meanEdgeLength += hmVec3Norm( e01 );
      meanEdgeLength += hmVec3Norm( e12 );
      meanEdgeLength += hmVec3Norm( e20 );

      nEdges += 3.;
   }
   meanEdgeLength /= nEdges;

   /* set t to square of mean edge length */
   distance->time = hmSquare( meanEdgeLength );
}

void hmTriDistanceSetBoundaryConditions( hmTriDistance* distance,
                                         double boundaryConditions )
{
   distance->boundaryConditions = boundaryConditions;
}

void hmTriDistanceBuild( hmTriDistance* distance )
{
   size_t nVertices = distance->surface->nVertices;

   if( distance->surface == NULL )
   {
      fprintf( stderr, "Error: hmTriDistanceBuild -- must specify a surface!\n" );
      exit( 1 );
   }

   hmTriDistanceDestroy( distance );

   hmVectorDoubleInitialize( &distance->startTimesProcessor );
   hmVectorDoubleInitialize( &distance->startTimesWallClock );
   hmTriDistanceStartTiming( distance );

   hmDenseMatrixInitialize( &distance->isSource,  nVertices, 1 );
   hmDenseMatrixInitialize( &distance->distance,  nVertices, 1 );
   hmDenseMatrixInitialize( &distance->potential, nVertices, 1 );

   /* only allocate space for both solutions if necessary */
   if( distance->boundaryConditions < 1. ) /* partial Neumann */
   {
      hmDenseMatrixInitialize( &distance->heatNeumann,   nVertices, 1 );
   }
   if( distance->boundaryConditions > 0. ) /* partial Dirichlet */
   {
      hmDenseMatrixInitialize( &distance->heatDirichlet, nVertices, 1 );
   }

   hmTriDistanceBuildMatrices( distance );
   hmTriDistanceFactorMatrices( distance );

   hmTriDistanceStopTiming( distance, "Total build time" );
}

void hmTriDistanceUpdate( hmTriDistance* distance )
{
   hmTriDistanceStartTiming( distance );

   hmTriDistanceStartTiming( distance );
   hmTriDistanceSolveHeatEquation( distance );
   hmTriDistanceStopTiming( distance, "Solve heat equation" );

   hmTriDistanceStartTiming( distance );
   hmTriDistanceComputePotential( distance );
   hmTriDistanceStopTiming( distance, "Compute potential" );

   hmTriDistanceStartTiming( distance );
   hmTriDistanceSolvePoissonEquation( distance );
   hmTriDistanceStopTiming( distance, "Solve Poisson equation" );

   hmTriDistanceStopTiming( distance, "Total update time" );
}

void hmTriDistanceSolveHeatEquation( hmTriDistance* distance )
{
   size_t nVertices = distance->surface->nVertices;
   const double BC = distance->boundaryConditions;
   int err;
   cl_program program;
   cl_kernel kernel;
   cl_context context;
   cl_command_queue queue;
   cl_device_id device_id;
   int gpu=1;
   cl_mem input1;
   cl_mem input2;
   cl_mem output;
   size_t global;
   size_t local;
   const char *KernelSource = "\n" \
"__kernel HeatEquation( \n" \
" __global double*input1, \n" \
" __global double*input2, \n" \
" __global double*output, \n" \
" const double BC, \n" \
" const unsigned int count) \n" \
"{ \n" \
" int i = get_global_id(0); \n" \
" if(i < count){ \n" \
" output[i]=(1-BC)*input1[i]+BC*input2[i]; \n" \
"} \n" \
"} \n" \
"\n";
   err=clGetDeviceIDs(NULL, gpu?CL_DEVICE_TYPE_GPU:CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
   context=clCreateContext(0, 1, &device_id, NULL, NULL, &err);
   queue=clCreateCommandQueue(context, device_id, 0, &err);
   program=clCreateProgramWithSource(context, 1, (const char **)&KernelSource, NULL, &err);
   err=clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
   kernel=clCreateKernel(program, "HeatEquation", &err);

   /* only compute both solutions if necessary */
   if( BC < 1. ) /* partial Neumann */
   {
      hmCholeskyFactorBacksolve( &distance->heatFlowNeumannFactor,
                                 &distance->heatNeumann,
                                 &distance->isSource );
   }
   if( BC > 0. ) /* partial Dirichlet */
   {
      hmCholeskyFactorBacksolve( &distance->heatFlowDirichletFactor,
                                 &distance->heatDirichlet,
                                 &distance->isSource );
   }

   /* store the final solution in hmTriDistance::heat,
    * combining the two solutions if necessary */
   if( BC > 0. && BC < 1. )
   {
   input1=clCreateBuffer(context, CL_MEM_READ_ONLY, nVertices * sizeof(double), NULL, NULL);
   input2=clCreateBuffer(context, CL_MEM_READ_ONLY, nVertices * sizeof(double), NULL, NULL);
   output=clCreateBuffer(context, CL_MEM_WRITE_ONLY, nVertices * sizeof(double), NULL, NULL);
   err=clEnqueueWriteBuffer(queue, input1, CL_TRUE, 0, nVertices * sizeof(double), distance->heatNeumann.values, 0, NULL, NULL);
   err=clEnqueueWriteBuffer(queue, input2, CL_TRUE, 0, nVertices * sizeof(double), distance->heatDirichlet.values, 0, NULL, NULL);
   err=clSetKernelArg(kernel, 0, sizeof(cl_mem), &input1);
   err=clSetKernelArg(kernel, 1, sizeof(cl_mem), &input2);
   err=clSetKernelArg(kernel, 2, sizeof(cl_mem), &output);
   err=clSetKernelArg(kernel, 3, sizeof(double), &BC);
   err=clSetKernelArg(kernel, 4, sizeof(unsigned int), &nVertices );
   err=clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(int), &local, NULL);
   global=nVertices;
   err=clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
   clFinish(queue);
   err=clEnqueueReadBuffer(queue, output, CL_TRUE, 0, nVertices * sizeof(double), distance->heatNeumann.values, 0, NULL, NULL);
   clReleaseMemObject(input1);
   clReleaseMemObject(input2);
   clReleaseMemObject(output);
   clReleaseProgram(program);
   clReleaseKernel(kernel);
   clReleaseCommandQueue(queue);
   clReleaseContext(context);
   }
   else if( BC == 0. ) /* pure Neumann */
   {
      distance->heat = distance->heatNeumann.values;
   }
   else /* pure Dirichlet */
   {
      distance->heat = distance->heatDirichlet.values;
   }
}

void hmTriDistanceComputePotential( hmTriDistance* distance )
{
   /* array counters */
   size_t i;

   /* local data handles */
   int nFaces    = distance->surface->nFaces;
   int nVertices = distance->surface->nVertices;
   const size_t*         f = distance->surface->faces;
   const double*         w = distance->surface->weights;
   const double*      heat = distance->heat;
         double* potential = distance->potential.values;

   /* current triangle data */
   double u0, u1, u2; /* heat values */
   double rMag; /* reciprocal of magnitude */
   double *t0, *t1, *t2; /* edge normals */
   double *e0, *e1, *e2; /* cotan-weighted edge vectors */
   hmVec3 X; /* normalized gradient */
   double e0DotX, e1DotX, e2DotX;

   /* initialize potential to zero */
   hmClearArrayDouble( potential, distance->surface->nVertices, 0. );

   /* get pointers to first three edge normals */
   t0 = &distance->surface->edgeNormals[0];
   t1 = &distance->surface->edgeNormals[3];
   t2 = &distance->surface->edgeNormals[6];

   /* get pointers to first three weighted edges */
   e0 = &distance->surface->weightedEdges[0];
   e1 = &distance->surface->weightedEdges[3];
   e2 = &distance->surface->weightedEdges[6];

   /* add contribution from each face */
   for( i = 0; i < nFaces; i++ )
   {
      /* get heat values at three vertices */
      u0 = fabs( heat[ f[0] ] );
      u1 = fabs( heat[ f[1] ] );
      u2 = fabs( heat[ f[2] ] );

      /* normalize heat values so that they have roughly unit magnitude */
      rMag = 1./hmMaxDouble( hmMaxDouble( u0, u1 ), u2 );
      if( !isinf(rMag) )
      {
         u0 *= rMag;
         u1 *= rMag;
         u2 *= rMag;

         /* compute normalized gradient */
         X[0] = u0*t0[0] + u1*t1[0] + u2*t2[0];
         X[1] = u0*t0[1] + u1*t1[1] + u2*t2[1];
         X[2] = u0*t0[2] + u1*t1[2] + u2*t2[2];
         hmVec3Normalize( X );

         /* add contribution to divergence */
         e0DotX = hmVec3Dot( e0, X );
         e1DotX = hmVec3Dot( e1, X );
         e2DotX = hmVec3Dot( e2, X );
         potential[ f[0] ] -= e1DotX - e2DotX;
         potential[ f[1] ] -= e2DotX - e0DotX;
         potential[ f[2] ] -= e0DotX - e1DotX;

         if( isnan( potential[f[0]] ) ||
             isnan( potential[f[1]] ) ||
             isnan( potential[f[2]] ) )
         {
            fprintf( stderr, "NaN\n============\n" );
            fprintf( stderr, "heat: %e %e %e\n", heat[f[0]], heat[f[1]], heat[f[2]] );
            fprintf( stderr, " mag: %e\n",  hmMaxDouble( hmMaxDouble( u0, u1 ), u2 ));
            fprintf( stderr, "rMag: %e\n", rMag );
            fprintf( stderr, "   u: %e %e %e\n", u0, u1, u2 );
            fprintf( stderr, "   X: %e %e %e\n", X[0], X[1], X[2] );
            fprintf( stderr, "ei*X: %e %e %e\n", e0DotX, e1DotX, e2DotX );
            exit( 1 );
         }
      }

      /* move to next face */
      f += 3;
      w += 3;
      t0 += 9; t1 += 9; t2 += 9;
      e0 += 9; e1 += 9; e2 += 9;
   }

   /* remove mean value so that the potential is
    * in the range of the Laplace operator */
   hmRemoveMean( potential, nVertices );
}

void hmTriDistanceSolvePoissonEquation( hmTriDistance* distance )
{
   size_t i;
   const size_t nVertices = distance->surface->nVertices;
   double minDistance = DBL_MAX;
   double* phi;

   hmCholeskyFactorBacksolve( &distance->laplacianFactor,
                              &distance->distance,
                              &distance->potential );

   /* subtract the minimum value */
   phi = distance->distance.values;
   for( i = 0; i < nVertices; i++ )
   {
      minDistance = hmMinDouble( minDistance, phi[i] );
   }
   #pragma omp parallel for
   for( i = 0; i < nVertices; i++ )
   {
      phi[i] -= minDistance;
   }
}

void hmTriDistanceBuildMatrices( hmTriDistance* distance )
{
   size_t i,j;
   size_t nz; /* current nonzero */
   size_t lastNeighborIndex;
   size_t count; /* number of times a given neighbor appears */
   size_t* columnStart;
   char* onBoundary;
   double* columnSum;
   double* x;
   size_t* n;
   const hmPairSizeTDouble *neighborsBegin, *neighborsEnd, *currentNeighbor;
   const double boundaryConditions = distance->boundaryConditions;
   const double time = distance->time;
   double A; /* vertex area */

   hmSparseMatrix* laplacian         = &distance->laplacian;
   hmSparseMatrix* heatFlowNeumann   = &distance->heatFlowNeumann;
   hmSparseMatrix* heatFlowDirichlet = &distance->heatFlowDirichlet;
   hmTriMesh* mesh = distance->surface;

   /* array counters */
   int k,iter;
   int j0, j1, j2;

   /* local data handles */
   size_t nFaces = mesh->nFaces;
   size_t nVertices = mesh->nVertices;
   const size_t* f;
   double* w; /* current weights */
   double* vertices = mesh->vertices;
   double* vertexAreas;
   /* current triangle data */
   double* p[3]; /* vertex positions */
   double *e[3]; /* edge vectors */
   double *t[3]; /* rotated edge vectors */

   hmVectorPairSizeTDouble *neighbors; /* temporary, redundant list of vertex neighbors */
   hmVectorPairSizeTDouble *uniqueNeighbors; /* final list of unique vertex neighbors */
   hmPairSizeTDouble neighbor; /* used to construct a record of the current neighbor */

   hmVec3 u, v; /* edge vectors */
   hmVec3 N; /* triangle normal */
   double uvSinTheta, uvCosTheta;
   hmDestroy( mesh->weights );
   mesh->weights = malloc( 3*nFaces * sizeof(double) );

   /* allocate storage for edge data */
   hmDestroy( mesh->edgeNormals );
   hmDestroy( mesh->weightedEdges );
   mesh->edgeNormals   = malloc( 9*nFaces * sizeof( double ));
   mesh->weightedEdges = malloc( 9*nFaces * sizeof( double ));

   /* initialize vertex areas to zero */
   hmDestroy( mesh->vertexAreas );
   mesh->vertexAreas = malloc( nVertices * sizeof( double ));
   hmClearArrayDouble( mesh->vertexAreas, nVertices, 0. );
   vertexAreas = mesh->vertexAreas;
      /* allocate a list of redundant neighbors for each vertex */
      neighbors = malloc( nVertices * sizeof( hmVectorPairSizeTDouble ));
      hmDestroy( mesh->vertexNeighbors );
      mesh->vertexNeighbors = malloc( nVertices * sizeof(hmVectorPairSizeTDouble) );
      uniqueNeighbors = mesh->vertexNeighbors; /* short name */
      columnStart=calloc((nVertices+1),sizeof(size_t));
      columnSum=calloc( nVertices,sizeof(double));

      #pragma omp parallel for
      for( i = 0; i < nVertices; i++ )
      {
         hmVectorPairSizeTDoubleInitialize( &neighbors[i] );
         hmVectorPairSizeTDoubleInitialize( &uniqueNeighbors[i] );
      }
      /* allocate an array of flags for boundary vertices */
      hmDestroy( mesh->onBoundary );
      mesh->onBoundary = malloc( nVertices * sizeof(char) );
   /* iterate over triangles */
   for(iter=0; iter < nFaces; iter++)
   {
      /* get vertex coordinates */
      f=mesh->faces+3*iter;
      w=mesh->weights+3*iter;

      p[0] = &vertices[ f[0]*3 ];
      p[1] = &vertices[ f[1]*3 ];
      p[2] = &vertices[ f[2]*3 ];
      /* iterate over triangle corners */
      for( k = 0; k < 3; k++ )
      {
         /* get outgoing edge vectors u, v at current corner */
         size_t v_index=f[k];
         t[k]=mesh->edgeNormals+9*iter+3*k;
         e[k]=mesh->weightedEdges+9*iter+3*k;
         j0 = (0+k) % 3;
         j1 = (1+k) % 3;
         j2 = (2+k) % 3;
         hmVec3Sub( u, p[j1], p[j0] );
         hmVec3Sub( v, p[j2], p[j0] );
         hmVec3Sub( e[k], p[j2], p[j1] );
         /* compute (one-half of) the cotangent weight */
         hmVec3Cross( N, u, v );
         hmVec3Cross( t[k], N, e[k] );
         uvSinTheta = hmVec3Norm( N );
         uvCosTheta = hmVec3Dot( u, v );
         vertexAreas[ v_index ] += uvSinTheta/6;
         w[k] = .5 * uvCosTheta / uvSinTheta;
         hmVec3Scale( e[k], w[k] );
      }
      for( k = 0; k < 3; k++ ){
      j1 = (1+k) % 3;
      j2 = (2+k) % 3;
      neighbor.n = f[j1]; neighbor.x = w[j2]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[k] ], neighbor );
      neighbor.n = f[j2]; neighbor.x = w[j1]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[k] ], neighbor );}
   }

   /* iterate over vertices */
   #pragma omp parallel for private(lastNeighborIndex,count,j)
   for( i = 0; i < nVertices; i++ )
   {

      /* sort neighbor list by index */
      hmVectorPairSizeTDoubleSort( &neighbors[i] );

      /* initially flag as an interior vertex */
      mesh->onBoundary[i] = 0;

      /* extract unique elements from neighbor list, summing weights */
      hmVectorPairSizeTDoubleResize( &uniqueNeighbors[i], 0 );
      lastNeighborIndex = -1;
      count = 0;
      for( j = 0; j < neighbors[i].size; j++ )
      {
         /* if we come across a new neighbor, add it to the list of unique neighbors */
         if( neighbors[i].entries[j].n != lastNeighborIndex )
         {
            /* if we encountered the previous neighbor only
             * once, this vertex must be on the surface boundary */
            if( count == 1 )
            {
               mesh->onBoundary[i] = 1;
            }
            count = 1;
            if(neighbors[i].entries[j].n>i){
              hmVectorPairSizeTDoublePushBack( &uniqueNeighbors[i], neighbors[i].entries[j] );
              columnStart[i+1]++;
            }
            lastNeighborIndex = neighbors[i].entries[j].n;
         }
         else
         {
            /* since we've seen this neighbor before, just accumulate its weight */
            if(neighbors[i].entries[j].n>i)
            uniqueNeighbors[i].entries[ uniqueNeighbors[i].size-1 ].x += neighbors[i].entries[j].x;
            count++;
         }
         columnSum[i]+=neighbors[i].entries[j].x;
      }

      /* if the final neighbor was encountered only once, this is a boundary vertex */
      if( count == 1 )
      {
         mesh->onBoundary[i] = 1;
      }
      hmVectorPairSizeTDoubleDestroy( &neighbors[i] );
   }
   free( neighbors );

   /* determine the starting entry of nonzeros in
    * each column, keeping the lower triangle only */
    for( i = 1; i < (nVertices+1); ++i)
     {
         columnStart[i]+=columnStart[i-1];
     }
   x=malloc( columnStart[nVertices] * sizeof(double));
   n=malloc( columnStart[nVertices] * sizeof(size_t));
   onBoundary = distance->surface->onBoundary;

   /* initialize matrices and copy column start pointers */
   hmSparseMatrixDestroy( laplacian );
   hmSparseMatrixInitialize( laplacian, nVertices, nVertices, (columnStart[nVertices]+nVertices) );
   if( boundaryConditions < 1. ) /* partial Neumann */
   {
      hmSparseMatrixDestroy( heatFlowNeumann );
      hmSparseMatrixInitialize( heatFlowNeumann,  nVertices, nVertices, (columnStart[nVertices]+nVertices) );}
   if( boundaryConditions > 0. ) /* partial Dirichlet */
   {
         hmSparseMatrixDestroy( heatFlowDirichlet );
         hmSparseMatrixInitialize( heatFlowDirichlet,  nVertices, nVertices, (columnStart[nVertices]+nVertices) );}

   for( i = 0; i < nVertices+1; i++ )
   {
      distance->laplacian.columnStart[i] = columnStart[i]+i;
      if( boundaryConditions < 1. ) /* partial Neumann */
      distance->heatFlowNeumann.columnStart[i] = columnStart[i]+i;
      if( boundaryConditions > 0. )/* partial Dirichlet */
      distance->heatFlowDirichlet.columnStart[i] = columnStart[i]+i;
      neighborsBegin = uniqueNeighbors[i].entries;
      neighborsEnd   = neighborsBegin + uniqueNeighbors[i].size;
      for( currentNeighbor  = neighborsBegin;
           currentNeighbor != neighborsEnd;
           currentNeighbor ++ )
      {
         x[columnStart[i]+currentNeighbor-neighborsBegin]=currentNeighbor->x;
         n[columnStart[i]+currentNeighbor-neighborsBegin]=currentNeighbor->n;
      }
   }

   /* fill nonzero entries */
   for( i = 0; i < nVertices; i++ )
   {
            /* set diagonal entry of Laplacian, adding a small
               regularization term in order to get strict
               positive-definiteness (needed for CHOLMOD) */
            laplacian->values[ columnStart[i]+i] = columnSum[i] + hmRegularization;
            laplacian->rowIndices[columnStart[i]+i] = i;

            A = distance->surface->vertexAreas[ i ];

            if( boundaryConditions < 1. ) /* partial Neumann */
            {
               heatFlowNeumann->values[columnStart[i]+i] = A + time*columnSum[i];
               heatFlowNeumann->rowIndices[columnStart[i]+i] = i;
            }
            if( boundaryConditions > 0. ) /* partial Dirichlet */
            {
               if( onBoundary[ i ] )
               {
                  /* use the identity (times the mass matrix) for boundary
                   * rows/columns to enforce zero-Dirichlet conditions */
                  heatFlowDirichlet->values[columnStart[i]+i] = A;
               }
               else
               {
                  heatFlowDirichlet->values[columnStart[i]+i] = A + time*columnSum[i];
               }
               heatFlowDirichlet->rowIndices[columnStart[i]+i] = i;
            }

         /* set off-diagonal entries below the diagonal */
        for(nz=columnStart[i]+i+1;nz<columnStart[i+1]+i+1;nz++)
         {
            laplacian->values[ nz ] = -x[nz-i-1];
            laplacian->rowIndices[ nz ] = n[nz-i-1];

            if( boundaryConditions < 1. ) /* partial Neumann */
            {
               heatFlowNeumann->values[ nz ] = -time*x[nz-i-1];
               heatFlowNeumann->rowIndices[ nz ] = n[nz-i-1];
            }
            if( boundaryConditions > 0. ) /* partial Dirichlet */
            {
               if( onBoundary[i] || onBoundary[n[nz-i-1]] )
               {
                  /* set off-diagonals to zero so that we retain
                   * the same sparsity pattern as other matrices */
                  heatFlowDirichlet->values[ nz ] = 0.;
               }
               else
               {
                  heatFlowDirichlet->values[ nz ] = -time*x[nz-i-1];
               }
               heatFlowDirichlet->rowIndices[ nz ] = n[nz-i-1];
            }
         }
      }
}

void hmTriDistanceFactorMatrices( hmTriDistance* distance )
{
   /* Laplacian */
   hmCholeskyFactorDestroy    ( &distance->laplacianFactor );
   hmCholeskyFactorInitialize ( &distance->laplacianFactor );
   hmCholeskyFactorReorder    ( &distance->laplacianFactor, &distance->laplacian );
   hmCholeskyFactorSymbolic   ( &distance->laplacianFactor, &distance->laplacian );
   hmCholeskyFactorNumerical  ( &distance->laplacianFactor, &distance->laplacian );

   /* only factor both heat flow operators if necessary */
   /* (note that the symbolic factorization for Laplace can be reused in both
    * cases since all three matrices have the same sparsity pattern) */
   if( distance->boundaryConditions < 1. ) /* partial Neumann */
   {
      hmCholeskyFactorDestroy    ( &distance->heatFlowNeumannFactor );
      hmCholeskyFactorInitialize ( &distance->heatFlowNeumannFactor );
      hmCholeskyFactorCopy       ( &distance->heatFlowNeumannFactor, &distance->laplacianFactor );
#ifdef HM_USE_HSLMA87 /* currently no way to copy symbolic factorization in HSL_MA87... */
      hmCholeskyFactorSymbolic   ( &distance->heatFlowNeumannFactor, &distance->heatFlowNeumann );
#endif
      hmCholeskyFactorNumerical  ( &distance->heatFlowNeumannFactor, &distance->heatFlowNeumann );
   }
   if( distance->boundaryConditions > 0. ) /* partial Dirichlet */
   {
      hmCholeskyFactorDestroy    ( &distance->heatFlowDirichletFactor );
      hmCholeskyFactorInitialize ( &distance->heatFlowDirichletFactor );
      hmCholeskyFactorCopy       ( &distance->heatFlowDirichletFactor, &distance->laplacianFactor );
#ifdef HM_USE_HSLMA87
#else
      hmCholeskyFactorSymbolic   ( &distance->heatFlowDirichletFactor, &distance->heatFlowDirichlet );
#endif
      hmCholeskyFactorNumerical  ( &distance->heatFlowDirichletFactor, &distance->heatFlowDirichlet );
   }
}

void hmTriDistanceUpdateTime( hmTriDistance* distance, double time )
{
   hmTriDistanceStartTiming( distance );

   distance->time = time;

   hmTriDistanceBuildMatrices( distance );

   hmCholeskyFactorNumerical( &distance->laplacianFactor, &distance->laplacian );

   if( distance->boundaryConditions < 1. ) /* partial Neumann */
   {
      hmCholeskyFactorNumerical( &distance->heatFlowNeumannFactor, &distance->heatFlowNeumann );
   }

   if( distance->boundaryConditions > 0. ) /* partial Dirichlet */
   {
      hmCholeskyFactorNumerical( &distance->heatFlowDirichletFactor, &distance->heatFlowDirichlet );
   }

   hmTriDistanceStopTiming( distance, "Update t parameter" );
}

void hmTriDistanceStartTiming( hmTriDistance* distance )
{
   struct timeval t;
   struct timezone tz;
   double startTimeProcessor;
   double startTimeWallClock;

   if( !distance->verbose )
   {
      return;
   }

   startTimeProcessor = (double) clock() / (double) CLOCKS_PER_SEC;
   hmVectorDoublePushBack( &distance->startTimesProcessor, startTimeProcessor );

   gettimeofday( &t, &tz );
   startTimeWallClock = (double) t.tv_sec + 1e-6*(double) t.tv_usec;
   hmVectorDoublePushBack( &distance->startTimesWallClock, startTimeWallClock );
}

void hmTriDistanceStopTiming( hmTriDistance* distance,
                              const char* label )
{
   int i;
   struct timeval t;
   struct timezone tz;
   double startTimeProcessor, stopTimeProcessor;
   double startTimeWallClock, stopTimeWallClock;
   int nTabs = distance->startTimesProcessor.size-1;

   if( !distance->verbose )
   {
      return;
   }

   startTimeProcessor = hmVectorDoublePopBack( &distance->startTimesProcessor );
    stopTimeProcessor = (double) clock() / (double) CLOCKS_PER_SEC;

   gettimeofday( &t, &tz );
   startTimeWallClock = hmVectorDoublePopBack( &distance->startTimesWallClock );
    stopTimeWallClock = (double) t.tv_sec + 1e-6*(double) t.tv_usec;

   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( "%s\n", label );
   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( "--------------------------------------------\n" );
   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( " processor time: %f seconds\n", stopTimeProcessor-startTimeProcessor );
   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( "wall-clock time: %f seconds\n", stopTimeWallClock-startTimeWallClock );
   printf( "\n" );
}
