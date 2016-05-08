#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hmTriMesh.h"
#include "hmVec3.h"
#include "hmUtility.h"

void hmTriMeshInitialize( hmTriMesh* mesh )
{
   mesh->nVertices = 0;
   mesh->vertices  = NULL;
   mesh->texCoords = NULL;

   mesh->nFaces = 0;
   mesh->faces  = NULL;

   mesh->referenced = 0;

   mesh->onBoundary      = NULL;
   mesh->vertexNeighbors = NULL;
   mesh->vertexAreas     = NULL;
   mesh->weights         = NULL;
   mesh->edgeNormals     = NULL;
   mesh->weightedEdges   = NULL;
}

void hmTriMeshDestroy( hmTriMesh* mesh )
{
   if( !mesh->referenced )
   {
      if( mesh->vertices != NULL )
      {
         free( mesh->vertices );
      }

      if( mesh->faces != NULL )
      {
         free( mesh->faces );
      }
   }

   hmDestroy( mesh->onBoundary );
   hmDestroy( mesh->vertexNeighbors );
   hmDestroy( mesh->vertexAreas );
   hmDestroy( mesh->weights );
   hmDestroy( mesh->edgeNormals );
   hmDestroy( mesh->weightedEdges );
   hmDestroy( mesh->texCoords );

   mesh->vertices = NULL;
   mesh->texCoords = NULL;
   mesh->faces = NULL;
   mesh->referenced = 0;
}

void hmTriMeshCopy( hmTriMesh* mesh1, const hmTriMesh* mesh2 )
{
   hmTriMeshDestroy( mesh1 );
   hmTriMeshCopyData( mesh1, mesh2->nVertices, mesh2->vertices,
                             mesh2->nFaces,    mesh2->faces );

   mesh1->texCoords = malloc( mesh2->nVertices*2 * sizeof(double) );
   memcpy( mesh1->texCoords, mesh2->texCoords, mesh2->nVertices*2 * sizeof(double) );
}

void hmTriMeshCopyData( hmTriMesh* mesh,
                        size_t nVertices, const double* vertices,
                        size_t nFaces,    const size_t* faces )
{
   hmTriMeshDestroy( mesh );

   mesh->nVertices = nVertices;
   mesh->vertices = malloc( nVertices*3 * sizeof(double) );
   memcpy( mesh->vertices, vertices, nVertices*3 * sizeof(double) );

   mesh->nFaces = nFaces;
   mesh->faces = malloc( nFaces*3 * sizeof(size_t) );
   memcpy( mesh->faces, faces, nFaces*3 * sizeof(size_t) );

   mesh->referenced = 0;

   /* initialize texture coords to zero */
   mesh->texCoords = malloc( nVertices*2 * sizeof(double) );
   hmClearArrayDouble( mesh->texCoords, nVertices*2, 0. );
}

void hmTriMeshReferenceData( hmTriMesh* mesh,
                             size_t nVertices, double* vertices,
                             size_t nFaces,    size_t* faces )
{
   hmTriMeshDestroy( mesh );

   mesh->nVertices = nVertices;
   mesh->vertices = vertices;

   mesh->nFaces = nFaces;
   mesh->faces = faces;

   mesh->referenced = 1;

   /* initialize texture coords to zero */
   mesh->texCoords = malloc( nVertices*2 * sizeof(double) );
   hmClearArrayDouble( mesh->texCoords, nVertices*2, 0. );
}

void hmTriMeshReadOBJ( hmTriMesh* mesh, const char* filename )
{
   FILE* in;
   char* line;
   char token[32];
   size_t size;
   double* v;
   size_t* f;
   unsigned long I0, I1, I2;

   hmTriMeshDestroy( mesh );
   hmTriMeshInitialize( mesh );

   if( !( in = fopen( filename, "r" )))
   {
      fprintf( stderr, "Error: could not read from file %s\n", filename );
      exit( 1 );
   }

   /* count the number of vertices, faces */
   while( !feof( in ))
   {
      line = fgetln( in, &size );
      if( line == NULL ) continue;

      sscanf( line, "%s", token );

      if( !strcmp( token, "v" ))
      {
         mesh->nVertices++;
      }
      else if( !strcmp( token, "f" ))
      {
         mesh->nFaces++;
      }
   }

   /* allocate storage */
   mesh->vertices  = malloc( mesh->nVertices*3 * sizeof(double) );
   mesh->texCoords = malloc( mesh->nVertices*2 * sizeof(double) );
   mesh->faces     = malloc(    mesh->nFaces*3 * sizeof(size_t) );

   /* read the mesh data */
   rewind( in );
   v = mesh->vertices;
   f = mesh->faces;
   while( !feof( in ))
   {
      line = fgetln( in, &size );
      if( line == NULL ) continue;

      sscanf( line, "%s", token );

      if( !strcmp( token, "v" ))
      {
         sscanf( line, "%*s %lf %lf %lf", &v[0], &v[1], &v[2] );
         v += 3;
      }
      else if( !strcmp( token, "f" ))
      {
         /* try reading triangle vertex indices, in several possible formats */
         if( 3 == sscanf( line, "%*s %lu %lu %lu", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu/%*d %lu/%*d %lu/%*d", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu/%*d/%*d %lu/%*d/%*d %lu/%*d/%*d", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu//%*d %lu//%*d %lu//%*d", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu// %lu// %lu//", &I0, &I1, &I2 ))
         {
            /* change 1-based indices to 0-based indices */
            f[0] = I0-1;
            f[1] = I1-1;
            f[2] = I2-1;
         }
         else
         {
            fprintf( stderr, "Error: could not parse face line %s\n", line );
            fprintf( stderr, "(in file %s)\n", filename );
            exit( 1 );
         }

         f += 3;
      }
   }

   fclose( in );
}

void hmTriMeshWriteOBJ( const hmTriMesh* mesh, const char* filename )
{
   FILE *out;
   double* v = mesh->vertices;
   double* vt = mesh->texCoords;
   size_t* f = mesh->faces;
   size_t i;
   unsigned long I0, I1, I2;

   if( !( out = fopen( filename, "w" )))
   {
      fprintf( stderr, "Warning: could not write to file %s\n", filename );
      return;
   }

   for( i = 0; i < mesh->nVertices; i++ )
   {
      fprintf( out, "v %.9f %.9f %.9f\n", v[0], v[1], v[2] );
      v += 3;
   }

   for( i = 0; i < mesh->nVertices; i++ )
   {
      fprintf( out, "vt %.9f %.9f\n", vt[0], vt[1] );
      vt += 2;
   }

   for( i = 0; i < mesh->nFaces; i++ )
   {
      I0 = 1+f[0];
      I1 = 1+f[1];
      I2 = 1+f[2];
      fprintf( out, "f %lu/%lu %lu/%lu %lu/%lu\n", I0, I0, I1, I1, I2, I2 );
      f += 3;
   }
}


void hmTriMeshComputeVertexAreas( hmTriMesh* mesh )
{
   size_t nVertices = mesh->nVertices;
   size_t nFaces = mesh->nFaces;

   const size_t* facesBegin = mesh->faces;
   const size_t* facesEnd = facesBegin + 3*nFaces;
   const size_t* f;

   double* vertexAreas;
   double* vertices = mesh->vertices;
   double *p0, *p1, *p2;
   hmVec3 u, v, w;

   double A;

   /* initialize vertex areas to zero */
   hmDestroy( mesh->vertexAreas );
   mesh->vertexAreas = malloc( nVertices * sizeof( double ));
   hmClearArrayDouble( mesh->vertexAreas, nVertices, 0. );
   vertexAreas = mesh->vertexAreas;

   /* iterate over faces */
   for( f = facesBegin; f != facesEnd; f += 3 )
   {
      /* get vertex coordinates p0, p1, p2 */
      p0 = &vertices[ f[0]*3 ];
      p1 = &vertices[ f[1]*3 ];
      p2 = &vertices[ f[2]*3 ];

      /* compute (one-third of) the triangle area A = |u x v| */
      hmVec3Sub( u, p1, p0 );
      hmVec3Sub( v, p2, p0 );
      hmVec3Cross( w, u, v );
      A = hmVec3Norm( w ) / 6.;

      /* add contribution to each of the three corner vertices */
      vertexAreas[ f[0] ] += A;
      vertexAreas[ f[1] ] += A;
      vertexAreas[ f[2] ] += A;
   }
}

double hmTriMeshL2Distance( hmTriMesh* mesh, double* phi1, double* phi2 )
{
   size_t i;
   double sum = 0.;

   hmTriMeshComputeVertexAreas( mesh );

   for( i = 0; i < mesh->nVertices; i++ )
   {
      sum += mesh->vertexAreas[i] * hmSquare( phi1[i] - phi2[i] );
   }

   return sqrt( sum );
}

double hmTriMeshLInfinityDistance( hmTriMesh* mesh, double* phi1, double* phi2 )
{
   size_t i;
   double maxDifference = 0.;

   for( i = 0; i < mesh->nVertices; i++ )
   {
      maxDifference = hmMaxDouble( maxDifference, fabs( phi1[i]-phi2[i] ));
   }

   return maxDifference;
}

double hmTriMeshMeanRelativeError( hmTriMesh* mesh, double* phi1, double* phi2 )
{
   size_t i;
   double mean = 0.;

   for( i = 0; i < mesh->nVertices; i++ )
   {
      if( phi1[i] != 0. )
      {
         mean += fabs( (phi2[i]-phi1[i]) / phi1[i] );
      }
   }

   mean /= (double) mesh->nVertices;

   return mean;
}

double hmTriMeshDiameter( hmTriMesh* mesh )
{
   size_t nVertices = mesh->nVertices;
   size_t nFaces = mesh->nFaces;

   const size_t* facesBegin = mesh->faces;
   const size_t* facesEnd = facesBegin + 3*nFaces;
   const size_t* f;

   const double* verticesBegin = mesh->vertices;
   const double* verticesEnd = verticesBegin + 3*nVertices;
   const double* v;

   double* vertices = mesh->vertices;
   double *p0, *p1, *p2;
   hmVec3 u1, u2, w;
   hmVec3 barycenter;

   double radius;
   double triangleArea, surfaceArea;

   hmVec3 centerOfMass;
   hmVec3Set( centerOfMass, 0., 0., 0. );

   /* iterate over faces */
   surfaceArea = 0.;
   for( f = facesBegin; f != facesEnd; f += 3 )
   {
      /* get vertex coordinates p0, p1, p2 */
      p0 = &vertices[ f[0]*3 ];
      p1 = &vertices[ f[1]*3 ];
      p2 = &vertices[ f[2]*3 ];

      /* compute barycenter (p0+p1+p2)/3 */
      hmVec3Set( barycenter, 0., 0., 0. );
      hmVec3Add( barycenter, barycenter, p0 );
      hmVec3Add( barycenter, barycenter, p1 );
      hmVec3Add( barycenter, barycenter, p2 );
      hmVec3Scale( barycenter, 1./3. );

      /* compute the triangle area A = |u x v| */
      hmVec3Sub( u1, p1, p0 );
      hmVec3Sub( u2, p2, p0 );
      hmVec3Cross( w, u1, u2 );
      triangleArea = hmVec3Norm( w ) / 2.;

      /* add contribution to total area */
      surfaceArea += triangleArea;

      /* add contribution to center of mass */
      hmVec3Scale( barycenter, triangleArea );
      hmVec3Inc( centerOfMass, barycenter );
   }

   hmVec3Scale( centerOfMass, 1./surfaceArea );

   /* iterate over vertices */
   radius = 0.;
   for( v = verticesBegin; v != verticesEnd; v+=3 )
   {
      hmVec3Sub( w, v, centerOfMass );
      radius = hmMaxDouble( radius, hmVec3Norm( w ));
   }

   return 2.*radius;
}
