struct FPOINT
{
	double x;
	double y;
};

struct LINESEGMENT
{
	struct FPOINT a;
	struct FPOINT b;
};

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define LineTest(one, two) (-((one.b.y - two.b.y)*(two.a.x - two.b.x) - (one.b.x - two.b.x)*(two.a.y - two.b.y))/((one.a.y - one.b.y)*(two.a.x - two.b.x) - (one.a.x - one.b.x)*(two.a.y - two.b.y)))
#define Parallel(one, two) (((one.a.y - one.b.y)*(two.a.x - two.b.x) - (one.a.x - one.b.x)*(two.a.y - two.b.y)) == 0.0)
#define GetX(line, A) (line.a.x * A + line.b.x * (1.0 - A))
#define GetY(line, A) (line.a.y * A + line.b.y * (1.0 - A))
#define DoubleEquals(one, two) (fabs(one - two) < 0.0000000001)
#define LineEquals(one, two) ((DoubleEquals(one.a.x, two.a.x) && DoubleEquals(one.a.y, two.a.y) && DoubleEquals(one.b.x, two.b.x) && DoubleEquals(one.b.y, two.b.y)) || (DoubleEquals(one.a.x, two.b.x) && DoubleEquals(one.a.y, two.b.y) && DoubleEquals(one.b.x, two.a.x) && DoubleEquals(one.b.y, two.a.y)))
#define PointEquals(one, two) (DoubleEquals(one.x, two.x) && DoubleEquals(one.y, two.y))
#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)

int pnt_in_poly_(int *nVert, double *pnt, double *Poly_x, double *Poly_y)
{
//-----------------------------------------------------------
// Is a point inside a polygon?
//-----------------------------------------------------------
  struct LINESEGMENT TestLine;
  struct LINESEGMENT PolyEdge;
  int    nVtx, inside, nXs;
  double minX,maxX,minY,maxY;
  double a;
  double aStore[4];
  double aArray[10][2];

  nVtx = *nVert;
//-----------------------------------------------------------
// Define the bounding box for the polygon
//-----------------------------------------------------------
  minX = Poly_x[0];
  maxX = Poly_x[0];
  minY = Poly_y[0];
  maxY = Poly_y[0];
  for( int i = 1 ; i < nVtx ; i++ )
  {
    minX = min(minX,Poly_x[i]);
    maxX = max(maxX,Poly_x[i]);
    minY = min(minY,Poly_y[i]);
    maxY = max(maxY,Poly_y[i]);
  }

//-----------------------------------------------------------
// is pnt outside box bounding polygon?
//-----------------------------------------------------------
  if( pnt[0] < minX || pnt[0] > maxX )
    return 0;
  else if( pnt[1] < minY || pnt[1] > maxY )
    return 0;
//-----------------------------------------------------------
// is pnt a vertex of the polygon?
//-----------------------------------------------------------
  struct FPOINT matchPnt, vertexPnt;
  matchPnt.x = pnt[0] ; matchPnt.y = pnt[1];
  for(int i = 0; i < nVtx; i++)
  {
   vertexPnt.x = Poly_x[i]; vertexPnt.y = Poly_y[i];
   if( PointEquals( matchPnt,vertexPnt) ) 
     return 1;
  }

//-----------------------------------------------------------
// form the "test" line
//-----------------------------------------------------------
  TestLine.a.x = pnt[0];
  TestLine.a.y = pnt[1];
  TestLine.b.x = maxX+1.;
  TestLine.b.y = pnt[1];

  inside = 0;
  nXs = 0;
  for( int i = 0 ; i < nVtx ; i++ )
  {
    PolyEdge.a.x = Poly_x[i];
    PolyEdge.a.y = Poly_y[i];
    PolyEdge.b.x = Poly_x[(i+1)%nVtx];
    PolyEdge.b.y = Poly_y[(i+1)%nVtx];
    aArray[i][1] = 0.0;
    if( !Parallel(PolyEdge,TestLine) )
    {
//    double a = LineTest(PolyEdge,TestLine);
      a = LineTest(PolyEdge,TestLine);
      aArray[i][0] = a;
      if( a >= 0.0 && a <= 1.0)
      {
        a = LineTest(TestLine,PolyEdge);
        aArray[i][1] = a;
        if( DoubleEquals(a,1.0) )
        {
          inside = 1;
          nXs = 0;
          break;
        }
        if( a >= 0.0 && a <= 1.0)
        {
          aStore[nXs++] = a;
          inside = !inside;
        }
      }
    }
  }

//-----------------------------------------------------------
// Check for point doubles
//-----------------------------------------------------------
  if( nXs > 4 )
  {
    printf("\nXsecting line count = %d\n",nXs);
    printf("\nnVtx = %d\n",nVtx);
    printf("\nPnt = (%f,%f)\n",pnt[0],pnt[1]);
    printf("\nTestLine\n");
    printf("(%f,%f)\n",TestLine.a.x,TestLine.a.y);
    printf("(%f,%f)\n",TestLine.b.x,TestLine.b.y);
    printf("\nVertices\n");
    for(int i = 0; i < nVtx; i++)
      printf("(%f,%f)\n",Poly_x[i],Poly_y[i]);
    printf("\nA values\n");
    for(int i = 0; i < nXs; i++)
      printf("%f\n",aStore[i]);
    printf("\nTerminating\n");
    exit( -1 );
  }
  else if( nXs > 1 )
  {
    for(int i = 0; i < nXs-1; i++)
      for(int j = i+1; j < nXs; j++)
        if( DoubleEquals(aStore[i],aStore[j]) )
          inside = !inside;
  }
  return(inside);
}

double poly_area_(int *nVtx, double *x, double *y)
{
   int i,im1,ip1,nv;
   double accum;
   double term, accump, accumm;

   nv    = *nVtx;
   accum = 0.0;
   accump = 0.0; accumm = 0.0;
   for(i = 0; i < nv; i++ )
   {
     ip1 = (i+1)%nv;
     im1 = i - 1;
     if( im1 < 0 ) im1 = nv - 1;
//   accum += (x[ip1] - x[im1])*y[i];
     term = (x[ip1] - x[im1])*y[i];
     if( term > 0.0 ) 
       accump += term;  
     else if( term < 0.0 )
       accumm += term;  
   }
// return( -.5*accum ); 
   return( -.5*(accump + accumm) ); 
   
}
