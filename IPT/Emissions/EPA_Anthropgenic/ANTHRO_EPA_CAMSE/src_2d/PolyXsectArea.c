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

int pnt_in_poly(int nVtx, double *pnt, double *Poly_x, double *Poly_y)
{
//-----------------------------------------------------------
// Is a point inside a polygon?
//-----------------------------------------------------------
  struct LINESEGMENT TestLine;
  struct LINESEGMENT PolyEdge;
  int inside, nXs;
  double minX,maxX,minY,maxY;
  double aStore[4];

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
  TestLine.b.x = maxX*1.1;
  TestLine.b.y = pnt[1];

  inside = 0;
  nXs = 0;
  for( int i = 0 ; i < nVtx ; i++ )
  {
    PolyEdge.a.x = Poly_x[i];
    PolyEdge.a.y = Poly_y[i];
    PolyEdge.b.x = Poly_x[(i+1)%nVtx];
    PolyEdge.b.y = Poly_y[(i+1)%nVtx];
    if( !Parallel(PolyEdge,TestLine) )
    {
      double a = LineTest(PolyEdge,TestLine);
      if( a >= 0.0 && a <= 1.0)
      {
        a = LineTest(TestLine,PolyEdge);
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

double poly_area(int nv, double *x, double *y)
{
   int i,im1,ip1;
   double accum;

   accum = 0.0;
   for(i = 0; i < nv; i++ )
   {
     ip1 = (i+1)%nv;
     im1 = i - 1;
     if( im1 < 0 ) im1 = nv - 1;
     accum += (x[ip1] - x[im1])*y[i];
   }
   return( -.5*accum ); 
   
}

void cleanup( struct LINESEGMENT *p1, struct LINESEGMENT *p2, struct LINESEGMENT *lines, int *mask )
{
   if( p1 != NULL) free(p1);
   if( p2 != NULL) free(p2);
   if( lines != NULL) free(lines);
   if( mask != NULL) free(mask);
}

double polyintersectarea_(int *nData_n_Mdl, int *nMdl_n_Data, int *vtxcnt1, int *vtxcnt2, double *Poly1_x, double *Poly1_y, double *Poly2_x, double *Poly2_y)
{
   struct LINESEGMENT *Poly1=NULL, *Poly2=NULL;
   struct LINESEGMENT *lines=NULL;
   int *mask=NULL;
   int lineIndex;
   int i, n, p, r;
   int cnt;
   int sLineSeg;

   int nVtx1 = *vtxcnt1;
   int nVtx2 = *vtxcnt2;

   double vertex[2];
//-----------------------------------------------------------
// check for complete inclusion, first Poly1 in Poly2
//-----------------------------------------------------------
   p = 0;
   r = 0;
   for (int i = 0; i < nVtx1; i++)    // Poly one in Poly two
   {
     vertex[0] = Poly1_x[i];
     vertex[1] = Poly1_y[i];
     n = pnt_in_poly(nVtx2,vertex,Poly2_x,Poly2_y);
     if( n == 0 ) 
       break;
     else
       p++;
   }
   r = (p == nVtx1) ? 0 : -1;
//-----------------------------------------------------------
// then Poly2 in Poly1
//-----------------------------------------------------------
   if( r != 0 ) {
     p = 0;
     for (int i = 0; i < nVtx2; i++)    // Poly two in Poly one
     {
       vertex[0] = Poly2_x[i];
       vertex[1] = Poly2_y[i];
       n = pnt_in_poly(nVtx1,vertex,Poly1_x,Poly1_y);
       if( n == 0 ) 
         break;
       else
         p++;
     }
     r = (p == nVtx2) ? 1 : -1;
   }

   if( r != -1 )
   {
     double area;
     if( r == 0 )
     {
       area = poly_area(nVtx1, Poly1_x, Poly1_y);
       (*nData_n_Mdl)++;
     }
     else
     {
       area = poly_area(nVtx2, Poly2_x, Poly2_y);
       (*nMdl_n_Data)++;
     }
     return area;
   }
//-----------------------------------------------------------
// Allocate polygons and associated vectors
//-----------------------------------------------------------
   {
     int nullCnt = 0;
     sLineSeg = sizeof(struct LINESEGMENT);
     Poly1 = (struct LINESEGMENT *)malloc((size_t)sLineSeg*nVtx1);
     if( Poly1 == NULL ) nullCnt++;
     Poly2 = (struct LINESEGMENT *)malloc((size_t)sLineSeg*nVtx2);
     if( Poly2 == NULL ) nullCnt++;
     lines = (struct LINESEGMENT *)malloc((size_t)sLineSeg*nVtx1*nVtx2);
     if( lines == NULL ) nullCnt++;
     mask  = (int *)malloc((size_t)sLineSeg*nVtx1*nVtx2);
     if( mask == NULL ) nullCnt++;
     if( nullCnt != 0 )
     {
       cleanup( Poly1, Poly2, lines, mask );
       return 0.0;
     }
   }

//-----------------------------------------------------------
// Convert from points to line segments...
//-----------------------------------------------------------
//    Poly one
   for (int i = 0; i < nVtx1; i++)
   {
      Poly1[i].a.x = Poly1_x[i];           // Point A
      Poly1[i].a.y = Poly1_y[i];
      Poly1[i].b.x = Poly1_x[(i+1)%nVtx1]; // Point B
      Poly1[i].b.y = Poly1_y[(i+1)%nVtx1];
   }
//    Poly two
   for (int i = 0; i < nVtx2; i++)
   {
      Poly2[i].a.x = Poly2_x[i];           // Point A
      Poly2[i].a.y = Poly2_y[i];
      Poly2[i].b.x = Poly2_x[(i+1)%nVtx2]; // Point B
      Poly2[i].b.y = Poly2_y[(i+1)%nVtx2];
   }

   lineIndex = 0;
//-----------------------------------------------------------
// Get the intersection line segments for the polygons ...
//-----------------------------------------------------------
   for (int polyIndex = 0; polyIndex < 2; polyIndex++)
   {
      struct LINESEGMENT MatchLine;
      struct LINESEGMENT *polyTst = NULL;
      struct LINESEGMENT *matchPoly = NULL;

      int nMatchVtx  = (polyIndex == 0) ? nVtx1 : nVtx2;
      int nTstVtx    = (polyIndex == 0) ? nVtx2 : nVtx1;
      matchPoly = (polyIndex == 0) ? Poly1 : Poly2;
      polyTst   = (polyIndex == 0) ? Poly2 : Poly1;
      
      for (int i = 0; i < nMatchVtx; i++)
      {
         struct FPOINT intersect[4];
         double aArray[4];
         int index = 0;
         MatchLine = matchPoly[i];
//-----------------------------------------------------------
// Check intersctions for every line in other quad...
//-----------------------------------------------------------
         for (int n = 0; n < nTstVtx; n++)
         {
            if (!Parallel(polyTst[n],MatchLine))
	    {
	       double a = LineTest(polyTst[n], MatchLine);
	       if (a >= 0.0 && a <= 1.0)
               {
		  a = LineTest(MatchLine, polyTst[n]);
		  intersect[index].x = GetX(MatchLine, a);
		  intersect[index].y = GetY(MatchLine, a);
		  aArray[index++] = a;
	       }
	    }
	 }
	 if (index == 0)
	    continue;
         else if( index > 4 )
         {
           printf("\nXsection count = %d\n",index);
           printf("\nTerminating\n");
           exit( -1 );
         }
//-----------------------------------------------------------
// Remove point doubles...
//-----------------------------------------------------------
	 for (int p = 0; p < index; p++) mask[p] = 0;

	 for (int p = 0; p < index-1; p++)
           if( !mask[p] )
           {
	    for (int r = p+1; r < index; r++)
	      if( !mask[r] )
	        mask[r] = PointEquals(intersect[r], intersect[p]) ? 1 : 0;
	   }

         cnt = 1;
	 for (int p = 1; p < index; p++)
         {
           if( !mask[p] )
           {
             intersect[cnt] = intersect[p];
             aArray[cnt++]  = aArray[p];
           }
         }
         index = cnt;
	 if (index == 0)
	    continue;
        
//-----------------------------------------------------------
// Sort the lists smallest a value to largest a value...
//-----------------------------------------------------------
	 struct FPOINT sortedIntersect[4];
	 double sortedA[4];
	 for (int n = 0; n < index; n++)
	 {
	   int inum;
	   for (inum = 0; inum < n; inum++)
	   {
	     if (aArray[n] < sortedA[inum])
	     {
//-----------------------------------------------------------
// shift list to fit segment in correct place...
//-----------------------------------------------------------
	       for (int p = n; p > inum; p--)
	       {
		 sortedA[p] = sortedA[p - 1];
		 sortedIntersect[p] = sortedIntersect[p - 1];
	       }
	       sortedA[inum] = aArray[n];
	       sortedIntersect[inum] = intersect[n];
	       break;
	     }
	   }
	   if (inum == n)
	   {
	     sortedA[inum] = aArray[n];
	     sortedIntersect[inum] = intersect[n];
	   }
	 }
//-----------------------------------------------------------
// create line segments and clip line segments
// to fit current quad and add to a list...
//-----------------------------------------------------------
	 for (int n = 0; n < index; n+=2)
	 {
           double sAn, sAnp1;
           sAn = sortedA[n];
           sAnp1 = sortedA[n+1];
	   if (sAn < 0.0)
	   {
	      if (sAnp1 > 1.0)
	        lines[lineIndex++] = MatchLine;
	      else if (sAnp1 >= 0.0)
	      {
	        lines[lineIndex].a.x = GetX(MatchLine, sAnp1);
		lines[lineIndex].a.y = GetY(MatchLine, sAnp1);
		lines[lineIndex++].b = MatchLine.b;
	      }
	    }
	    else if (sAn > 1.0)
	    {
	      if (sAnp1 <= 1.0 && sAnp1 >= 0.0)
	      {
	        lines[lineIndex].b.x = GetX(MatchLine, sAnp1);
		lines[lineIndex].b.y = GetY(MatchLine, sAnp1);
		lines[lineIndex++].a = MatchLine.a;
	      }
	    }
	    else
	    {
	      if (sAnp1 > 1.0)
	      {
		lines[lineIndex].b.x = GetX(MatchLine, sAn);
		lines[lineIndex].b.y = GetY(MatchLine, sAn);
		lines[lineIndex++].a = MatchLine.a;
	      }
	      else if (sAnp1 < 0.0)
	      {
		lines[lineIndex].a.x = GetX(MatchLine, sAn);
		lines[lineIndex].a.y = GetY(MatchLine, sAn);
		lines[lineIndex++].b = MatchLine.b;
	      }
	      else
	      {
		lines[lineIndex].a.x = GetX(MatchLine, sAn);
		lines[lineIndex].a.y = GetY(MatchLine, sAn);
		lines[lineIndex].b.x = GetX(MatchLine, sAnp1);
		lines[lineIndex].b.y = GetY(MatchLine, sAnp1);
		lineIndex++;
	      }
	    }
	  }
        }
      }

//------------------------------------------------------
// if there are no intersections area is zero...
//------------------------------------------------------
      if (lineIndex == 0)
      {
        cleanup( Poly1, Poly2, lines, mask );
        return 0.0;
      }
//------------------------------------------------------
// Remove lines with zero length
//------------------------------------------------------
      for (i = 0; i < lineIndex; i++)
      {
        mask[i] = PointEquals(lines[i].a, lines[i].b) ? 1 : 0;
      }

      r = 0;
      for (p = 0; p < lineIndex; p++)
      {
         if( !mask[p] )
            lines[r++] = lines[p];
      }
      lineIndex = r;
      if (lineIndex == 0)
      {
        cleanup( Poly1, Poly2, lines, mask );
        return 0.0;
      }
//------------------------------------------------------
// Remove line doubles
//------------------------------------------------------
   for (i = 0; i < lineIndex; i++) mask[i] = 0;
   for (i = 0; i < lineIndex; i++)
   {
     if( !mask[i] )
     {
       for (n = i+1; n < lineIndex; n++)
         if( !mask[n] )
           mask[n] = LineEquals(lines[i], lines[n]) ? 1 : 0;
     }
   }
   r = 0;
   for (p = 0; p < lineIndex; p++)
   {
      if( !mask[p] )
         lines[r++] = lines[p];
   }
   lineIndex = r;
//------------------------------------------------------
// if there are no intersections area is zero...
//------------------------------------------------------
   if (lineIndex == 0)
   {
     cleanup( Poly1, Poly2, lines, mask );
     return 0.0;
   }
   else if (lineIndex > 16)
   {
     printf("\nXsection line count = %d\n",lineIndex);
     printf("\nTerminating\n");
     exit( -1 );
   }
//------------------------------------------------------
// find adjacent lines to points...
//------------------------------------------------------
	struct FPOINT points[16];
	double poly_x[16];
	double poly_y[16];
	struct FPOINT matchPnt;
	int pointNumber = 0;
	int vtxcnt = 0;
	for (int i = 0; i < lineIndex; i++)
	{
          int alreadyDone = 0;
	  for (int n = 0; n < pointNumber; n++)
          {
	    if (PointEquals(points[n], lines[i].a))
	      alreadyDone = 1;
	  }
	  if (!alreadyDone)
	    points[pointNumber++] = lines[i].a;

	  alreadyDone = 0;
	  for (int n = 0; n < pointNumber; n++)
	  {
	    if (PointEquals(points[n], lines[i].b))
	      alreadyDone = 1;
	  }
          if (!alreadyDone)
	    points[pointNumber++] = lines[i].b;
	}
//------------------------------------------------------
// remove unconnected line segments...
//------------------------------------------------------
        for( r = 0; r < lineIndex; r++ ) mask[r] = 1;
        for( r = 0; r < lineIndex; r++ )
        {
          for( p = 0; p < lineIndex; p++ )
          {
            if( p != r )
            {
              if( PointEquals(lines[r].a,lines[p].a) ||
                  PointEquals(lines[r].a,lines[p].b) )
                break;
            }
          }
          if( p != lineIndex )
          {
            for( p = 0; p < lineIndex; p++ )
            {
              if( p != r )
              {
                if( PointEquals(lines[r].b,lines[p].a) ||
                    PointEquals(lines[r].b,lines[p].b) )
                  break;
              }
            }
            mask[r] = (p != lineIndex) ? 1 : 0;
          }
          else
            mask[r] = 0;
        }

        vtxcnt = 0;
        for(p = 0; p < pointNumber; p++)
          vtxcnt += mask[p];
        if( vtxcnt < 3 )
        {
          cleanup( Poly1, Poly2, lines, mask );
          return 0.0;
        }
        
//------------------------------------------------------
// setup vertices of intersection polygon
//------------------------------------------------------
        cnt = 2;
        for(p = 0; p < pointNumber; p++)
          if( mask[p] )
          {
            poly_x[0] = lines[p].a.x;
            poly_y[0] = lines[p].a.y;
            poly_x[1] = lines[p].b.x;
            poly_y[1] = lines[p].b.y;
            matchPnt  = lines[p].b;
            mask[p] = 0;
            break;
          }

        while(cnt < vtxcnt)
        {
          for(p = 0; p < pointNumber; p++)
          {
            if( mask[p] )
            {
              if( PointEquals(lines[p].a,matchPnt) )
              { 
                poly_x[cnt] = lines[p].b.x;
                poly_y[cnt++] = lines[p].b.y;
                matchPnt  = lines[p].b;
                mask[p] = 0;
                break;
              }
              else if( PointEquals(lines[p].b,matchPnt) )
              { 
                poly_x[cnt] = lines[p].a.x;
                poly_y[cnt++] = lines[p].a.y;
                matchPnt  = lines[p].a;
                mask[p] = 0;
                break;
              }
            }
          }
        }
        double Area = fabs(poly_area(vtxcnt, poly_x, poly_y));
      
        cleanup( Poly1, Poly2, lines, mask );
        return Area;
}
