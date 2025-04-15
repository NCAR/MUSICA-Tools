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
#include <stdbool.h>

#define LineTest(one, two) (-((one.b.y - two.b.y)*(two.a.x - two.b.x) - (one.b.x - two.b.x)*(two.a.y - two.b.y))/((one.a.y - one.b.y)*(two.a.x - two.b.x) - (one.a.x - one.b.x)*(two.a.y - two.b.y)))
#define Parallel(one, two) (((one.a.y - one.b.y)*(two.a.x - two.b.x) - (one.a.x - one.b.x)*(two.a.y - two.b.y)) == 0.0)
#define GetX(line, A) (line.a.x * A + line.b.x * (1.0 - A))
#define GetY(line, A) (line.a.y * A + line.b.y * (1.0 - A))
#define DoubleEquals(one, two) (fabs(one - two) < 0.0000000001)
#define LineEquals(one, two) ((DoubleEquals(one.a.x, two.a.x) && DoubleEquals(one.a.y, two.a.y) && DoubleEquals(one.b.x, two.b.x) && DoubleEquals(one.b.y, two.b.y)) || (DoubleEquals(one.a.x, two.b.x) && DoubleEquals(one.a.y, two.b.y) && DoubleEquals(one.b.x, two.a.x) && DoubleEquals(one.b.y, two.a.y)))
#define PointEquals(one, two) (DoubleEquals(one.x, two.x) && DoubleEquals(one.y, two.y))
#define LineConnects(one, two) (PointEquals(one.a, two.a) || PointEquals(one.a, two.b) || PointEquals(one.b, two.a) || PointEquals(one.b, two.b))
#define Length(one, two) (sqrt((one.x - two.x) * (one.x - two.x) + (one.y - two.y) * (one.y - two.y)))
#define max(a, b) (a>b?a:b)
#define min(a, b) (a<b?a:b)

int pnt_in_triangle(double *Tri_x, double *Tri_y)
{
  double a, b, c;
  int retval;

  retval = 0;
  a = (Tri_x[0] - Tri_x[3])*(Tri_y[1] - Tri_y[3]) - (Tri_x[1] - Tri_x[3])*(Tri_y[0] - Tri_y[3]);
  b = (Tri_x[1] - Tri_x[3])*(Tri_y[2] - Tri_y[3]) - (Tri_x[2] - Tri_x[3])*(Tri_y[1] - Tri_y[3]);
  if( a*b >= 0.0 )
  {
    c = (Tri_x[2] - Tri_x[3])*(Tri_y[0] - Tri_y[3]) - (Tri_x[0] - Tri_x[3])*(Tri_y[2] - Tri_y[3]);
    retval = (b*c >= 0.0) ? 1 : 0;
  }
  return retval;
}

int pnt_in_quad(double *pnt, double *Quad_x, double *Quad_y)
{
  int n;
  double Tri_x[4], Tri_y[4];

  for(n = 0; n < 3; n++)
  {
    Tri_x[n] = Quad_x[n];
    Tri_y[n] = Quad_y[n];
  }
  Tri_x[3] = pnt[0];
  Tri_y[3] = pnt[1];
  n = pnt_in_triangle(Tri_x,Tri_y);
  if( !n )
  {
    for(n = 1; n < 3; n++)
    {
      Tri_x[n] = Quad_x[n+1];
      Tri_y[n] = Quad_y[n+1];
    }
    n = pnt_in_triangle(Tri_x,Tri_y);
    return n;
  }
  else
    return 1;
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

double quadintersectarea_(double *Quad1_x, double *Quad1_y, double *Quad2_x, double *Quad2_y)
{
//-----------------------------------------------------------
// calculate intersection area between two quadrilaterals
//-----------------------------------------------------------
   struct LINESEGMENT quad[2][4];
   struct LINESEGMENT lines[16];
   int lineIndex;
   int i, n, p, r;
   int cnt;
   int mask[16];

//-----------------------------------------------------------
// check for complete inclusion, first Quad1 in Quad2
//-----------------------------------------------------------
   p = 0;
   r = 0;
   for (int i = 0; i < 4; i++)    // Quad one in Quad two
   {
     double vertex[2];
     vertex[0] = Quad1_x[i];
     vertex[1] = Quad1_y[i];
     n = pnt_in_quad(vertex,Quad2_x,Quad2_y);
     if( n == 0 ) 
       break;
     else
       p += n;
   }

//-----------------------------------------------------------
// then Quad2 in Quad1
//-----------------------------------------------------------
   if( p != 4 ) {
     p = 0;
     for (int i = 0; i < 4; i++)    // Quad two in Quad one
     {
       double vertex[2];
       vertex[0] = Quad2_x[i];
       vertex[1] = Quad2_y[i];
       n = pnt_in_quad(vertex,Quad1_x,Quad1_y);
       if( n == 0 ) 
         break;
       else
         p += n;
     }
     if( p == 4 ) r = 1;
   }
   if( p == 4 )
   {
     double area;
     if( r == 0 )
       area = poly_area(p, Quad1_x, Quad1_y);
     else
       area = poly_area(p, Quad2_x, Quad2_y);
     return area;
   }

//-----------------------------------------------------------
// Convert from points to line segments...
//-----------------------------------------------------------
   for (int i = 0; i < 4; i++)
   {
//    Quad one
      quad[0][i].a.x = Quad1_x[i];           //Point A
      quad[0][i].a.y = Quad1_y[i];
      quad[0][i].b.x = Quad1_x[(i + 1) % 4]; //Point B
      quad[0][i].b.y = Quad1_y[(i + 1) % 4];
//    Quad two
      quad[1][i].a.x = Quad2_x[i];           //Point A
      quad[1][i].a.y = Quad2_y[i];
      quad[1][i].b.x = Quad2_x[(i + 1) % 4]; //Point B
      quad[1][i].b.y = Quad2_y[(i + 1) % 4];
   }

   lineIndex = 0;
//-----------------------------------------------------------
// Get the intersection line segments for the quads...
//-----------------------------------------------------------
   for (int quadIndex = 0; quadIndex < 2; quadIndex++)
   {
      struct LINESEGMENT MatchLine;
      struct LINESEGMENT quadTst[4];
      for( int p = 0; p < 4; p++ ) 
        quadTst[p] = quad[1-quadIndex][p];
      for (int i = 0; i < 4; i++)
      {
         struct FPOINT intersect[4];
         double aArray[4];
         int index = 0;
         MatchLine = quad[quadIndex][i];
//-----------------------------------------------------------
// Check intersctions for every line in other quad...
//-----------------------------------------------------------
         for (int n = 0; n < 4; n++)
         {
            if (!Parallel(quadTst[n],MatchLine))
	    {
	       double a = LineTest(quadTst[n], MatchLine);
	       if (a >= 0.0 && a <= 1.0)
               {
		  a = LineTest(MatchLine, quadTst[n]);
		  intersect[index].x = GetX(MatchLine, a);
		  intersect[index].y = GetY(MatchLine, a);
		  aArray[index++] = a;
	       }
	    }
	 }
	 if (index == 0)
	    continue;
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
/*       struct FPOINT sortedIntersect[4]; sortedIntersect[0];
	 double sortedA[4]; sortedA[0]; */
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
      if (lineIndex == 0) return 0.0;
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
      if (lineIndex == 0) return 0.0;
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
   if (lineIndex == 0) return 0.0;
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
        if( vtxcnt < 3 ) return 0.0;
        
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
        return Area;
}
