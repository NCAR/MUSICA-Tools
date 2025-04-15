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

int pnt_in_triangle(double *Tri_x, double *Tri_y)
{
  double a, b, c;

  a = (Tri_x[0] - Tri_x[3])*(Tri_y[1] - Tri_y[3]) - (Tri_x[1] - Tri_x[3])*(Tri_y[0] - Tri_y[3]);
  b = (Tri_x[1] - Tri_x[3])*(Tri_y[2] - Tri_y[3]) - (Tri_x[2] - Tri_x[3])*(Tri_y[1] - Tri_y[3]);
  c = (Tri_x[2] - Tri_x[3])*(Tri_y[0] - Tri_y[3]) - (Tri_x[0] - Tri_x[3])*(Tri_y[2] - Tri_y[3]);
  if( a*b >= 0.0 && b*c >= 0.0 )
    return 1;
  else
    return 0;
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

void set_alphabeta(double *Quad_x, double *Quad_y, double *alpha, double *beta)
{

   alpha[0] = Quad_x[0];
   alpha[1] = Quad_x[1] - Quad_x[0];
   alpha[2] = Quad_x[3] - Quad_x[0];
   alpha[3] = (Quad_x[0] + Quad_x[2]) - (Quad_x[1] + Quad_x[3]);

   beta[0] = Quad_y[0];
   beta[1] = Quad_y[1] - Quad_y[0];
   beta[2] = Quad_y[3] - Quad_y[0];
   beta[3] = (Quad_y[0] + Quad_y[2]) - (Quad_y[1] + Quad_y[3]);

}

int xy_n_quad(double *Quad_x, double *Quad_y, double *xy)
{
   double aa, bb, cc;
   double l, m;
   double alpha[4], beta[4];
   int    unity_map;

   unity_map = (Quad_x[0] == 0. && Quad_x[1] == 1. && Quad_x[2] == 1. && Quad_x[3] == 0.);
   if( unity_map )
     unity_map = (Quad_y[0] == 0. && Quad_y[1] == 0. && Quad_y[2] == 1. && Quad_y[3] == 1.);

   if( !unity_map ) {
     set_alphabeta(Quad_x, Quad_y, alpha, beta);

     aa = alpha[3]*beta[2] - alpha[2]*beta[3];
     bb = alpha[3]*beta[0] + alpha[1]*beta[2] + xy[0]*beta[3];
     bb -= (alpha[0]*beta[3] + alpha[2]*beta[1] + xy[1]*alpha[3]);
     cc = alpha[1]*beta[0] + xy[0]*beta[1] - (alpha[0]*beta[1] + xy[1]*alpha[1]);

     m = .5*(sqrt(bb*bb - 4.*aa*cc) -bb)/aa;
     if( m < 0.0 || m > 1.0 ) 
       return 0;
     else 
     {
        l = (xy[0] - (alpha[0] + alpha[2]*m))/(alpha[1] + alpha[3]*m);
        if( l < 0.0 || l > 1.0 ) 
          return 0;
        else
          return 1;
     }
   } else {
     if( xy[1] < 0.0 || xy[1] > 1.0 ) 
       return 0;
     else 
       if( xy[0] < 0.0 || xy[0] > 1.0 )
         return 0;
       else
         return 1;
   }
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
// Convert from points to line segments...
//-----------------------------------------------------------
   struct LINESEGMENT quad[2][4];
   struct LINESEGMENT shapes[16];
   int shapeIndex = 0;
   int i, n, p, r;
   int cnt;
   int mask[16];
   bool gotit;

//-----------------------------------------------------------
// first a few diagnostics
//-----------------------------------------------------------
   {
     double area;
     area = poly_area(4, Quad1_x, Quad1_y);
     area = poly_area(4, Quad2_x, Quad2_y);
   }

//-----------------------------------------------------------
// check for complete inclusion
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

   for (int i = 0; i < 4; i++)    // Quad one
   {
      quad[0][i].a.x = Quad1_x[i]; //Point A
      quad[0][i].a.y = Quad1_y[i];
      quad[0][i].b.x = Quad1_x[(i + 1) % 4]; //Point B
      quad[0][i].b.y = Quad1_y[(i + 1) % 4];
      quad[1][i].a.x = Quad2_x[i]; //Point A
      quad[1][i].a.y = Quad2_y[i];
      quad[1][i].b.x = Quad2_x[(i + 1) % 4]; //Point B
      quad[1][i].b.y = Quad2_y[(i + 1) % 4];
   }


//-----------------------------------------------------------
// Get the intersection line segments for the quads...
//-----------------------------------------------------------
   for (int quadIndex = 0; quadIndex < 2; quadIndex++)
   {
      for (int i = 0; i < 4; i++)
      {
         struct FPOINT intersect[4];
         double aArray[4];
         int index = 0;
//-----------------------------------------------------------
// Check intersctions for every line in other quad...
//-----------------------------------------------------------
         for (int n = 0; n < 4; n++)
         {
            if (!Parallel(quad[1 - quadIndex][n], quad[quadIndex][i]))
	    {
	       double a = LineTest(quad[1 - quadIndex][n], quad[quadIndex][i]);
	       if (a >= 0.0 && a <= 1.0)
               {
		  a = LineTest(quad[quadIndex][i], quad[1 - quadIndex][n]);
		  intersect[index].x = GetX(quad[quadIndex][i], a);
		  intersect[index].y = GetY(quad[quadIndex][i], a);
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
	      mask[r] = (!mask[r] && PointEquals(intersect[r], intersect[p])) ? 1 : 0;
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
// create line semgents and clip line segments to fit current quad and add to a list...
//-----------------------------------------------------------
	 for (int n = 0; n < index; n+=2)
	 {
	   if (sortedA[n] < 0.0)
	   {
	      if (sortedA[n + 1] > 1.0)
	      {
	        shapes[shapeIndex] = quad[quadIndex][i];
		shapeIndex++;
	      }
	      else if (sortedA[n + 1] >= 0.0)
	      {
	        shapes[shapeIndex].a.x = GetX(quad[quadIndex][i], sortedA[n + 1]);
		shapes[shapeIndex].a.y = GetY(quad[quadIndex][i], sortedA[n + 1]);
		shapes[shapeIndex].b = quad[quadIndex][i].b;
		shapeIndex++;
	      }
	    }
	    else if (sortedA[n] > 1.0)
	    {
	      if (sortedA[n + 1] <= 1.0 && sortedA[n + 1] >= 0.0)
	      {
	        shapes[shapeIndex].b.x = GetX(quad[quadIndex][i], sortedA[n + 1]);
		shapes[shapeIndex].b.y = GetY(quad[quadIndex][i], sortedA[n + 1]);
		shapes[shapeIndex].a = quad[quadIndex][i].a;
		shapeIndex++;
	      }
	    }
	    else
	    {
	      if (sortedA[n + 1] > 1.0)
	      {
		shapes[shapeIndex].b.x = GetX(quad[quadIndex][i], sortedA[n]);
		shapes[shapeIndex].b.y = GetY(quad[quadIndex][i], sortedA[n]);
		shapes[shapeIndex].a = quad[quadIndex][i].a;
		shapeIndex++;
	      }
	      else if (sortedA[n + 1] < 0.0)
	      {
		shapes[shapeIndex].a.x = GetX(quad[quadIndex][i], sortedA[n]);
		shapes[shapeIndex].a.y = GetY(quad[quadIndex][i], sortedA[n]);
		shapes[shapeIndex].b = quad[quadIndex][i].b;
		shapeIndex++;
	      }
	      else
	      {
		shapes[shapeIndex].a.x = GetX(quad[quadIndex][i], sortedA[n]);
		shapes[shapeIndex].a.y = GetY(quad[quadIndex][i], sortedA[n]);
		shapes[shapeIndex].b.x = GetX(quad[quadIndex][i], sortedA[n + 1]);
		shapes[shapeIndex].b.y = GetY(quad[quadIndex][i], sortedA[n + 1]);
		shapeIndex++;
	      }
	    }
	  }
        }
      }

//------------------------------------------------------
// if there are no intersections area is zero...
//------------------------------------------------------
      if (shapeIndex == 0) return 0.0;
//------------------------------------------------------
// Remove lines with zero length
//------------------------------------------------------
      for (i = 0; i < shapeIndex; i++) mask[i] = 0;
        mask[i] = PointEquals(shapes[i].a, shapes[i].b) ? 1 : 0;

      r = 0;
      for (p = 0; p < shapeIndex; p++)
         if( !mask[p] )
            shapes[r++] = shapes[p];
      shapeIndex = r;
      if (shapeIndex == 0) return 0.0;
//------------------------------------------------------
// Remove line doubles
//------------------------------------------------------
   for (i = 0; i < shapeIndex; i++) mask[i] = 0;
   for (i = 0; i < shapeIndex; i++)
     if( !mask[i] )
     {
       for (n = i+1; n < shapeIndex; n++)
         mask[n] = (!mask[n] && LineEquals(shapes[i], shapes[n])) ? 1 : 0;
     }
   r = 1;
   for (p = 1; p < shapeIndex; p++)
      if( !mask[p] )
         shapes[r++] = shapes[p];
   shapeIndex = r;
//------------------------------------------------------
// if there are no intersections area is zero...
//------------------------------------------------------
   if (shapeIndex == 0) return 0.0;
//------------------------------------------------------
// find adjacent lines to points...
//------------------------------------------------------
	struct FPOINT points[16];
	double poly_x[16];
	double poly_y[16];
	struct FPOINT matchPnt;
	int pointNumber = 0;
	int vtxcnt = 0;
	for (int i = 0; i < shapeIndex; i++)
	{
          int alreadyDone = 0;
	  for (int n = 0; n < pointNumber; n++)
          {
	    if (PointEquals(points[n], shapes[i].a))
	      alreadyDone = 1;
	  }
	  if (!alreadyDone)
	    points[pointNumber++] = shapes[i].a;

	  alreadyDone = 0;
	  for (int n = 0; n < pointNumber; n++)
	  {
	    if (PointEquals(points[n], shapes[i].b))
	      alreadyDone = 1;
	  }
          if (!alreadyDone)
	    points[pointNumber++] = shapes[i].b;
	}
//------------------------------------------------------
// remove unconnected line segments...
//------------------------------------------------------
        for( r = 0; r < shapeIndex; r++ ) mask[r] = 1;
        for( r = 0; r < shapeIndex; r++ )
        {
          for( p = 0; p < shapeIndex; p++ )
          {
            if( p != r )
            {
              if( PointEquals(shapes[r].a,shapes[p].a) ||
                  PointEquals(shapes[r].a,shapes[p].b) )
                break;
            }
          }
          if( p != shapeIndex )
          {
            for( p = 0; p < shapeIndex; p++ )
            {
              if( p != r )
              {
                if( PointEquals(shapes[r].b,shapes[p].a) ||
                    PointEquals(shapes[r].b,shapes[p].b) )
                  break;
              }
            }
            mask[r] = (p != shapeIndex) ? 1 : 0;
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
            poly_x[0] = shapes[p].a.x;
            poly_y[0] = shapes[p].a.y;
            poly_x[1] = shapes[p].b.x;
            poly_y[1] = shapes[p].b.y;
            matchPnt  = shapes[p].b;
            mask[p] = 0;
            break;
          }

        while(cnt < vtxcnt)
        {
          for(p = 0; p < pointNumber; p++)
          {
            if( mask[p] )
            {
              if( PointEquals(shapes[p].a,matchPnt) )
              { 
                poly_x[cnt] = shapes[p].b.x;
                poly_y[cnt++] = shapes[p].b.y;
                matchPnt  = shapes[p].b;
                mask[p] = 0;
                break;
              }
              else if( PointEquals(shapes[p].b,matchPnt) )
              { 
                poly_x[cnt] = shapes[p].a.x;
                poly_y[cnt++] = shapes[p].a.y;
                matchPnt  = shapes[p].a;
                mask[p] = 0;
                break;
              }
            }
          }
        }
        double Area = fabs(poly_area(vtxcnt, poly_x, poly_y));
        return Area;
}
