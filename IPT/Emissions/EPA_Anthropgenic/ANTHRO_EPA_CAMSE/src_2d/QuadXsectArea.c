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

#define LineTest(one, two) (-((one.b.y - two.b.y)*(two.a.x - two.b.x) - (one.b.x - two.b.x)*(two.a.y - two.b.y))/((one.a.y - one.b.y)*(two.a.x - two.b.x) - (one.a.x - one.b.x)*(two.a.y - two.b.y)))
#define Parallel(one, two) (((one.a.y - one.b.y)*(two.a.x - two.b.x) - (one.a.x - one.b.x)*(two.a.y - two.b.y)) == 0.0)
#define GetX(line, A) (line.a.x * A + line.b.x * (1.0 - A))
#define GetY(line, A) (line.a.y * A + line.b.y * (1.0 - A))
#define DoubleEquals(one, two) (fabs(one - two) < 0.0000000001)
#define LineEquals(one, two) ((DoubleEquals(one.a.x, two.a.x) && DoubleEquals(one.a.y, two.a.y) && DoubleEquals(one.b.x, two.b.x) && DoubleEquals(one.b.y, two.b.y)) || (DoubleEquals(one.a.x, two.b.x) && DoubleEquals(one.a.y, two.b.y) && DoubleEquals(one.b.x, two.a.x) && DoubleEquals(one.b.y, two.a.y)))
#define PointEquals(one, two) (DoubleEquals(one.x, two.x) && DoubleEquals(one.y, two.y))
#define LineConnects(one, two) (PointEquals(one.a, two.a) || PointEquals(one.a, two.b) || PointEquals(one.b, two.a) || PointEquals(one.b, two.b))
#define Length(one, two) (sqrt((one.x - two.x) * (one.x - two.x) + (one.y - two.y) * (one.y - two.y)))

double quadintersectarea_(double *Quad1_x, double *Quad1_y, double *Quad2_x, double *Quad2_y)
{
   // Convert from points to line segments...
   struct LINESEGMENT quad[2][4];
   struct LINESEGMENT shapes[16];
   int shapeIndex = 0;

   for (int i = 0; i < 4; i++)
   {
      // Quad one
      quad[0][i].a.x = Quad1_x[i]; // Point A
      quad[0][i].a.y = Quad1_y[i];
      quad[0][i].b.x = Quad1_x[(i + 1) % 4]; // Point B
      quad[0][i].b.y = Quad1_y[(i + 1) % 4];
      // Quad two
      quad[1][i].a.x = Quad2_x[i]; // Point A
      quad[1][i].a.y = Quad2_y[i];
      quad[1][i].b.x = Quad2_x[(i + 1) % 4]; // Point B
      quad[1][i].b.y = Quad2_y[(i + 1) % 4];
   }

//------------------------------------------------------
// Get the intersection line segments for the quads
//------------------------------------------------------
   for (int quadIndex = 0; quadIndex < 2; quadIndex++)
   {
      for (int i = 0; i < 4; i++)
      {
         struct FPOINT intersect[4];
         double aArray[4];
         int index = 0;
   // Check intersctions for every line in other quad...
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
		  aArray[index] = a;
		  index++;
	       }
	    }
	 }
	 if (index == 0)
	    continue;

   // Remove point doubles...
         for (int r = 0; r < index; r++)
	 {
	    for (int n = 0; n < index; n++)
	    {
	       if (r != n && PointEquals(intersect[r], intersect[n]))
	       {
	          index--;
	          for (int p = n; p < index; p++)
	          {
		     intersect[p] = intersect[p + 1];
		     aArray[p] = aArray[p + 1];
		  }
		  if (r <= n) r--;
		  n--;
	       }
	    }
	 }
   // Sort the lists smallest a value to largest a value...
//       struct FPOINT sortedIntersect[4]; sortedIntersect[0];
	 struct FPOINT sortedIntersect[4];
//       double sortedA[4]; sortedA[0];
	 double sortedA[4];
	 for (int n = 0; n < index; n++)
	 {
	    int inum;
	    for (inum = 0; inum < n; inum++)
	    {
	       if (aArray[n] < sortedA[inum])
	       {
		// Shift list to fit segment in correct place...
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

   // Create line semgents and clip line segments to fit current quad and add to a list...
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
// If there are no intersections area is zero
//------------------------------------------------------
   if (shapeIndex == 0)
     return 0.0;

//------------------------------------------------------
// Remove line doubles and lines with zero length
//------------------------------------------------------
   for (int i = 0; i < shapeIndex; i++)
   {
      for (int n = 0; n < shapeIndex; n++)
      {
         if ((LineEquals(shapes[i], shapes[n]) && (i != n)) || PointEquals(shapes[n].a, shapes[n].b))
         {
	    shapeIndex--;
	    for (int p = n; p < shapeIndex; p++)
	       shapes[p] = shapes[p + 1];
	    if (i <= n) i--;
	    n--;
	 }
      }
   }

//------------------------------------------------------
// Find adjacent lines to points
//------------------------------------------------------
   struct FPOINT points[16];
   struct FPOINT lineadjacency[16];
   int pointNumber = 0;
   for (int i = 0; i < shapeIndex; i++)
   {
      int alreadyDone = 0;
      for (int n = 0; n < pointNumber; n++)
      {
         if (PointEquals(points[n], shapes[i].a))
	    alreadyDone = 1;
      }
      if (!alreadyDone)
      {
         lineadjacency[pointNumber].x = i;
         points[pointNumber] = shapes[i].a;
	 for (int n = 0; n < shapeIndex; n++)
         {
	    if (n != i && (PointEquals(points[pointNumber], shapes[n].a) || PointEquals(points[pointNumber], shapes[n].b)))
	       lineadjacency[pointNumber].y = n;
	 }
	 pointNumber++;
      }
      alreadyDone = 0;
      for (int n = 0; n < pointNumber; n++)
      {
         if (PointEquals(points[n], shapes[i].b))
	    alreadyDone = 1;
      }
      if (!alreadyDone)
      {
         lineadjacency[pointNumber].x = i;
         points[pointNumber] = shapes[i].b;
	 for (int n = 0; n < shapeIndex; n++)
	 {
	    if (n != i && (PointEquals(points[pointNumber], shapes[n].a) || PointEquals(points[pointNumber], shapes[n].b)))
	       lineadjacency[pointNumber].y = n;
	 }
	 pointNumber++;
      }
   }

//------------------------------------------------------
// Find connected points
//------------------------------------------------------
   struct FPOINT adjacency[16];
   int lineUsed[16];
   for( int i = 0; i < 16; i++)
      lineUsed[i] = 0;
      for (int i = 0; i < pointNumber; i++)
      {
         if (!lineUsed[(int)lineadjacency[i].x])
         {
	    for (int n = 0; n < pointNumber; n++)
	    {
	       if (n != i)
	       {
	          if (lineadjacency[i].x == lineadjacency[n].x)
		  {
		     adjacency[i].x = n;
		     adjacency[n].x = i;
		     lineUsed[(int)lineadjacency[i].x] = 1;
		     break;
		  }
	          else if (lineadjacency[i].x == lineadjacency[n].y)
		  {
		     adjacency[i].x = n;
		     adjacency[n].y = i;
		     lineUsed[(int)lineadjacency[i].x] = 1;
		     break;
		  }
	       }
	     }
	  }
	  if (!lineUsed[(int)lineadjacency[i].y])
	  {
	     for (int n = 0; n < pointNumber; n++)
	     {
	        if (n != i)
		{
		   if (lineadjacency[i].y == lineadjacency[n].x)
		   {
		      adjacency[i].y = n;
		      adjacency[n].x = i;
		      lineUsed[(int)lineadjacency[i].y] = 1;
		      break;
		   }
		   else if (lineadjacency[i].y == lineadjacency[n].y)
		   {
		      adjacency[i].y = n;
		      adjacency[n].y = i;
		      lineUsed[(int)lineadjacency[i].y] = 1;
		      break;
		   }

		}
	     }
	  }
       }

//------------------------------------------------------
// Spit up complex figure into smaller triangles
//------------------------------------------------------
       struct FPOINT triangles[16][3];
       int triUsed[16];
       for( int i = 0; i < 16; i++ )
         triUsed[i] = 0;
       int currentTri = 0;
       int pointsRemain = pointNumber;

       while (pointsRemain > 2)
       {
          for (int i = 0; i < pointNumber; i++)
	  {
	    if (!triUsed[i])
	    {
	       int n;
	       for (n = 0; n < 2; n++)
	       {
	          struct LINESEGMENT test;
		  test.a = points[i];
		  if (n == 0)
		     test.b = points[(int)adjacency[i].x];
		  else
		     test.b = points[(int)adjacency[i].y];

		  int less = 0;
		  int greater = 0;
		  for (int p = 0; p < shapeIndex; p++)
		  {
		     if (!Parallel(test, shapes[p]))
		     {
		        double a = LineTest(test, shapes[p]);
			if (a < 0.0)
			   less++;
			if (a > 1.0)
			   greater++;
		     }
		  }
		  if (less % 2 == 1 || greater % 2 == 1)
		  {
		  // Concave section, trainge not in shape so don't use it
		     n = 3;
		     continue;
		  }
	       }
	       if (n == 3)
	          continue;
	       struct LINESEGMENT line;
	       line.a = points[(int)adjacency[i].x];
	       line.b = points[(int)adjacency[i].y];
				
	       double a;
	       for (n = 0; n < shapeIndex; n++)
	       {
	          if (!Parallel(line, shapes[n]))
		     a = LineTest(line, shapes[n]);
		     if (a > 0.0 && a < 1.0)
		        n = shapeIndex + 1;
	       }
	       if (n == shapeIndex + 1)
                  continue; // Traingle intersects shape outline, don't use it

               // Add triagle to list...
	       triangles[currentTri][0] = points[i];
	       triangles[currentTri][1] = points[(int)adjacency[i].x];
	       triangles[currentTri][2] = points[(int)adjacency[i].y];
	       currentTri++;
	       pointsRemain--;
	       triUsed[i] = 1;
				
   // Remove point from list and readjust adjacencies...
	       if (adjacency[(int)adjacency[i].x].x == i)
	       {
		  if (adjacency[(int)adjacency[i].y].x == i)
		  {
		     adjacency[(int)adjacency[i].x].x = adjacency[i].y;
		     adjacency[(int)adjacency[i].y].x = adjacency[i].x;
		  }
		  else
		  {
		     adjacency[(int)adjacency[i].x].x = adjacency[i].y;
		     adjacency[(int)adjacency[i].y].y = adjacency[i].x;
		  }
	       }
	       else
	       {
	          if (adjacency[(int)adjacency[i].y].x == i)
		  {
		     adjacency[(int)adjacency[i].x].y = adjacency[i].y;
		     adjacency[(int)adjacency[i].y].x = adjacency[i].x;
		  }
		  else
		  {
		     adjacency[(int)adjacency[i].x].y = adjacency[i].y;
		     adjacency[(int)adjacency[i].y].y = adjacency[i].x;
		  }
	       }
	       break;
	    }
         }
      }
	
//------------------------------------------------------
// Loop through all of the triangles created
// and add up the area
//------------------------------------------------------
      double area = 0.0; // Sum of the area of the triagles
      for (int i = 0; i < currentTri; i++)
      {
         // Length of the sides of the triangle
         double a = Length(triangles[i][0], triangles[i][1]);
         double b = Length(triangles[i][1], triangles[i][2]);
         double c = Length(triangles[i][0], triangles[i][2]);

         double s = .5 * (a + b + c); // Heron's area of an SSS triagle formula
	 if( s != 0. )
	    area += sqrt(s * (s - a) * (s - b) * (s - c)); // Calculate area and add it to the sum of areas
      }
      return area; // Area found (hurray!)
}
