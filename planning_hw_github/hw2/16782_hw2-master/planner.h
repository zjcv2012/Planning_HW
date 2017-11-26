#include <math.h>
#include "mex.h"
#include <unordered_map>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <time.h>
#include <unordered_set>
#include <set>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654
#define Epsilon 1.2

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

using namespace std;


class RRTTree {
public:
	RRTTree(int num, double* armstart_anglesV_rad)
	{
		numofDOFs = num;
		vertices.clear();
		edges.clear();
		vertices.push_back(armstart_anglesV_rad);
	}

	int AddVertex(double* sample_angelsV_rad)
	{
		int vid = vertices.size();
		vertices.push_back(sample_angelsV_rad);

		return vid;
	}

	int GetRootId()
	{
		return 0;
	}

	int GetLastId()
	{
		return vertices.size()-1;
	}

	void AddEdge(int sid, int eid)
	{
		edges[eid] = sid;
	}

	int GetNearestNeibour(double* sample_angelsV_rad)
	{
		int nearest_id = 0;
		double nearest_dist = numeric_limits<double>::infinity();
		for(int i=0; i<vertices.size(); i++)
		{
			int tmp_dist = 0;
			for(int j=0; j<numofDOFs; j++)
			{
				tmp_dist+=(vertices[i][j]-sample_angelsV_rad[j])*(vertices[i][j]-sample_angelsV_rad[j]);
			}
			tmp_dist=sqrt(tmp_dist);
			if(tmp_dist<nearest_dist)
			{
				nearest_dist = tmp_dist;
				nearest_id = i;
			}
		}

		return nearest_id;
	}

	double* GetVertice(int id)
	{
		return vertices[id];
	}

	int numofVertices()
	{
		return vertices.size();
	}

	int GetParent(int id)
	{
		return edges[id];
	}


private:
	vector<double* > vertices;
	unordered_map<int,int> edges;
	int numofDOFs;
};


typedef struct {
  int X1, Y1;
  int X2, Y2;
  int Increment;
  int UsingYIndex;
  int DeltaX, DeltaY;
  int DTerm;
  int IncrE, IncrNE;
  int XIndex, YIndex;
  int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
  params->UsingYIndex = 0;

  if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
    (params->UsingYIndex)++;

  if (params->UsingYIndex)
    {
      params->Y1=p1x;
      params->X1=p1y;
      params->Y2=p2x;
      params->X2=p2y;
    }
  else
    {
      params->X1=p1x;
      params->Y1=p1y;
      params->X2=p2x;
      params->Y2=p2y;
    }

   if ((p2x - p1x) * (p2y - p1y) < 0)
    {
      params->Flipped = 1;
      params->Y1 = -params->Y1;
      params->Y2 = -params->Y2;
    }
  else
    params->Flipped = 0;

  if (params->X2 > params->X1)
    params->Increment = 1;
  else
    params->Increment = -1;

  params->DeltaX=params->X2-params->X1;
  params->DeltaY=params->Y2-params->Y1;

  params->IncrE=2*params->DeltaY*params->Increment;
  params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
  params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

  params->XIndex = params->X1;
  params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
  if (params->UsingYIndex)
    {
      *y = params->XIndex;
      *x = params->YIndex;
      if (params->Flipped)
        *x = -*x;
    }
  else
    {
      *x = params->XIndex;
      *y = params->YIndex;
      if (params->Flipped)
        *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
  if (params->XIndex == params->X2)
    {
      return 0;
    }
  params->XIndex += params->Increment;
  if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
    params->DTerm += params->IncrE;
  else
    {
      params->DTerm += params->IncrNE;
      params->YIndex += params->Increment;
    }
  return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
		   int x_size,
 		   int y_size)

{
	bresenham_param_t params;
	int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
    
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
            return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
		   int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
    
 	//iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
    y1 = 0;
	for(i = 0; i < numofDOFs; i++)
	{
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
				return 0;
	}    
}

void GenerateRandomConfig(double* random_angles, double *map, int numofDOFs, int x_size, int y_size)
{
	while(1)
	{
		for(int i=0; i<numofDOFs; i++)
		{
			int irand = rand() % 360;
			random_angles[i]= (double)irand / 360 *2*PI;
		}
		if(IsValidArmConfiguration(random_angles, numofDOFs, map, x_size, y_size))
		{
			break;
		}
	}
}

double dist_angles(double* angle1, double* angle2, int numofDOFs)
{
	double dist = 0;
	for(int i=0; i<numofDOFs; i++)
	{
		dist = dist + (angle1[i]-angle2[i]) * (angle1[i]-angle2[i]);
	}
	return sqrt(dist);
}

double* Extend(double* start_angles_rad, double* end_angles_rad, double* map, int numofDOFs, int x_size, int y_size, double ep)
{
	double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(start_angles_rad[j] - end_angles_rad[j]))
            distance = fabs(start_angles_rad[j] - end_angles_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));

    double* tmp_angles_rad = (double*)malloc(numofDOFs*sizeof(double));
    double* saved_angles_rad = NULL;

    for (i = 1; i < numofsamples; i++)
    {
    	for(j = 0; j<numofDOFs; j++)
    	{
    		tmp_angles_rad[j] = start_angles_rad[j] + ((double)(i)/(numofsamples-1))*(end_angles_rad[j] - start_angles_rad[j]);
    	}
    	if(IsValidArmConfiguration(tmp_angles_rad, numofDOFs, map, x_size, y_size) && 
    		dist_angles(tmp_angles_rad, start_angles_rad, numofDOFs) < ep)
    	{
    		if(i==1)
    		{
    			saved_angles_rad = (double*)malloc(numofDOFs*sizeof(double));
    		}
    		memcpy(saved_angles_rad, tmp_angles_rad,numofDOFs*sizeof(double));
    	}
    	else
    	{break;}
    }

    free(tmp_angles_rad);
    return saved_angles_rad;
}



bool angles_equal(double* angle1, double* angle2, int numofDOFs)
{
	if(angle1 == NULL || angle2 == NULL)
		return false;
	for(int i=0; i<numofDOFs; i++)
	{
		if(angle1[i] != angle2[i])
		{
			return false;
			break;
		}
	}
	return true;
}


class PRMGraph{
public:
	PRMGraph(int num, double* armstart_anglesV_rad, double* armgoal_anglesV_rad)
	{
		numofDOFs = num;
		vertices.clear();
		edges.clear();
		g.clear();
		//bool(PRMGraph::* fn_pt)(int,int) = compareVertice;
		//set<int,bool(PRMGraph::*)(int,int)> CloseSet (fn_pt);

		int start_id = AddVertex(armstart_anglesV_rad);
		start_angles_rad = armstart_anglesV_rad;
		goal_angles_rad = armgoal_anglesV_rad;
	}

	int AddVertex(double* sample_angelsV_rad)
	{
		int vid = vertices.size();
		vertices.push_back(sample_angelsV_rad);
		vector<int> tmp;
		edges.push_back(tmp);

		return vid;
	}

	void AddEdge(int sid, int eid)
	{
		edges[sid].push_back(eid);
		edges[eid].push_back(sid);
	}

	int numofVertices()
	{
		return vertices.size();
	}

	double* GetVertice(int id)
	{
		return vertices[id];
	}

	vector<int> GetSuccessors(int id)
	{
		return edges[id];
	}



	vector<double*> CopyVertices()
	{
		vector<double*> copy;
		double* tmp = (double*)malloc(numofDOFs*sizeof(double));
		for(auto i : vertices)
		{
			memcpy(tmp, i,numofDOFs*sizeof(double));
			copy.push_back(tmp);
		}
		free(tmp);

		return copy;

	}

	unordered_set<int> GetKNN(double* sample_angelsV_rad, int k)
	{
		unordered_set<int> KNN;

		int n = min(k, (int)vertices.size());

/*		if(vertices.size()>5)
		{
			mexPrintf("test: %f, %f, %f \n",sample_angelsV_rad[0], vertices[2][0], vertices[3][0]);
		}*/


		while(n>0)
		{
			int nearest_id = 0;
			double nearest_dist = numeric_limits<double>::infinity();
			for(int i=0; i<vertices.size()-1; i++)
			{
				if(KNN.count(i)==1) continue;
				double tmp_dist = 0;
				for(int j=0; j<numofDOFs; j++)
				{
					//mexPrintf("%f\n",vertices[i][j]-sample_angelsV_rad[j]);
					tmp_dist = tmp_dist + (vertices[i][j]-sample_angelsV_rad[j])*(vertices[i][j]-sample_angelsV_rad[j]);
				}
				tmp_dist=sqrt(tmp_dist);
				if(tmp_dist<nearest_dist)
				{
					nearest_dist = tmp_dist;
					nearest_id = i;
					
				}
			}
			KNN.insert(nearest_id);	
			n--;			
		}


		return KNN;
	}

	void Preprocess(double* map, double* armgoal_anglesV_rad,int N, int x_size, int y_size) //N is the number of sampling, k is the maximum number of connected neighbours.
	{
		
		int i = 0;
		int goal_id = AddVertex(armgoal_anglesV_rad);

		int num_edges = 0;

		double* random_angles_rad =(double*)malloc(numofDOFs*sizeof(double));

		while(i<N)
		{
			i++;
			GenerateRandomConfig(random_angles_rad, map, numofDOFs, x_size, y_size);

			double* tmp_angles_rad = (double*)malloc(numofDOFs*sizeof(double));
			memcpy(tmp_angles_rad, random_angles_rad, numofDOFs*sizeof(double));

			int random_id = AddVertex(tmp_angles_rad);

		//	mexPrintf("%d\n",random_id);
			unordered_set<int> nb_id = GetKNN(tmp_angles_rad, 5);
		/*	for(unordered_set<int>::iterator itr = nb_id.begin(); itr!=nb_id.end(); itr++)
			{
				mexPrintf("%d ",*itr);
			}
			mexPrintf("\n");*/

			for(unordered_set<int>::iterator itr = nb_id.begin(); itr!=nb_id.end(); itr++)
			{
				double* extend_angles_rad = Extend(vertices[*itr], tmp_angles_rad, map, numofDOFs, x_size, y_size, numeric_limits<double>::infinity());

				if(angles_equal(extend_angles_rad, tmp_angles_rad, numofDOFs))
				{
					AddEdge(random_id,*itr);
					num_edges++;
				}

			}


		}


		mexPrintf("Edges connected to goal: %d\n", edges[goal_id].size());
		mexPrintf("Edges connected to start: %d\n", edges[0].size());
		mexPrintf("Number of edges: %d\n", num_edges);
	}

	void PreprocessR(double* map, double* armgoal_anglesV_rad,int N, double R, int k, int x_size, int y_size) //N is the number of sampling, k is the maximum number of connected neighbours.
	{
		
		int i = 0;
		int goal_id = AddVertex(armgoal_anglesV_rad);

		int num_edges = 0;

		while(i<N)
		{

			i++;
			int num_connected = 0;

			double* random_angles_rad =(double*)malloc(numofDOFs*sizeof(double));


			GenerateRandomConfig(random_angles_rad, map, numofDOFs, x_size, y_size);

			double* tmp_angles_rad = (double*)malloc(numofDOFs*sizeof(double));
			memcpy(tmp_angles_rad, random_angles_rad, numofDOFs*sizeof(double));

			int random_id = AddVertex(tmp_angles_rad);

			for(int j=0; j<vertices.size()-1; j++)
			{
				if(dist_angles(vertices[j],	tmp_angles_rad, numofDOFs)<R)
				{
					double* extend_angles_rad = Extend(vertices[j], tmp_angles_rad, map, numofDOFs, x_size, y_size, R);
					if(angles_equal(extend_angles_rad, tmp_angles_rad, numofDOFs))
					{
						AddEdge(j, random_id);
						num_connected++;
						num_edges++;
						if(num_connected>=k) break;
					}
				}
			}
		}	

		mexPrintf("Edges connected to goal: %d\n", edges[1].size());
		mexPrintf("Edges connected to start: %d\n", edges[0].size());
		mexPrintf("Number of edges: %d\n", num_edges);
	}

	bool compareVertice(const int l_id, const int r_id)
	{
	   double l_dist = g[l_id] + 3*dist_angles(vertices[l_id], goal_angles_rad, numofDOFs);
	   double r_dist = g[r_id] + 3*dist_angles(vertices[r_id], goal_angles_rad, numofDOFs);

	   return l_dist < r_dist;
	}

	unordered_map<int, double> g;
	//set<int> OpenSet;

private:
	vector<double* > vertices;
	vector<vector<int>> edges;
	int numofDOFs;
	double* start_angles_rad;
	double* goal_angles_rad;	
};