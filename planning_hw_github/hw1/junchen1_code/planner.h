#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>


//Add for planning algorithms
#include <vector> 
#include <unordered_set>
#include <unordered_map>
#include <utility> 
#include <queue>


/* Input Arguments */
#define MAP_IN      prhs[0]
#define ROBOT_IN    prhs[1]
#define GOAL_IN     prhs[2]


/* Output Arguments */
#define ACTION_OUT  plhs[0]

/*access to the map is shifted to account for 0-based indexing in the map, whereas
1-based indexing in matlab (so, robotpose and goalpose are 1-indexed)*/
#define GETMAPINDEX(X, Y, XSIZE, YSIZE) ((Y-1)*XSIZE + (X-1))

#if !defined(MAX)
#define MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

/* Primitives Information */
#define NUMOFDIRS 8
#define NUMOFPRIMS 5
#define NUMOFINTERSTATES 10
#define NUMOFDIM 3

#define RES 0.1
#define thresh 0.5
