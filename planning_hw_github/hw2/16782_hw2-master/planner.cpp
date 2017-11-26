/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include "planner.h"
#include <time.h>

static void plannerPRM( 
        double*  map,
        int x_size,
        int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
     int numofDOFs,
     double*** plan,
     int* planlength)
{
  *plan = NULL;
  *planlength = 0;
  srand (time(NULL)); //seed the time

  PRMGraph graph(numofDOFs, armstart_anglesV_rad, armgoal_anglesV_rad);
  //graph.Preprocess(map, armgoal_anglesV_rad, 10000, x_size, y_size);

  clock_t t1=clock();
  graph.PreprocessR(map, armgoal_anglesV_rad, 30000, Epsilon, 15, x_size, y_size);
  mexPrintf("Constructing graph takes %f\n", (float)(clock()-t1)/CLOCKS_PER_SEC);
 
  unordered_set<int> CloseSet;
  unordered_set<int> OpenSet;
  unordered_map<int, int> parent;

  bool found = false;
  OpenSet.insert(0);
  graph.g[0] = 0;

  int goal_id = 1;

  int num_expanded = 0;

  while(!OpenSet.empty() && !found)
  {
    unordered_set<int>::iterator itr;
    unordered_set<int>::iterator close_index;
    double min_f =numeric_limits<double>::infinity();

    num_expanded++;

    for(itr=OpenSet.begin(); itr!=OpenSet.end(); ++itr)
    {
      double f = graph.g[*itr]+ 3*dist_angles(graph.GetVertice(*itr), graph.GetVertice(goal_id), numofDOFs);
      if (f<min_f)
      {
          min_f = f;
          close_index = itr;
      }
    }
    int current_id = *close_index;
    //mexPrintf("test2: %d\n", current_id);

    CloseSet.insert(current_id);
    OpenSet.erase(close_index);

    vector<int> successors = graph.GetSuccessors(current_id);
    for(auto id :successors)
    {
      if(id == goal_id)
      {
        mexPrintf("Found! \n");
        mexPrintf("number of expanded node: %d\n", num_expanded);
        found = true;
        parent[goal_id] = current_id;
        break;
      }

      if(CloseSet.count(id) == 0)
      {
        if(OpenSet.count(id) == 0)
        {
          OpenSet.insert(id);
          graph.g[id] = graph.g[current_id] + dist_angles(graph.GetVertice(id), graph.GetVertice(current_id), numofDOFs);
          parent[id] = current_id;
        }
        else
        {
          if(graph.g[id]> graph.g[current_id] + dist_angles(graph.GetVertice(id), graph.GetVertice(current_id), numofDOFs))
          {
            graph.g[id] = graph.g[current_id] + dist_angles(graph.GetVertice(id), graph.GetVertice(current_id), numofDOFs);
            parent[id] = current_id;
          }
        }
      }
    }

  }


  if(found)
  {
    vector<int> path;
    int next_id = goal_id;
    while(next_id!=0)
    {
      path.insert(path.begin(),next_id);
      next_id = parent[next_id];
    }
    path.insert(path.begin(),0);

    *planlength = path.size();

    *plan = (double**) malloc(path.size()*sizeof(double*));

    for(int i=0; i<path.size(); i++)
    {
      (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
      memcpy((*plan)[i], graph.GetVertice(path[i]), numofDOFs*sizeof(double));
    }
  }

  else
  {
    mexPrintf("Path Not Found!\n");
  }
  return;
}


static void plannerRRT(
       double*  map,
       int x_size,
       int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
     int numofDOFs,
     double*** plan,
     int* planlength)
{
  //no plan by default
  *plan = NULL;
  *planlength = 0;

  int num_samples = 0;

  srand (time(NULL)); //seed the time

  RRTTree tree(numofDOFs,armstart_anglesV_rad);

  bool found = false; //path is not found in default

  double* random_angles_rad =(double*)malloc(numofDOFs*sizeof(double));

  while(!found && num_samples < 50000)
  {
    num_samples++;
    GenerateRandomConfig(random_angles_rad, map, numofDOFs, x_size, y_size);

    int nb_id = tree.GetNearestNeibour(random_angles_rad);
    double* extend_angles_rad = Extend(tree.GetVertice(nb_id), random_angles_rad, map, numofDOFs, x_size, y_size, Epsilon);
    if(extend_angles_rad != NULL) // if it is not NULL
    {
      //mexPrintf("Test2\n");
      int extend_id = tree.AddVertex(extend_angles_rad);
      tree.AddEdge(nb_id, extend_id);

      extend_angles_rad = Extend(extend_angles_rad,armgoal_anglesV_rad,map, numofDOFs, x_size, y_size, Epsilon);
      if(angles_equal(extend_angles_rad, armgoal_anglesV_rad,numofDOFs))
      {
        int goal_id = tree.AddVertex(armgoal_anglesV_rad);
        tree.AddEdge(extend_id,goal_id);
        found = true;
        mexPrintf("Caught\n");
        break;
      }

    }
  }

  if(found)
  {
    vector<int> path;
    int next_id = tree.GetLastId();
    while(next_id!=0)
    {
      path.insert(path.begin(),next_id);
      next_id = tree.GetParent(next_id);
    }
    path.insert(path.begin(),tree.GetRootId());

    *planlength = path.size();

    *plan = (double**) malloc(path.size()*sizeof(double*));

    for(int i=0; i<path.size(); i++)
    {
      (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
      memcpy((*plan)[i], tree.GetVertice(path[i]), numofDOFs*sizeof(double));
    }
  }

  else
  {
    mexPrintf("Path Not Found!\n");
  }

  mexPrintf("Number of Samples: %d\n", num_samples);
  return;

}


static void plannerRRTSTAR(
       double*  map,
       int x_size,
       int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
     int numofDOFs,
     double*** plan,
     int* planlength)
{
  //no plan by default
  *plan = NULL;
  *planlength = 0;

  double radius = 2*Epsilon;

  int num_samples = 0;

  srand (time(NULL)); //seed the time

  RRTTree tree(numofDOFs,armstart_anglesV_rad);

  bool found = false; //path is not found in default

  double* random_angles_rad =(double*)malloc(numofDOFs*sizeof(double));

  unordered_map<int, double> cost;
  cost[0] = 0;
  clock_t t1;

  int goal_id;

  //while(!found && num_samples < 50000)
  while(num_samples<50000)
  {
    if(found)
    {
      if((float)(clock()-t1)/CLOCKS_PER_SEC > 1)
        break;
    }
    num_samples++;
    //mexPrintf("test: %d\n", num_samples);

    GenerateRandomConfig(random_angles_rad, map, numofDOFs, x_size, y_size);

    int nb_id = tree.GetNearestNeibour(random_angles_rad);
    double* extend_angles_rad = Extend(tree.GetVertice(nb_id), random_angles_rad, map, numofDOFs, x_size, y_size, Epsilon);
    if(extend_angles_rad != NULL) // if it is not NULL
    {
      int extend_id;
      
      //Try to improve the path
      if(angles_equal(extend_angles_rad, random_angles_rad, numofDOFs))
      {
        extend_id = tree.AddVertex(extend_angles_rad);
        vector<int> extend_nb;
        double min_cost = cost[nb_id] +  dist_angles(tree.GetVertice(nb_id), tree.GetVertice(extend_id), numofDOFs);
        int parent = nb_id;
        for(int i=0; i<tree.numofVertices()-1; i++) // 
        {
          if(dist_angles(tree.GetVertice(extend_id),tree.GetVertice(i), numofDOFs) < radius)
          {
            extend_nb.push_back(i);
            if(cost[i] + dist_angles(tree.GetVertice(extend_id),tree.GetVertice(i), numofDOFs) < min_cost)
            {
              min_cost = cost[i] + dist_angles(tree.GetVertice(extend_id),tree.GetVertice(i), numofDOFs);
              parent = i; 
            }
          }
        }
        cost[extend_id] = min_cost;
        tree.AddEdge(parent, extend_id);

        for(auto nb: extend_nb)
        {
          if(cost[nb] > cost[extend_id] + dist_angles(tree.GetVertice(extend_id),tree.GetVertice(nb), numofDOFs))
          {
            cost[nb] = cost[extend_id] + dist_angles(tree.GetVertice(extend_id),tree.GetVertice(nb), numofDOFs);
            tree.AddEdge(extend_id, nb);
          }
        }
      }
      else
      {
        extend_id = tree.AddVertex(extend_angles_rad);
        tree.AddEdge(nb_id, extend_id);
        cost[extend_id] = cost[nb_id] + dist_angles(tree.GetVertice(nb_id), tree.GetVertice(extend_id), numofDOFs);
      }

      if(!found)
      {
        extend_angles_rad = Extend(extend_angles_rad,armgoal_anglesV_rad,map, numofDOFs, x_size, y_size, Epsilon);
        if(angles_equal(extend_angles_rad, armgoal_anglesV_rad,numofDOFs))
        {
          goal_id = tree.AddVertex(armgoal_anglesV_rad);
          tree.AddEdge(extend_id,goal_id);
          cost[goal_id] = cost[extend_id] + dist_angles(tree.GetVertice(extend_id), tree.GetVertice(goal_id),numofDOFs);
          mexPrintf("Caught\n");
          t1=clock();
          found = true;

          //break;
        }
      }

    }

  }

  if(found)
  {
    vector<int> path;
    int next_id = goal_id;
    while(next_id!=0)
    {
     // mexPrintf("test2: %d\n", next_id);
      path.insert(path.begin(),next_id);
      next_id = tree.GetParent(next_id);
    }
    path.insert(path.begin(),tree.GetRootId());

    *planlength = path.size();

    *plan = (double**) malloc(path.size()*sizeof(double*));

    for(int i=0; i<path.size(); i++)
    {
      (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
      memcpy((*plan)[i], tree.GetVertice(path[i]), numofDOFs*sizeof(double));
    }
  }

  else
  {
    mexPrintf("Path Not Found!\n");
  }

  mexPrintf("Number of Samples: %d\n", num_samples);
  return;

}

static void plannerRRTCONNECT(
       double*  map,
       int x_size,
       int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
     int numofDOFs,
     double*** plan,
     int* planlength)
{
  //no plan by default
  *plan = NULL;
  *planlength = 0;

  int num_samples = 0;

  srand (time(NULL)); //seed the time

  RRTTree start_tree(numofDOFs,armstart_anglesV_rad);
  RRTTree goal_tree(numofDOFs, armgoal_anglesV_rad);

  bool found = false; //path is not found in default

  double* random_angles_rad =(double*)malloc(numofDOFs*sizeof(double));

  bool extend_stree = true; //boolean for switch


  while(!found && num_samples < 50000)
  {
    num_samples++;
    if(extend_stree)
    {
      //Extend start_tree to random_angles
      GenerateRandomConfig(random_angles_rad, map, numofDOFs, x_size, y_size);
      int start_nb_id = start_tree.GetNearestNeibour(random_angles_rad);
      double* start_extend_angles_rad = Extend(start_tree.GetVertice(start_nb_id), random_angles_rad, map, numofDOFs, x_size, y_size, Epsilon);
      if(start_extend_angles_rad == NULL) continue;
      int start_extend_id = start_tree.AddVertex(start_extend_angles_rad);
      start_tree.AddEdge(start_nb_id, start_extend_id);

      //Extend goal_tree to extend_angles
      int goal_nb_id = goal_tree.GetNearestNeibour(start_extend_angles_rad);
      double* goal_extend_angles_rad = Extend(goal_tree.GetVertice(goal_nb_id), start_extend_angles_rad, map, numofDOFs, x_size, y_size, Epsilon);
      if(goal_extend_angles_rad == NULL) continue;
      int goal_extend_id = goal_tree.AddVertex(goal_extend_angles_rad);
      goal_tree.AddEdge(goal_nb_id,goal_extend_id);
      if(angles_equal(start_extend_angles_rad, goal_extend_angles_rad, numofDOFs))
      {
        mexPrintf("Caught\n");
        found = true;
        break;
      }
    }

    else
    {
      //Extend goal tree to random_angles
      GenerateRandomConfig(random_angles_rad, map, numofDOFs, x_size, y_size);
      int goal_nb_id = goal_tree.GetNearestNeibour(random_angles_rad);
      double* goal_extend_angles_rad = Extend(goal_tree.GetVertice(goal_nb_id), random_angles_rad, map, numofDOFs, x_size, y_size, Epsilon);
      if(goal_extend_angles_rad == NULL) continue;
      int goal_extend_id = goal_tree.AddVertex(goal_extend_angles_rad);
      goal_tree.AddEdge(goal_nb_id, goal_extend_id);

      //Extend start tree to extend_angles
      int start_nb_id = start_tree.GetNearestNeibour(goal_extend_angles_rad);
      double* start_extend_angles_rad = Extend(start_tree.GetVertice(start_nb_id), goal_extend_angles_rad, map, numofDOFs, x_size, y_size, Epsilon);
      if(start_extend_angles_rad == NULL) continue;
      int start_extend_id = start_tree.AddVertex(start_extend_angles_rad);
      start_tree.AddEdge(start_nb_id,start_extend_id);
      if(angles_equal(goal_extend_angles_rad, start_extend_angles_rad, numofDOFs))
      {
        mexPrintf("Caught\n");
        found = true;
        break;
      }      
    }

    extend_stree = !extend_stree;
  }

  if(found)
  {
    vector<int> path;
    int next_id = start_tree.GetLastId();
    while(next_id!=0)
    {
      path.insert(path.begin(),next_id);
      next_id = start_tree.GetParent(next_id);
    }
    path.insert(path.begin(),start_tree.GetRootId());

    int start_tree_size = path.size();

    next_id = goal_tree.GetLastId();

    while(next_id!=0)
    {
      next_id = goal_tree.GetParent(next_id);
      path.push_back(next_id);
    }

    *planlength = path.size();

    *plan = (double**) malloc(path.size()*sizeof(double*));

    for(int i=0; i<path.size(); i++)
    {
      (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
      if(i<start_tree_size)
        memcpy((*plan)[i], start_tree.GetVertice(path[i]), numofDOFs*sizeof(double));
      else
        memcpy((*plan)[i], goal_tree.GetVertice(path[i]), numofDOFs*sizeof(double));
    }    

  }

  else
  {
    mexPrintf("Path Not Found!\n");
  }

  mexPrintf("Number of Samples: %d\n", num_samples);
}


//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    
    /* Check for proper number of arguments */    
    if (nrhs != 4) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Four input arguments required."); 
    } else if (nlhs != 2) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))){
        	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);
 
    //get the planner id
    int planner_id = (int)*mxGetPr(PLANNER_ID_IN);
    if(planner_id < 0 || planner_id > 3){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidplanner_id",
                "planner id should be between 0 and 3 inclusive");         
    }
    

    //call the planner
    double** plan = NULL;
    int planlength = 0;

    bool valid_start_goal = true;

    //Check the start and goal configuration.
    if(!IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size) ||
       !IsValidArmConfiguration(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size))
    {
      mexPrintf("Invalid start or goal configuration!\n");
      valid_start_goal = false;
    }
      
    
    //you can may be call the corresponding planner function here
    //if (planner_id == RRT)
    //{
    //    plannerRRT(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    //}
    
    //dummy planner which only computes interpolated path

    if(valid_start_goal)
    {
      clock_t t1=clock();

      switch(planner_id)
      {
        case RRT:
        {
          mexPrintf("The planner is RRT\n");
          plannerRRT(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
          break;
        }
        case RRTCONNECT:
        {
          mexPrintf("The planner is RRTCONNECT\n");
          plannerRRTCONNECT(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
          break;
        }
        case RRTSTAR:
        {
          mexPrintf("The planner is RRTSTAR\n");
          plannerRRTSTAR(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
          break;
        }
        case PRM:
        {
          mexPrintf("The planner is PRM\n");
          plannerPRM(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
          break;
        }
        default: mexPrintf("Not valid id\n");
        break;
      }
     
      mexPrintf("The planning takes %f\n", (float)(clock()-t1)/CLOCKS_PER_SEC);
      
      printf("planner returned plan of length=%d\n", planlength); 
    }
    
    /* Create return values */
    if(planlength > 0)
    {
        double plan_dist = 0;

        for(int i=0; i<planlength-1; i++)
        {
          plan_dist = plan_dist + dist_angles(plan[i], plan[i+1],numofDOFs);
        }

        mexPrintf("Total dist of the plan = %f\n", plan_dist);

        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
                plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    int* planlength_out = (int*) mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    
    return;
    
}





