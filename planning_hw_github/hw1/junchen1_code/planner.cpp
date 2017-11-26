/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include "planner.h"
#include <time.h>
using namespace std;

typedef float PrimArray[NUMOFDIRS][NUMOFPRIMS][NUMOFINTERSTATES][NUMOFDIM];

int temp = 0;
float stored_goalposeX=0;
float stored_goalposeY=0;
vector<int> motion_vector;
int motion_num = 0;

bool first_time = true;

unordered_map<int,float> h; //the key of the map only considers poseX and poseY


//Added functions for planner

//transform global position to index
int pose2index(float poseX, float poseY, int dir, int x_size, int y_size)
{
    int gridposx = (int)(poseX / RES + 0.5);
    int gridposy = (int)(poseY / RES + 0.5);

    return ((gridposy-1)*x_size + (gridposx-1))*NUMOFDIRS + dir;
}

//transform index to global position
pair<float,float> index2pose(int index, int x_size, int y_size)
{
    index = index / NUMOFDIRS;
    int gridposy = index/x_size+1;
    int gridposx = index-x_size*(gridposy-1)+1;
    return make_pair(gridposx*RES, gridposy*RES);
}

int index2dir(int index)
{
    return index % NUMOFDIRS;
}

int pose2index_h(float poseX, float poseY, int x_size, int y_size)
{
    int gridposx = (int)(poseX / RES + 0.5);
    int gridposy = (int)(poseY / RES + 0.5);

    return (gridposy-1)*x_size + (gridposx-1);
}

pair<float,float> index2pose_h(int index, int x_size, int y_size)
{
    int gridposy = index/x_size+1;
    int gridposx = index-x_size*(gridposy-1)+1;
    return make_pair(gridposx*RES, gridposy*RES);   
}


//heuristic funtion(euclidean distance)
float H(float poseX, float poseY, float goalX, float goalY, int x_size, int y_size)
{
    int tmp_index = pose2index_h(poseX,poseY,x_size,y_size);

    if(h.count(tmp_index)>0)
        return h[tmp_index];
    else
        return sqrt((poseX-goalX)*(poseX-goalX)+(poseY-goalY)*(poseY-goalY));
}

float dist(float poseX, float poseY, float goalX, float goalY)
{
    return sqrt((poseX-goalX)*(poseX-goalX)+(poseY-goalY)*(poseY-goalY));
}

//recover motion from start_index and end_index
int recover_motion(int start_index, int end_index, int x_size, int y_size)
{
    int start_dir = index2dir(start_index);
    int end_dir = index2dir(end_index);

    if(start_dir!=end_dir)
    {
        if(start_dir<end_dir)
        {
            if(end_dir-start_dir==1) 
                return 3;
            else
                return 4;
        }
        else
        {
            if(start_dir-end_dir==1) 
                return 4;
            else 
                return 3;
        }
    }
    else
    {
        pair<float,float> start_pose = index2pose(start_index,x_size,y_size);
        pair<float,float> end_pose = index2pose(end_index,x_size,y_size);
        if (dist(start_pose.first,start_pose.second,end_pose.first,end_pose.second)>0.5)
            return 1;
        else
        {
            if(start_dir==0||start_dir==7||start_dir==1)
            {
                if(end_pose.first>start_pose.first) 
                    return 0;
                else 
                    return 2;
            }
            else if(start_dir==3||start_dir==4||start_dir==5)
            {
                if(end_pose.first<start_pose.first) 
                    return 0;
                else 
                    return 2;
            }
            else if(start_dir==2)
            {
                if(end_pose.second>start_pose.second) 
                    return 0;
                else 
                    return 2;
            }
            else
            {
                if(end_pose.second<start_pose.second) 
                    return 0;
                else 
                    return 2;
            }
        }
    }
}




bool applyaction(double *map, int x_size, int y_size, float robotposeX, float robotposeY, float robotposeTheta,
                 float *newx, float *newy, float *newtheta, PrimArray mprim, int dir, int prim)
{
    int i;
    for (i = 0; i < NUMOFINTERSTATES; i++) {
        *newx = robotposeX + mprim[dir][prim][i][0];
        *newy = robotposeY + mprim[dir][prim][i][1];
        *newtheta = mprim[dir][prim][i][2];
        
        int gridposx = (int)(*newx / RES + 0.5);
        int gridposy = (int)(*newy / RES + 0.5);

        /* check validity */
        if (gridposx < 1 || gridposx > x_size || gridposy < 1 || gridposy > y_size){
            return false;
        }

        if ((int)map[GETMAPINDEX(gridposx, gridposy, x_size, y_size)] != 0){
            return false;
        }

    }


    return true;
}

int getPrimitiveDirectionforRobotPose(float angle)
{
    /* returns the direction index with respect to the PrimArray */
    /* normalize bw 0 to 2pi */
    if (angle < 0.0) {
        angle += 2 * M_PI;
    }
    int dir = (int)(angle / (2 * M_PI / NUMOFDIRS)+0.5);
    if (dir == 8) {
        dir = 0;
    }
    return dir;
}

//Implement Dijkstra's algorithm to calculate heuristic 
void generateH(float startposeX, float startposeY, float goalposeX, float goalposeY, int x_size, int y_size, double* map)
{

    mexPrintf("Start to build the map\n");
    unordered_set<int> open_set;
    unordered_set<int> closed_set;    
    h.clear();

    vector<pair<float,float>> neighbours = {{-1*RES,0}, {RES,0}, {0,-1*RES}, {0,RES}};
    int start_index = pose2index_h(startposeX, startposeY, x_size, y_size);
    int goal_index = pose2index_h(goalposeX, goalposeY, x_size, y_size);

    h[goal_index] = 0;
    open_set.insert(goal_index);


    clock_t t1=clock();



    while(!open_set.empty() && (float)(clock()-t1)/CLOCKS_PER_SEC<0.9)
    {
        unordered_set<int>::iterator itr;
        unordered_set<int>::iterator close_index;
        

        double min_f =100000;

        for(itr=open_set.begin(); itr!=open_set.end(); ++itr)
        {
            double f = h[*itr];
            if (f<min_f)
            {
                min_f = f;
                close_index = itr;
            }
        }

        int current_index = *close_index;
        open_set.erase(close_index);
        closed_set.insert(current_index);


        for(auto i:neighbours)
        {
            float new_x =index2pose_h(current_index,x_size,y_size).first + i.first;
            float new_y =index2pose_h(current_index,x_size,y_size).second+ i.second;
            int gridposx = (int)(new_x / RES + 0.5);
            int gridposy = (int)(new_y / RES + 0.5);

            if (!(gridposx < 1 || gridposx > x_size || gridposy < 1 || gridposy > y_size))
            {
                if((int)map[GETMAPINDEX(gridposx, gridposy, x_size, y_size)] == 0)
                {
                    int new_index = GETMAPINDEX(gridposx, gridposy, x_size, y_size);
                    if(closed_set.count(new_index)==0)
                    {
                        if(open_set.count(new_index)==0)
                        {
                            open_set.insert(new_index);
                            h[new_index] = h[current_index] + (float)RES*(float)sqrt(0.5);
                        }
                        else
                        {
                            h[new_index] = min(h[new_index], h[current_index] + (float)RES*(float)sqrt(0.5));
                        }
                    }
                }
            }
        }

    }

    return;
}

static void planner(
		   double*	map,
		   int x_size,
 		   int y_size,
           float robotposeX,
            float robotposeY,
            float robotposeTheta,
            float goalposeX,
            float goalposeY,
            PrimArray mprim,
            int *prim_id)
{   
    mexPrintf("temp=%d\n", temp);
    temp = temp+1;

    *prim_id = 0; /* arbitrary action */
   if(!motion_vector.empty()) 
    {
        *prim_id=motion_vector[motion_num];
        motion_num = motion_num + 1;

        // If the robot is close enough to the goal or all motions are implemented
        // clear the motion_vector and replan        
        if(motion_num>=motion_vector.size() || (fabs(robotposeX-goalposeX)<5 && fabs(robotposeY-goalposeY)<5))
        {
            motion_vector.clear();
            motion_num = 0;
            return ;
        }
    }
    else
    {
        //For A* algorithm
        unordered_set<int> open_set;
        unordered_set<int> closed_set;
        unordered_map<int,int> parent;
        unordered_map<int,float> g;

        clock_t t1=clock();

        if(stored_goalposeX!= goalposeX) stored_goalposeX = goalposeX;
        if(stored_goalposeY!= goalposeY) stored_goalposeY = goalposeY;

        if(first_time)
        {
            generateH(robotposeX, robotposeY, stored_goalposeX, stored_goalposeY, x_size, y_size, map);
            mexPrintf("The mapping takes %f\n", (float)(clock()-t1)/CLOCKS_PER_SEC);
            t1 = clock();
        }



        int start_index = pose2index(robotposeX,robotposeY,getPrimitiveDirectionforRobotPose(robotposeTheta),x_size,y_size);
        int goal_index;
        open_set.insert(start_index);
        g[start_index]= 0 ;

        bool caught= false; 


        while(!caught && !open_set.empty())
        {

            unordered_set<int>::iterator itr;
            unordered_set<int>::iterator close_index;
            double min_f =100000;

            for(itr=open_set.begin(); itr!=open_set.end(); ++itr)
            {
                double f = g[*itr]+ 3*H(index2pose(*itr,x_size,y_size).first,index2pose(*itr,x_size,y_size).second, 
                    stored_goalposeX, stored_goalposeY, x_size, y_size);
                if (f<min_f)
                {
                    min_f = f;
                    close_index = itr;
                }
            }

            int current_index = *close_index;
        
            closed_set.insert(current_index);
            open_set.erase(close_index);

            float robotposeX = index2pose(current_index,x_size,y_size).first;
            float robotposeY = index2pose(current_index,x_size,y_size).second;
            int dir = index2dir(current_index);

            //mexPrintf("current pose is: %f, %f, %d\n",robotposeX,robotposeY, dir);

            
            for(int prim=0; prim<NUMOFPRIMS; prim++)
            {
                float new_x,new_y,new_theta;
                bool ret = applyaction(map, x_size, y_size, robotposeX, robotposeY, robotposeTheta, &new_x, &new_y, &new_theta, mprim, dir, prim);
                if(ret)
                {
                    //mexPrintf("new pose is: %f, %f, %d\n",new_x, new_y, getPrimitiveDirectionforRobotPose(new_theta));
                    
                    if(fabs(new_x-stored_goalposeX)<thresh && fabs(new_y-stored_goalposeY)<thresh)
                    {
                        mexPrintf("Caught!\n");
                        caught=true;
                        goal_index = pose2index(new_x, new_y, getPrimitiveDirectionforRobotPose(new_theta), x_size, y_size);
                        parent[goal_index]= current_index;
                        mexPrintf("The planning takes %f\n", (float)(clock()-t1)/CLOCKS_PER_SEC);
                        break;
                    }   
                    int new_index = pose2index(new_x, new_y, getPrimitiveDirectionforRobotPose(new_theta), x_size, y_size);
                    if(closed_set.count(new_index)==0)
                    {
                        if(open_set.count(new_index)==0)
                        {
                            open_set.insert(new_index);
                            g[new_index] = g[current_index]+dist(robotposeX,robotposeY,new_x,new_y);
                            parent[new_index] = current_index;
                        }
                        else
                        {
                            if(g[*close_index]+dist(robotposeX,robotposeY,new_x,new_y)<g[new_index])
                            {
                                g[new_index] = g[current_index]+dist(robotposeX,robotposeY,new_x,new_y);
                                parent[new_index] = current_index;
                            }
                        }
                    }
                }
            }
        }

        if(caught)
        {

            h.clear();

            int next_index=goal_index;
            while(next_index!=start_index)
            {
              
                //Save the motion 
                if(first_time)
                {
                    motion_vector.insert(motion_vector.begin(), recover_motion(parent[next_index],next_index, x_size,y_size));                    
                }
                else 
                {
                    if(parent[next_index]==start_index)
                        motion_vector.insert(motion_vector.begin(), recover_motion(parent[next_index],next_index, x_size,y_size));                        
                }                
                next_index = parent[next_index];

            }
            *prim_id=motion_vector[0];
            motion_num=1;


            first_time=false;

            //If there is only one motion, then replan 
            if(motion_vector.size() == 1)
                motion_vector.clear();
        }
    }


    mexPrintf("action %d\n", *prim_id);
    return;
}

/*prhs contains input parameters (3): 
1st is matrix with all the obstacles
2nd is a row vector <x,y> for the robot pose
3rd is a row vector <x,y> for the target pose
plhs should contain output parameters (1): 
1st is a row vector <dx,dy> which corresponds to the action that the robot should make*/

void parseMotionPrimitives(PrimArray mprim)
{
    FILE * fp;
    fp = fopen ("unicycle_8angles.mprim", "r+");
    char skip_c[100];
    int skip_f;
    float resolution;
    int num_angles;
    int num_mprims;
    fscanf(fp, "%s %f", skip_c, &resolution);
    fscanf(fp, "%s %d", skip_c, &num_angles);
    fscanf(fp, "%s %d", skip_c, &num_mprims);

    int i, j, k;
    for (i = 0; i < NUMOFDIRS; ++i) {
        for (j = 0; j < NUMOFPRIMS; ++j) {
            fscanf(fp, "%s %d", skip_c, &skip_f);
            for (k = 0; k < NUMOFINTERSTATES; ++k) {
                fscanf(fp, "%f %f %f", &mprim[i][j][k][0], &mprim[i][j][k][1], &mprim[i][j][k][2]);
            }

        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{

    /* Read motion primtives */
    PrimArray motion_primitives;
    parseMotionPrimitives(motion_primitives);

    /* Check for proper number of arguments */    
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Three input arguments required."); 
    } else if (nlhs != 1) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = mxGetM(MAP_IN);
    int y_size = mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the dimensions of the robotpose and the robotpose itself*/     
    int robotpose_M = mxGetM(ROBOT_IN);
    int robotpose_N = mxGetN(ROBOT_IN);
    if(robotpose_M != 1 || robotpose_N != 3){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidrobotpose",
                "robotpose vector should be 1 by 3.");         
    }
    double* robotposeV = mxGetPr(ROBOT_IN);
    float robotposeX = (float)robotposeV[0];
    float robotposeY = (float)robotposeV[1];
    float robotposeTheta = (float)robotposeV[2];
    
    /* get the dimensions of the goalpose and the goalpose itself*/     
    int goalpose_M = mxGetM(GOAL_IN);
    int goalpose_N = mxGetN(GOAL_IN);
    if(goalpose_M != 1 || goalpose_N != 3){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidgoalpose",
                "goalpose vector should be 1 by 3.");         
    }
    double* goalposeV = mxGetPr(GOAL_IN);
    float goalposeX = (float)goalposeV[0];
    float goalposeY = (float)goalposeV[1];
        
    /* Create a matrix for the return action */ 
    ACTION_OUT = mxCreateNumericMatrix( 1, 1, mxINT8_CLASS, mxREAL); 
    int* action_ptr = (int*) mxGetData(ACTION_OUT);

    /* Do the actual planning in a subroutine */
    planner(map, x_size, y_size, robotposeX, robotposeY, robotposeTheta, goalposeX, goalposeY, motion_primitives, &action_ptr[0]);

    return ;
    
}





