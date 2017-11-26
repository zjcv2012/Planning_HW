#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <functional>

using namespace std;

/* Input Arguments */
#define	BLOCKSV_IN      prhs[0]
#define	TRIANGLESV_IN      prhs[1]
#define	TABLEINDEX_IN      prhs[2]
#define	ONVSTART_IN      prhs[3]
#define	CLEARVSTART_IN      prhs[4]
#define	ONVGOAL_IN      prhs[5]
#define	CLEARVGOAL_IN      prhs[6]
#define	MOVEACTIONINDEX_IN      prhs[7]
#define	MOVETOTABLEACTIONINDEX_IN      prhs[8]


/* Output Arguments */
#define	PLAN_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

struct symbState
{
	unordered_set<int> clearV;
	unordered_map<int, int> onV;

	symbState(){}

	symbState(const symbState& state)
	{
		clearV = state.clearV;
		onV = state.onV;
	}

	symbState(int** onV_start, int onV_start_length, int* clearV_start, int numofclear_start)
	{
		for(int i=0; i<onV_start_length; i++)
		{
			onV[onV_start[i][0]]=onV_start[i][1];
		}

		for(int i=0; i<numofclear_start; i++)
		{
			clearV.insert(clearV_start[i]);
		}
	}
};

struct Action
{
	int action_index;
	int action_first;
	int action_second;
	int action_third;

	Action(){}

	Action(int index, int first, int second, int third): action_index(index), action_first(first), action_second(second), action_third(third){}
};


class symbPlanner
{
public:
	symbPlanner(int* blocksV, int numofblocks, int* trianglesV, int numoftriangles, int TableIndex, 
            int** onV_start, int onV_start_length, int* clearV_start, int numofclear_start, 
            int** onV_goal, int onV_goal_length, int* clearV_goal, int numofclear_goal, 
            int moveActionIndex, int moveToTableActionIndex)
	{

		symbState StartState(onV_start, onV_start_length, clearV_start, numofclear_start);
		StateSet.push_back(StartState);

		this->TableIndex = TableIndex;
		this->onV_goal = onV_goal;
		this->onV_goal_length = onV_goal_length;
		this->clearV_goal = clearV_goal;
		this->numofclear_goal = numofclear_goal;
		this->moveActionIndex = moveActionIndex;
		this->moveToTableActionIndex = moveToTableActionIndex;

		for(int i=0; i<numofblocks; i++)
		{
			Blocks.insert(blocksV[i]);
		}

		for(int i=0; i<numoftriangles; i++)
		{
			Triangles.insert(trianglesV[i]);
		}

	}

	int AddState(symbState state)
	{
		StateSet.push_back(state);
		return StateSet.size()-1;
	}

	bool stateEqual(symbState state1, symbState state2)
	{
		if(state1.onV.size()!=state2.onV.size() || state1.clearV.size()!=state2.clearV.size())
			return false;
		for(auto itr = state1.clearV.begin(); itr != state1.clearV.end(); itr++)
		{
			if(state2.clearV.count(*itr)==0)
				return false;
		}
		for(auto itr = state1.onV.begin(); itr != state1.onV.end(); itr++)
		{
			if(state2.onV.count(itr->first)==0)
				return false;
			else
			{
				if(state2.onV[itr->first]!=itr->second)
					return false;
			}
		}

		return true;
	}

	int checkExistence(symbState state)
	{
		for(int i=0; i<StateSet.size(); i++)
		{
			if (stateEqual(state, StateSet[i])==1)
				return i;
		}
		return -1;

	}

	int computeHeuristic(symbState state)
	{
		int h=0;
		for(int i=0; i<onV_goal_length; i++)
		{
			if(state.onV.count(onV_goal[i][0])==0)
				h++;
			else
			{
				if(state.onV[onV_goal[i][0]] != onV_goal[i][1])
					h++;
			}
		}

		for(int i=0; i<numofclear_goal; i++)
		{
			if(state.clearV.count(clearV_goal[i])==0)
				h++;
		}

		return h;
	}

	vector<int> getValidSuccessors(int id)
	{
		symbState currentState = StateSet[id];
		vector<int> ValidSuccessors;
		for(auto itr1=currentState.clearV.begin(); itr1!=currentState.clearV.end(); ++itr1)
		{
			for(auto itr2 = currentState.clearV.begin(); itr2!=currentState.clearV.end();++itr2)
			{
				if(*itr2 == *itr1 || Triangles.count(*itr2)>0) continue;
				symbState newState(currentState);

				if(newState.onV[*itr1]!=TableIndex)
				newState.clearV.insert(newState.onV[*itr1]);
				newState.onV[*itr1] = *itr2;
				newState.clearV.erase(*itr2);

				int c_id = checkExistence(newState); 


				if(c_id == -1)
				{
					int new_id = AddState(newState);
					action_map[new_id] = Action(moveActionIndex, *itr1, currentState.onV[*itr1],*itr2);
					ValidSuccessors.push_back(new_id);
					parent[new_id] = id;
				}
				else
				{
					tmp_action_map[c_id] = Action(moveActionIndex, *itr1, currentState.onV[*itr1],*itr2);
					ValidSuccessors.push_back(c_id);
				}
			}
			if(currentState.onV[*itr1]!=TableIndex)
			{
				symbState newState(currentState);
				newState.clearV.insert(newState.onV[*itr1]);
				newState.onV[*itr1] = TableIndex;

				int c_id = checkExistence(newState);

				if(c_id==-1)
				{
					int new_id = AddState(newState);
					action_map[new_id] = Action(moveToTableActionIndex, *itr1, currentState.onV[*itr1],TableIndex);
					ValidSuccessors.push_back(new_id);
					parent[new_id] = id;
				}
				else
				{
					tmp_action_map[c_id] = Action(moveToTableActionIndex, *itr1, currentState.onV[*itr1],TableIndex);
					ValidSuccessors.push_back(c_id);
				}
			}
		}

		return ValidSuccessors;

	}

	vector<Action> Plan()
	{
		unordered_set<int> openSet;
		unordered_set<int> closeSet;

		auto cmp = [this](int left, int right) 
		{ 
			int f_left = (this->g)[left]+1.5*(this->h)[left];
			int f_right = (this->g)[right]+1.5*(this->h)[right];
			//return (((this->g)[left])+(this->h[left]))<((this->g[right])+(this->h[right]));
			return f_left>f_right;
		};
		priority_queue<int, vector<int>, decltype(cmp)> pq(cmp);
		
		g[0] = 0;
		h[0] = computeHeuristic(StateSet[0]);
		//openSet.insert(0);

		pq.push(0);
		openSet.insert(0);

		int test_num = 0; 

		int goal_id;
		vector<Action> plan;

		bool reached = false;

		while(!pq.empty() && !reached && test_num<2000)
		{
			test_num++;
		/*	unordered_set<int>::iterator itr;
			int min_dist = 100000;
			int current_id;
			for(itr = openSet.begin(); itr != openSet.end(); itr++)
			{
				int tmp_dist = g[*itr] + h[*itr];
				//int tmp_dist = g[*itr];
				if(tmp_dist<min_dist)
				{
					current_id = *itr;
					min_dist = tmp_dist;
				}
				/*
				else if(tmp_dist == min_dist && computeHeuristic(StateSet[*itr])<computeHeuristic(StateSet[current_id]))
				{
					current_id = *itr;
					min_dist = tmp_dist;
				}*/
		//	}

		//	openSet.erase(current_id);*/

			int current_id = pq.top();
			pq.pop();

			openSet.erase(current_id);
			closeSet.insert(current_id);
			
			//mexPrintf("current_id: %d\n", current_id);			
			vector<int> ValidSuccessors = getValidSuccessors(current_id);

			//mexPrintf("test: %d\n", ValidSuccessors.size());
			for(auto successor: ValidSuccessors)
			{
				if(closeSet.count(successor)==1) continue;
				if(openSet.count(successor)==0)
				{
					g[successor] = g[current_id] + 1;
					h[successor] = computeHeuristic(StateSet[successor]);
					openSet.insert(successor);
					pq.push(successor);
					if(h[successor]==0)
					{
						goal_id = successor;
						reached = true;
						mexPrintf("reached!, goalId: %d\n", goal_id);
						mexPrintf("test_num: %d\n",test_num);
						break;
					}
				}
				else
				{
					if(g[successor]>g[current_id]+1)
					{
						g[successor] = g[current_id]+1;
						parent[successor] = current_id;
						action_map[successor] = tmp_action_map[successor];
					}
				}
			}
		}

		if(reached)
		{
			int next_id = goal_id;
			while(next_id!=0)
			{
				plan.insert(plan.begin(),action_map[next_id]);
				next_id = parent[next_id];
			}

		}
		return plan;

	}


private:
	symbState StartState;
	unordered_set<int> Blocks;
	unordered_set<int> Triangles;
	vector<symbState> StateSet;

	unordered_map<int, Action> action_map;
	unordered_map<int, Action> tmp_action_map;
	unordered_map<int, int> parent;


	int TableIndex;
	int moveActionIndex;
	int moveToTableActionIndex;

	int** onV_goal;
	int onV_goal_length;
	int* clearV_goal;
	int numofclear_goal;

	unordered_map<int, int> g;
	unordered_map<int, int> h;

};