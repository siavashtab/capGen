
/*****************************************************************************************\
**
** CAP_header class
**
** This file contains the routines needed to generate cap instances
**
**
**
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
\******************************************************************************************/

#ifndef CAP_h
#define CAP_h

#include"Solver_CPLEX.h"
#include"Config.h"

struct CAP_data {
	int ORIG;
	int DEST;
	vec_f fixed_cost;
	vec_f capacity;
	vec_f demand;
	vec2_f cost;
	vec_f sigma;
};


class CAP {

public:

	void Generate_CAP(string instance);  //generate CAP model 
	                                     //from cap41-capa,capb,capc instances

	Prob CAP_prob;
	Prob master_prob;
	Prob ph_sub_prob;
	vector<Prob> stage_sub_prob;

private:

	void initialize();
	void read_data(string instance);  //read the parameters
	void create_tim(string instance);
	void create_stoc(string instance);
	void create_core(string instance);

	void Create_Vars();
	void Create_RNGs();
	void Create_OBJ();
	void Create_Prob(string instance);
	void MakeRawProb();
	void MakeMasterProb();
	void MakeSubProb();

	CAP_data data;

	Solver_CPLEX solver;
	vec_var vars_FAC;
	vec_var2 vars_SHIP;
	vec_rng rngs;
	vector<var_vec> vars;
	OBJ_struct obj;
	Prob prob;
	vec_d sigma;

	vector<vector<string>> time_info;

};

#endif // !CAP_h
