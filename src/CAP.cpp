

/*****************************************************************************************\
**
** CAP Class
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

#include"CAP.h"
#define 	min(X, Y) 			((X) <= (Y) ? (X) : (Y))
#define 	max(X, Y) 			((X) >= (Y) ? (X) : (Y))


void CAP::Generate_CAP(string instance)
{
	read_data(instance);
	time_info.resize(2);
	create_core(instance);
	create_tim(instance);
	create_stoc(instance);
}

void CAP::read_data(string instance)
{
	string format = ".txt";
	string dirr = ".\\data\\CAP\\";
	string filename = dirr + instance + format;

	ifstream file(filename);

	size_t i = 0;
	for (; i < instance.length(); i++) { if (isdigit(instance[i])) break; }

	// remove the first chars, which aren't digits
	instance = instance.substr(i, instance.length() - i);

	// convert the remaining text to an integer
	int id = atoi(instance.c_str());

	file >> data.ORIG >> data.DEST;
	data.demand.resize(data.DEST);
	data.sigma.resize(data.DEST);
	data.fixed_cost.resize(data.ORIG);
	data.capacity.resize(data.ORIG);
	data.cost.resize(data.ORIG);
	for (int i = 0; i < data.ORIG; i++) data.cost[i].resize(data.DEST);

	string capacity3;
	if (filename == "capa.txt") {
		for (int i = 0; i < data.ORIG; i++) {
			file >> capacity3 >> data.fixed_cost[i];
			data.capacity[i] = 1200;
		}
	}
	else {
		for (int i = 0; i < data.ORIG; i++)
			file >> data.capacity[i] >> data.fixed_cost[i];
	}
	for (int i = 0; i < data.DEST; i++) {
		file >> data.demand[i];
		for (int j = 0; j < data.ORIG; j++)
			file >> data.cost[j][i];
	}
	double max_dem = 0;
	double tot_cap = 0;
	for (int i = 0; i < data.DEST; i++) {
		srand(i + id*1000);
		data.sigma[i] = 0.1*data.demand[i] + rand()%1 * 0.2 * data.demand[i];
		printf("\nsigma %d: %0.4f", i, data.sigma[i]);
		max_dem += data.demand[i] + 3 * data.sigma[i];
	}
	for (int i = 0; i < data.ORIG; i++) {
		srand(i+ data.DEST + id*1000);
		float center = max_dem / (data.ORIG / 2);
		float dev = rand()%10 * 0.05 * center;
		data.capacity[i] = 0.5*center + dev;
		printf("\ncapacity %d: %0.4f", i, data.capacity[i]);
		tot_cap += data.capacity[i];
	}
	printf("\ntotal capacity: %0.4f - maximum demand: %0.4f", tot_cap, max_dem);
	cout << "*************************** Data have beed read ***************************" << endl;

}

void CAP::create_core(string instance)
{
	solver.open_solver();

	Create_Vars(); //Create CAP variables vars: FAC_i(binary) + SHIP_i_j(Num)
	Create_RNGs(); //Create Constraints
	Create_OBJ(); //Create Objective Function

	vars[0].id = "FAC";
	vars[0].dim = 1;
	vars[0].var1 = vars_FAC;
	vars[1].id = "SHIP";
	vars[1].dim = 2;
	vars[1].var2 = vars_SHIP;

	Create_Prob(instance);
}

void CAP::create_tim(string instance)
{
	string type = ".txt";
	string dirr = ".\\models\\CAP\\";
	string name = dirr + instance + ".tim";

	ofstream file_tim(name);

	file_tim << "TIME" <<
		setw(10) << " " <<
		setw(10) << "CAP" <<
		setw(10) << " " << 
		setw(10) << " " <<'\n';
	file_tim << "PERIODS" << endl;
	file_tim << setw(10) << " " <<
		setw(10) << time_info[0][0] <<
		setw(10) << time_info[0][1] <<
		setw(10) << " " <<
		setw(10) << "TIME1" << '\n';
	file_tim << setw(10) << " " <<
		setw(10) << time_info[1][0] <<
		setw(10) << time_info[1][1] <<
		setw(10) << " " <<
		setw(10) << "TIME2" << '\n';
	file_tim << "ENDATA" << endl;

}

void CAP::create_stoc(string instance)
{
	string type = ".txt";
	string dirr = ".\\models\\CAP\\";
	string name = dirr + instance + ".sto";

	ofstream file_stoc(name);

	if (CAP_Random == 0) //If randomness is DISCRETE Distribution
	{
		file_stoc << "STOCH" <<
			setw(20) << "CAP" <<'\n';
		file_stoc << "INDEP" <<
			setw(20) << "DISCRETE" << '\n';
		for (int j = 0; j < data.DEST; j++) //create sections for demand
		{
			float val;
			float prob;
			float sum_prob = 0;
			for (int i = 3; i > 0; i--) 
			{
				val = max(0,data.demand[j] - data.sigma[j] * i);
				prob = 1 / (float)(2 * (i * 4));
				sum_prob += prob;
				file_stoc << setw(10) << "RHS" <<
					         setw(20) << "DEM_D" + to_string(j+1) <<
					         setw(10) << val <<
					         setw(20)  << prob <<'\n';
			}
			file_stoc << setw(10) << "RHS" <<
				setw(20) << "DEM_D" + to_string(j + 1) <<
				setw(10) << data.demand[j] <<
				setw(20) << 1 - sum_prob*2 << '\n';
			for (int i = 1; i < 4; i++)
			{
				val = data.demand[j] + data.sigma[j] * i;
				prob = 1 / (float)(2 * (i * 4));
				file_stoc << setw(10) << "RHS" <<
					setw(20) << "DEM_D" + to_string(j + 1) <<
					setw(10) << val <<
					setw(20) << prob  << '\n';
			}
			if(j < data.DEST - 1)  file_stoc << "*" << endl;
		}
		file_stoc << "ENDATA" << endl;
	}
}

void CAP::Create_Vars()
{
	vars.resize(2);
	//create binary variables for locating facilities
	vars_FAC.resize(data.ORIG);
	for (int i = 0; i < data.ORIG; i++)
	{
		vars_FAC[i].id = "FAC";
		char varName[100];
		sprintf(varName, "FAC_O%d", i+1);
		vars_FAC[i].name = varName;
		vars_FAC[i].lb = 0;
		vars_FAC[i].type = 0;
		vars_FAC[i].ub = 1;
		solver.Create_Var(vars_FAC[i]);
	}
	time_info[0].push_back("FAC_O1");

	//create continuous variables
	vars_SHIP.resize(data.ORIG);
	for (int i = 0; i < data.ORIG; i++) vars_SHIP[i].resize(data.DEST);
	for (int i = 0; i < data.ORIG; i++)
	{
		for (int j = 0; j < data.DEST; j++)
		{
			vars_SHIP[i][j].id = "SHIP";
			char varName[100];
			sprintf(varName, "SHIP_O%d_D%d", i+1, j+1);
			vars_SHIP[i][j].name = varName;
			vars_SHIP[i][j].lb = 0;
			vars_SHIP[i][j].type = 1;
			vars_SHIP[i][j].ub = IloInfinity;
			solver.Create_Var(vars_SHIP[i][j]);
		}
	}
	time_info[1].push_back("SHIP_O1_D1");

	cout << "*************************** Vars have beed created ***************************" << endl;

}

void CAP::Create_RNGs()
{
	
	//Create Surrogate \sum_i x_i cap_i >= max \sum_j d_j
	IloExpr exprs(solver.env);
	RNG_struct rngsu;
	for (int i = 0; i < data.ORIG; i++)
	{
		exprs += vars_FAC[i].x * data.capacity[i];
		vars_FAC[i].rng_coefs.push_back(2 * data.capacity[i]);
		string name = "SURROGATE";
		vars_FAC[i].rng_names.push_back(name);
	}
	IloNum rhs_surr = 0;
	for (int j = 0; j < data.DEST; j++)
	{
		rhs_surr += data.demand[j] + 3 * data.sigma[j];
	}
	rngsu.expr = exprs;
	rngsu.id = "SURROGATE";
	rngsu.israndom = false;
	rngsu.lb = rhs_surr;
	char varName[100];
	sprintf(varName, "SURROGATE");
	rngsu.name = varName;
	rngsu.rand_type = 2;
	rngsu.ub = INFINITY;
	rngsu.type = 1;
	solver.Create_RNG(rngsu);
	rngs.push_back(rngsu);
	//expr.clear();
	

	//Create Second Set of Constraints \sum_j y_i_j <= cap_i x_i \forall i
	for (int i = 0; i < data.ORIG; i++)
	{
		IloExpr expr(solver.env);
		RNG_struct rng;
		for (int j = 0; j < data.DEST; j++)
		{
			expr += vars_SHIP[i][j].x;
			vars_SHIP[i][j].rng_coefs.push_back(1.0);
			string name = "CAP_O" + to_string(i);
			vars_SHIP[i][j].rng_names.push_back(name);
		}
		expr -= 2 * data.capacity[i] * vars_FAC[i].x;
		vars_FAC[i].rng_coefs.push_back(-1 * 2 * data.capacity[i]);
		string name = "CAP_O" + to_string(i + 1);
		vars_FAC[i].rng_names.push_back(name);
		rng.expr = expr;
		rng.id = "Capacity";
		rng.israndom = false;
		rng.lb = -INFINITY;
		char varName[100];
		sprintf(varName, "CAP_O%d", i+1);
		rng.name = varName;
		rng.rand_type = 0;
		rng.ub = 0.0;
		rng.type = 0;
		solver.Create_RNG(rng);
		rngs.push_back(rng);
		//expr.clear();
	}
	time_info[1].push_back("CAP_O1");

	//Create First Set of Constraints \sum_i y_i_j >= demand_j \forall j
	for (int j = 0; j < data.DEST; j++)
	{
		IloExpr expr(solver.env);
		RNG_struct rng;
		for (int i = 0; i < data.ORIG; i++)
		{
			expr += vars_SHIP[i][j].x;
			vars_SHIP[i][j].rng_coefs.push_back(1.0);
			string name = "DEM_D" + to_string(j+1);
			vars_SHIP[i][j].rng_names.push_back(name);
		}
		rng.expr = expr;
		rng.id = "Demand";
		rng.israndom = true;
		rng.lb = data.demand[j];
		char varName[100];
		sprintf(varName, "DEM_D%d", j+1);
		rng.name = varName;
		rng.rand_type = 0;
		rng.ub = INFINITY;
		rng.type = 1;
		solver.Create_RNG(rng);
		rngs.push_back(rng);
		//expr.clear();
	}
	

	cout << "*************************** RNGs have beed created ***************************" << endl;

}

void CAP::Create_OBJ()
{
	IloExpr expr(solver.env);
	//Binary Fixed Cost
	for (int i = 0; i < data.ORIG; i++)
	{
		expr += data.fixed_cost[i] * vars_FAC[i].x;
		vars_FAC[i].obj_coef = data.fixed_cost[i];
	}

	//SHIP Costs
	for (int i = 0; i < data.ORIG; i++)
	{
		for (int j = 0; j < data.DEST; j++)
		{
			expr += data.cost[i][j] * vars_SHIP[i][j].x;
			vars_SHIP[i][j].obj_coef = data.cost[i][j];
		}
	}
	
	obj.expr = expr;
	obj.constant = 0.0;
	obj.id = "COST";
	obj.type = 1;
	char varName[100];
	sprintf(varName, "Cost");
	obj.name = varName;
	solver.Create_OBJ(obj);
	time_info[0].push_back("Cost");
	cout << "*************************** Obj have beed created ***************************" << endl;

}

void CAP::MakeRawProb()
{
	//Objective Function:
	CAP_prob.obj_raw = obj.obj;
	
	
	//RNGs
	for (int i = 0; i < rngs.size(); i++)
	{
		CAP_prob.range_raw.add(rngs[i].rng);
	}
	
	//Vars
	CAP_prob.rng_coefs_raw.resize(CAP_prob.range_raw.getSize());
	for (int i = 0; i < vars_FAC.size(); i++)
	{
		CAP_prob.vars_raw.add(vars_FAC[i].x);
	}
	for (int i = 0; i < vars_SHIP.size(); i++)
	{
		for (int j = 0; j < vars_SHIP[i].size(); j++)
		{
			CAP_prob.vars_raw.add(vars_SHIP[i][j].x);
		}
	}
	cout << "*************************** RAW Prob have beed created ***************************" << endl;
}

void CAP::Create_Prob(string instance)
{
	string type = ".mps";
	string dirr = ".\\models\\CAP\\";
	string name = dirr + instance + ".lp";
	string name2 = dirr + instance + ".mps";
	
	CAP_prob.env = solver.env;
	CAP_prob.name = "CAP";
	CAP_prob.model_name = name;
	CAP_prob.id = "CAP";
	CAP_prob.type = 2;
	solver.open_prob(CAP_prob);
	MakeRawProb();		  
	solver.Create_Prob(CAP_prob);
	CAP_prob.env = solver.env;
	CAP_prob.name = "CAP";
	CAP_prob.model_name = name2;
	CAP_prob.id = "CAP";
	CAP_prob.type = 2;
	solver.open_prob(CAP_prob);
	MakeRawProb();
	solver.Create_Prob(CAP_prob);				  



	cout << "*************************** Prob have beed created ***************************" << endl;


	solver.Print_Prob(CAP_prob);
	cout << "*************************** Prob have beed printed ***************************" << endl;

}