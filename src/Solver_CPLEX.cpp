

/*****************************************************************************************\
**
** Solver Class
**
** This file contains the routines needed to access the CPLEX solver
**
**
**
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
\******************************************************************************************/

#include"Solver_CPLEX.h"


void Solver_CPLEX::open_solver()
{
	env = IloEnv();
}

void Solver_CPLEX::end_solver()
{
	env.end();
}

void Solver_CPLEX::open_prob(Prob& prob)
{
	prob.env = env;
	prob.model = IloModel(prob.env);
	prob.cplex = IloCplex(prob.env);
	prob.obj_raw = IloObjective(prob.env);
	prob.vars_raw = IloNumVarArray(prob.env);
	prob.surro_vars_raw = IloNumVarArray(prob.env);
	prob.R_vars_raw = IloNumVarArray(prob.env);
	prob.range_raw = IloRangeArray(prob.env);
}

void Solver_CPLEX::Clear_Prob(Prob& prob, string type)
{
	if (type == "optimality")
	{
		prob.model.remove(prob.benders_opt);
		prob.benders_opt.clear();
	}
	else if (type == "feasibility")
	{
		prob.model.remove(prob.benders_feas);
		prob.benders_feas.clear();
	}
}

void Solver_CPLEX::Create_Var(Var_struct& var)
{
	if (var.type == 0)
	{
		var.x = IloNumVar(env, var.lb, var.ub, ILOBOOL, var.name);
	}
	else if (var.type == 1)
	{
		var.x = IloNumVar(env, var.lb, var.ub, ILOFLOAT, var.name);
	}
	else {
		var.x = IloNumVar(env, var.lb, var.ub, ILOINT, var.name);
	}
	
}



void Solver_CPLEX::Create_RNG(RNG_struct& rng)
{
	if (rng.type == 0)
	{
		rng.rng = rng.expr <= rng.ub;
		rng.rng.setName(rng.name);
	}
	else if (rng.type == 1)
	{
		rng.rng = rng.lb <= rng.expr;
		rng.rng.setName(rng.name);
	}
	else 
	{
		rng.rng = IloRange(env, rng.lb, rng.expr, rng.ub, rng.name);
	}
}



void Solver_CPLEX::Create_OBJ(OBJ_struct& obj)
{
	if (obj.type == 0)
	{
		obj.obj = IloObjective(env, obj.expr, IloObjective::Maximize, obj.name);
	}
	else 
	{
		obj.obj = IloObjective(env, obj.expr, IloObjective::Minimize, obj.name);
	}
}



void Solver_CPLEX::get_Linear_obj_coeffs(Prob& prob)
{
	//Get OBJ coefficients
	for (IloExpr::LinearIterator it = IloExpr(prob.obj_raw.getExpr()).getLinearIterator(); it.ok(); ++it) {
		Coeff_Sparse empty_coeff;
		empty_coeff.row = 0;
		empty_coeff.row_name = prob.obj_raw.getName();
		empty_coeff.col = prob.vars_raw.find(it.getVar());
		empty_coeff.col_name = it.getVar().getName();
		empty_coeff.val = it.getCoef();
		prob.obj_coef_raw.push_back(empty_coeff);
	}

	//std::sort(prob.obj_coef_raw.begin(), prob.obj_coef_raw.end());

}

void Solver_CPLEX::get_Linear_rng_coeffs(Prob& prob)
{
	 //Get RNG coefficients
	prob.rng_coefs_raw.resize(prob.range_raw.getSize());
	for (int r = 0; r < prob.range_raw.getSize(); r++)
	{
		for (IloExpr::LinearIterator it = IloExpr(prob.range_raw[r].getExpr()).getLinearIterator(); it.ok(); ++it) {
			Coeff_Sparse empty_coeff;
			empty_coeff.row = r;
			empty_coeff.row_name = prob.range_raw[r].getName();
			empty_coeff.col = prob.vars_raw.find(it.getVar());
			empty_coeff.col_name = it.getVar().getName();
			empty_coeff.val = it.getCoef();
			prob.rng_coefs_raw[r].push_back(empty_coeff);
		}
		
		//std::sort(prob.rng_coefs_raw[r].begin(), prob.rng_coefs_raw[r].end());
	}
}

void Solver_CPLEX::Create_Prob(Prob& prob)
{
	//prob.num_var = prob.vars_raw.getSize();
	//prob.num_rng = prob.range_raw.getSize();
	prob.model.add(prob.vars_raw);
	//prob.model.add(prob.surro_vars_raw);
	prob.model.add(prob.obj_raw);
	prob.model.add(prob.range_raw);
	prob.cplex.extract(prob.model);
	//prob.cplex.setName(prob.name);
	//prob.cplex.setOut(env.getNullStream());
	//prob.cplex.setWarning(env.getNullStream());
	//prob.cplex.setParam(IloCplex::TiLim, TimeLim);
	//prob.cplex.setParam(IloCplex::EpRHS, 1e-4);
	//prob.cplex.setParam(IloCplex::EpOpt, 1e-4);
	//prob.cplex.setParam(IloCplex::EpGap, EPGAP);
	//prob.cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);

}


string Solver_CPLEX::getProbtype(Prob& prob) 
{

	if (prob.cplex.isMIP()) 
	{
		if (prob.cplex.isQC())       return "MIPQP";
		else if (prob.cplex.isQO())  return "MIQP";
		else                         return "MILP";
	}
	else 
	{
		if (prob.cplex.isQC())         return "QCP";
		else if (prob.cplex.isQO())    return "QP";
		else                           return "LP";
	}

	return "Prob Type ERROR";

}

void Solver_CPLEX::setDefault(Prob& prob) 
{

	/*  SET PARAMETERS TO DEFAULT (CONFIG FILE IS USED)  */
	prob.cplex.setOut(env.getNullStream());
	prob.cplex.setWarning(env.getNullStream());
	//prob.cplex.setParam(IloCplex::RootAlg, IloCplex::Algorithm::Dual);
	//prob.cplex.setParam(IloCplex::EpGap, EPGAP);
	//prob.cplex.setParam(IloCplex::TiLim, TimeLim);
	//cplex.setParam(IloCplex::Threads, 1);
}

void Solver_CPLEX::open_rng(IloRangeArray& rng, int iter)
{
	if (iter == 0)
	{
		rng = IloRangeArray(env);
	}
}

IloExpr Solver_CPLEX::set_QP_obj(Prob& probl, vector<SOL_str> xhat)
{
	
	IloExpr empty_expr(probl.env);

	for (int v = 0; v < probl.vars_raw.getSize(); v++)
	{
		empty_expr += probl.lambda * (probl.vars_raw[v] - xhat[v].value);
		empty_expr += (probl.sigma / 2) * probl.vars_raw[v] * probl.vars_raw[v];
		empty_expr -= probl.sigma * probl.vars_raw[v] * xhat[v].value;
		empty_expr += (probl.sigma / 2) * xhat[v].value * xhat[v].value;
	}


	return empty_expr;
}

#pragma region Usual Functions: Name, ...

void Solver_CPLEX::ImportModel(Prob& prob)
{

	try {
		prob.cplex.importModel(prob.model, prob.name, prob.obj_raw, prob.vars_raw, prob.range_raw);
		prob.isRaw = true; //the model is not distributed in detailed structs, rather in raw structs
		prob.probType = getProbtype(prob);
		prob.num_var = prob.vars_raw.getSize();
		prob.num_rng = prob.range_raw.getSize();
	}
	catch (IloException &except) {
		cout << except << endl;
		exit(1);
	}

}

const char *Solver_CPLEX::getName(IloExtractable &any_object)
{
	return any_object.getName();
}

int Solver_CPLEX::getSize(IloArray<IloExtractable> &any_object)
{
	return any_object.getSize();
}

void Solver_CPLEX::add_to_array(IloArray<IloExtractable> &any_object_array, IloExtractable &any_object)
{
	any_object_array.add(any_object);
}

bool Solver_CPLEX::Solve_Prob(Prob& prob, bool isEval)
{
	//cout << prob.model << endl;
	//int empty;
	//cin >> empty;

	if (prob.cplex.solve())
	{
		prob.zstar = prob.cplex.getObjValue();
		prob.sol.clear();
		double ext = 0;

		if (!isEval)
		{
			prob.duals.clear();
			prob.sol.clear();
			prob.sol_surrogate.clear();

			for (int v = 0; v < prob.vars_raw.getSize(); v++)
			{
				SOL_str empt;
				empt.col = prob.vars_raw[v].getName();
				empt.value = prob.cplex.getValue(prob.vars_raw[v]);
				empt.col_num = v;
				empt.var = &prob.vars_raw[v];
				prob.sol.push_back(empt);
				if (prob.isReg)
				{
					ext -= (prob.sigma / 2) * empt.value * empt.value;
					ext += prob.sigma * empt.value * prob.xhat[v].value;
					prob.zstar -= (prob.sigma / 2) * prob.xhat[v].value * prob.xhat[v].value;
				}
			}
			
			if (prob.isReg)
			{
				//cout << " " << ext << " , " << prob.zstar << " ";
				prob.zstar += ext;
				prob.zstar_wout_surro = prob.zstar;
			}

			prob.zstar_wout_surro = prob.zstar;

			if (prob.has_surrogate)
			{
				prob.sol_surrogate.clear();
				for (int v = 0; v < prob.surro_vars_raw.getSize(); v++)
				{
					SOL_str empt;
					empt.col = prob.surro_vars_raw[v].getName();
					empt.value = prob.cplex.getValue(prob.surro_vars_raw[v]);
					prob.zstar_wout_surro -= empt.value;
					prob.sol_surrogate.push_back(empt);
				}

			}

			prob.duals.clear();
			for (int r = 0; r < prob.range_raw.getSize(); r++)
			{
				IloNum dual = 0;
				dual = GetDual(prob.range_raw[r], prob.cplex);
				SOL_str empt;
				empt.row = prob.range_raw[r].getName();
				empt.value = dual;
				prob.duals.push_back(empt);
			}

			if (prob.opt)
			{
				prob.opt_cut_duals.clear();
				for (int r = 0; r < prob.benders_opt.getSize(); r++)
				{
					IloNum dual = 0;
					dual = GetDual(prob.benders_opt[r], prob.cplex);
					SOL_str empt;
					empt.row = prob.benders_opt[r].getName();
					empt.value = dual;
					prob.opt_cut_duals.push_back(empt);
				}
			}
		}

		return true;

	}
	else
	{
		cout << "ERROR: Solver: Problem is not Solved!" << endl;
		cout << "Status: " << prob.cplex.getStatus() << endl;
		cout << prob.model << endl;
		system("pause");
		return false;
	}
}

void Solver_CPLEX::Print_Prob(Prob& prob)
{
	prob.cplex.exportModel(prob.model_name.c_str());
}

void Solver_CPLEX::set_rhs_val(IloRange& rng, IloNum val)
{
	if (rng.getUB() <
		INFINITY && rng.getLB() > -INFINITY) {
		rng.setBounds(val, val);
	}
	else if (rng.getUB() < INFINITY) {
		rng.setUB(val);
	}
	else {
		rng.setLB(val);
	}

}

IloNum Solver_CPLEX::get_rhs(IloRange& rng)
{
	//cout << rng.getUB() << "  " << rng.getLB() << endl;
	//cout << INFINITY << "  " << -INFINITY << endl;
	if (rng.getUB() <
		INFINITY && rng.getLB() > -INFINITY) {
		return rng.getLB();
	}
	else if (rng.getUB() < INFINITY) {
		return rng.getUB();
	}
	else {
		return rng.getLB();
	}
}

//decomposition of the range based on its variables and the range reference from mean_prob
void Solver_CPLEX::decompose_range(IloRangeArray& mean_rng, Prob& prob, IloNumVarArray vars, int start, int end)
{
	
	IloExpr empty_expr(prob.env);
	int expr_size = 0;
	int range_idx = 0;

	for (int r = start; r < end; r++)
	{

		vector<Coeff_Sparse> empty_coef;
		vector<Coeff_Sparse> empty_prev_coef;
		IloNumVarArray empt_prev_vars(prob.env);

		for (IloExpr::LinearIterator it = mean_rng[r].getLinearIterator(); it.ok(); ++it)
		{
			if (vars.find(it.getVar()) <= vars.getSize() && vars.find(it.getVar()) >= 0)
			{
				empty_expr += it.getVar() * it.getCoef();
				Coeff_Sparse empt;
				empt.col = (int)vars.find(it.getVar());
				empt.col_name = it.getVar().getName();
				empt.row = r;
				empt.row_name = mean_rng[r].getName();
				empt.val = it.getCoef();
				empty_coef.push_back(empt);
				expr_size++;
			}
			else
			{
				//cout << "prev var: " << it.getVar() << endl;
				Coeff_Sparse empt;
				empt.col_name = it.getVar().getName();
				empt.row = r;
				empt.row_name = mean_rng[r].getName();
				empt.val = it.getCoef();
				empty_prev_coef.push_back(empt);
				empt_prev_vars.add(it.getVar());
			}
		}//loop over variables

		if (expr_size != 0)
		{
			IloRange empty_range = Create_RNG_explicit(mean_rng[r].getLB(), mean_rng[r].getUB(),
				rng_type(mean_rng[r]), empty_expr, mean_rng[r].getName());

			prob.range_raw.add(empty_range);

			prob.prev_rng_coefs_raw.push_back(empty_prev_coef);
			prob.rng_coefs_raw.push_back(empty_coef);
			//cout << empt_prev_vars.getSize() << " ";
			prob.prev_vars_raw.push_back(empt_prev_vars);
			//cout << prob.prev_vars_raw[prob.prev_vars_raw.size() - 1].getSize() << endl;
			//cout << prob.prev_vars_raw.size() << endl;
		}

		expr_size = 0;

		empty_coef.clear();
		empty_prev_coef.clear();
		empty_expr.clear();
		//cout << prob.prev_vars_raw[0].getSize() << endl;
	}//loop over ranges

	//for (int i = 0; i < prob.prev_vars_raw.size(); i++)
	//{
	//	for (int j = 0; j < prob.prev_vars_raw[i].getSize(); j++)
	//	{
	//		cout << i << " " << j << endl;
	//	}
	//}
	
}

void Solver_CPLEX::set_rhs_var(Prob& prob, vector<SOL_str>& var, vector<SOL_str>& rhs)
{
	for (int r = 0; r < prob.range_raw.getSize(); r++)
	{
		IloNum Cx = 0;
		//Calculate C.x
		
		for (int i = 0; i < var.size(); i++)
		{
			int pos = prob.prev_vars_raw[r].find(*var[i].var);
			//cout << r << " " << prob.prev_vars_raw[r].getSize() << " ";
			//cout << *var[i].var << " " << pos << " " << endl;
			if (pos != -1)
			{
				//cout << prob.prev_rng_coefs_raw[r][pos].val << " ";
				//cout << var[i].value << " ";
				Cx += prob.prev_rng_coefs_raw[r][pos].val * var[i].value;
			}	
		}
		//cout << endl;
		//cout << Cx << endl;
		//set rhs
		set_rhs_val(prob.range_raw[r], rhs[r].value - Cx);
	}
}

void Solver_CPLEX::set_rhs(Prob& prob, vector<SOL_str> rhs)
{
	for (int r = 0; r < prob.range_raw.getSize(); r++)
	{
		//set rhs
		set_rhs_val(prob.range_raw[r], rhs[r].value);
	}
}

//// Creating CPLEX objects
IloNumVar Solver_CPLEX::Create_Var_explicit(IloNum lb, IloNum ub, int type, const char *name)
{
	if (type == 0) //creating a binary variable
	{
		return IloNumVar(env, lb, ub, ILOBOOL, name);
	}
	else if (type == 1) //creating a continuous variable
	{
		return IloNumVar(env, lb, ub, ILOFLOAT, name);
	}
	else //create an integer variable
	{
		return IloNumVar(env, lb, ub, ILOINT, name);
	}
}

int Solver_CPLEX::rng_type(IloRange& rng)
{
	if (rng.getUB() <
		INFINITY && rng.getLB() > -INFINITY) {
		return 2;
	}
	else if (rng.getUB() < INFINITY) {
		return 0;
	}
	else {
		return 1;
	}
}

IloRange Solver_CPLEX::Create_RNG_explicit(IloNum lb, IloNum ub, int type, IloExpr expr, const char *name)
{
	IloRange rng;

	if (type == 0)       // <= type of constraints
	{
		rng = expr <= ub;
		rng.setName(name);
		return rng;
	}
	else if (type == 1) // >= type of constraints
	{
		rng = lb <= expr;
		rng.setName(name);
		return rng;
	}
	else               //equality constraints
	{
		rng = IloRange(env, lb, expr, ub, name);
		return rng;
	}
}

IloObjective Solver_CPLEX::Create_OBJ_explicit(IloObjective::Sense sense, IloExpr expr, const char *name)
{
	return IloObjective(env, expr, sense, name);
}
//*****************************
#pragma endregion collection of functions directly related to solvers like getting name and ...

#pragma region Optimization related Functions Duals, Rays, ...
IloNum Solver_CPLEX::GetDual(IloRange& rng, IloCplex& cplex)
{
	IloNum dual = 0;
	dual = cplex.getDual(rng);
	return dual;
}

IloNum Solver_CPLEX::GetDual_Var_Bound(IloNum val, IloNumVar var)
{
	if (val == var.getUB())
	{

	}
	else if (val == var.getLB())
	{

	}
	else
	{
		return 0.0;
	}
}
#pragma endregion functions related to optimization theory and complementary slackness


#pragma region Add Cutting Planes
void Solver_CPLEX::AddCutToModel(Prob& prob)
{
	prob.model.add(prob.benders_opt);
}
#pragma endregion //Adding Cutting Planes
