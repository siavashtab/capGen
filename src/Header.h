#ifndef HEADER_H
#define HEADER_H
#pragma warning(disable:4786)
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <ilcplex/ilocplex.h>
#define INFINITY 99999999999
using namespace std;

typedef vector<IloBoolVar> vec_bool;
typedef vector<vec_bool> vec_2bool;

typedef vector<IloNumVar> vec_num;
typedef vector<vec_num> vec_2num;

typedef vector<IloIntVar> vec_int;
typedef vector<vec_int> vec_2int;

typedef vector<int> vec_i;
typedef vector<double> vec_d;
typedef vector<float> vec_f;
typedef vector<vec_i> vec2_i;
typedef vector<vec_d> vec2_d;
typedef vector<vec2_d> vec3_d;
typedef vector<vec_f> vec2_f;


typedef IloArray<IloNumVarArray>  IloNumVarArray2;
typedef IloArray<IloBoolVarArray>  IloBoolVarArray2;
typedef IloArray<IloNumVarArray2>  IloNumVarArray3;
typedef IloArray<IloRangeArray>  IloRangeArray2;
typedef IloArray<IloNumArray> TwoDMatrix;
typedef IloArray<IloExprArray> IloExprArray2;
typedef IloArray<IloExprArray2> IloExprArray3;
typedef IloArray<IloExprArray3> IloExprArray4;


#endif;
