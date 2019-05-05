#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include<stdio.h>
#include <chrono> ////


#define n 387
#define m 388


double z_opt;
double x_opt[n];


using namespace std;
using namespace std::chrono;

namespace operations_research { ///
	void Simplex(double A[m][n], double B[m], double C[n], double min_max[2*n], double& OPT, double* OPT_VAL)
	{
		MPSolver::OptimizationProblemType opt_prob_type = MPSolver::GLOP_LINEAR_PROGRAMMING;////
	    MPSolver solver("Glop", opt_prob_type);

	    MPVariable* variables[n];

	    for(int i=0; i<n; i++)
	    {
	        variables[i] = solver.MakeNumVar(min_max[i], min_max[i+n], to_string(i));
	        //variables[i] = solver.MakeNumVar(0, 2, to_string(i));
	        //objective->SetCoefficient(variables[i], C[i]);
	    }

	   
	    int infinity = solver.infinity(); ///
	    MPConstraint* constraints[m];
	    for(int i = 0; i < m; i++)
	    {
	        constraints[i] = solver.MakeRowConstraint(-infinity, B[i]);
	        for(int j=0; j<n; j++) {
	        	constraints[i]->SetCoefficient(variables[j], A[i][j]);
	        }
	    }


	    printf("\nNumber of variables = %d", solver.NumVariables());

	    printf("\nNumber of constraints = %d", solver.NumConstraints());

	    MPObjective* objective = solver.MutableObjective();
	   	for(int i=0; i<n; i++)
	    {
	        //variables[i] = solver.MakeNumVar(min_max[i], min_max[i+n], "var");
	        objective->SetCoefficient(variables[i], C[i]);
	    }

		objective->SetMaximization(); //Can be changed


	    solver.Solve(); //return


	    OPT = objective->Value();

	    for(int i=0; i<n; i++) {
	        OPT_VAL[i] = variables[i]->solution_value();
	    }
	}
}


void BB(double A[m][n], double B[m], double C[n], double min_max[2*n]) {

	double loc_opt;
	double loc_val[n];

	operations_research::Simplex(A, B, C, min_max, loc_opt, loc_val);
	//if (p == -1) return; /////infeasibility condition we don't know
	if (z_opt >= loc_opt) return;



	int i = 0;
	for(i = 0; i < n; i++)
	    if(!(floorf(loc_val[i])==ceilf(loc_val[i])))
	        break;

	if(i==n){
	    if(loc_opt > z_opt){
	      z_opt = loc_opt;
		  copy(begin(loc_val), end(loc_val), begin(x_opt));
	    }
	    return;
	}
	
	double temp = min_max[i+n];
	min_max[i+n] = floorf(loc_val[i]);
	BB(A, B, C, min_max);

	min_max[i] = ceilf(loc_val[i]);
	min_max[i+n] = temp;
	BB(A, B, C, min_max);
	
	return;
}



int main(int argc, char** argv) {

	double A[m][n], B[m], C[n];
    double min_max[2*n]; //Initialize with zeros and infinity
	int K,N;
	cin>>K>>N;
	B[0] = K;
	for(int i=0;i<n;i++)
		cin>>C[i]>>A[0][i];
	for(int j=1;j<m;j++){
		B[j] = 1;
		for(int k=0;k<n;k++){
		if(j==k+1) A[j][k] = 1;
		else A[j][k] = 0;
	}
	}
    /*A[0][0] = 1;
    A[0][1] = 2;
    A[1][0] = 2;
    A[1][1] = 1;

    B[0] = 5;
    B[1] = 5;

    C[0] = 1;
    C[1] = 1;*/


/*
    MPSolver::OptimizationProblemType opt_prob_type = MPSolver::GLOP_LINEAR_PROGRAMMING;////
    MPSolver solver("MIPProblem", opt_prob_type);
*/
    const double inf = numeric_limits<double>::max();


    const int int_infinity = numeric_limits<int>::max();

    z_opt = 0;

    for(int i=0; i<n; i++)
    {
        min_max[i] = 0;//-infinity;
        min_max[i+n] = inf;
    }

    //operations_research::Simplex(A, B, C, min_max, z_opt, x_opt);

auto start = high_resolution_clock::now();
    BB(A, B, C, min_max);

auto stop = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(stop - start);
    cout<<endl;
    cout<<z_opt<<endl;
    /*for(int i = 0; i < n; i++)
    	cout<<x_opt[i]<<endl;*/
cout << duration.count() << endl;

  }
