#include <iostream>
#include <mpi.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <mpi.h>
#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"
 
using namespace std;

#define n 2
#define m 3

bool double_is_int(double num) {
   double absolute = abs(num);
   return absolute == floor(absolute);
}

int is_array_int(double LOCAL_VAL[n]) {
	int i = 0;
	for(i=0; i<n; i++)
	{
		if(double_is_int(LOCAL_VAL[i]))
			continue;
		else
			break;
	}
	return i;
}

namespace operations_research { ///
	void MIPSolver(double A[m][n], double B[m], double C[n], double A_range[2*n + n+1], double OPT[1], double OPT_VAL[n])
	{
		MPSolver::OptimizationProblemType opt_prob_type = MPSolver::GLOP_LINEAR_PROGRAMMING;////
	    MPSolver solver("Glop", opt_prob_type);

	    MPVariable* variables[n];

	    for(int i=0; i<n; i++)
	    {
	        variables[i] = solver.MakeNumVar(A_range[i], A_range[i+n], "var");
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

	    MPObjective* objective = solver.MutableObjective();
	   	for(int i=0; i<n; i++)
	    {
	        //variables[i] = solver.MakeNumVar(min_max[i], min_max[i+n], "var");
	        objective->SetCoefficient(variables[i], C[i]);
	    }

		objective->SetMaximization(); //Can be changed

	    printf("\nNumber of variables = %d", solver.NumVariables());

	    printf("\nNumber of constraints = %d", solver.NumConstraints());

	    solver.Solve(); //return

	    OPT[0] = objective->Value();

	    for(int i=0; i<n; i++) {
	        OPT_VAL[i] = variables[i]->solution_value();
	    }
	}
}


int main(int argc, char **argv)
{

	double A[m][n], B[m], C[n], LOCAL_VAL[n], GLOBAL_VAL[n], LOCAL_OPT[1];
    double A_range[2*n + n+1]; //Initialize with zeros and infinity
    int first_double_index, color;

	//Initializing all variables in LP
    C[0] = 3, C[1] = 4;
    A[0][0] = 1, A[0][1] = 2, A[1][0] = -3, A[1][1] = 1, A[2][0] = 1, A[2][1] = -1;
    B[0] = 14, B[1] = 0, B[2] = 2;

    //3x + 4y
    // x + 2y <= 14.
  	// -3x + y <= 0.
	// x - y <= 2.

    const double infinity = numeric_limits<double>::max();


    const int int_infinity = numeric_limits<int>::max();
    double GLOBAL_OPT = -infinity;
    LOCAL_OPT[0] = -infinity;

    //Initializing X Range
    for(int i=0; i<n; i++)
    {
        A_range[i] = 0;
        A_range[i+n] = infinity;
    }

    int pid, total_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &total_proc);

	if (pid == 0)
		color = 56;
	else
		color = log2(pid); // Determine color based on depth

	MPI_Comm depth_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, pid, &depth_comm);

	// Process 0 is not required.
    if (pid!=0)
    {

	// Split the communicator based on the depth and use the
	// original PID for ordering

	

   if(pid!=1)
    {
    	MPI_Recv( A_range, 2*n + (n+1), MPI_DOUBLE, pid/2, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	GLOBAL_OPT = A_range[2*n];
    	for (int i=0; i<n; i++)
    		GLOBAL_VAL[i] = A_range[2*n + 1 + i];
    }

    // printf("%d \n", color);
    if(A_range[0] != -1.0)
    {
    	operations_research::MIPSolver(A, B, C, A_range, LOCAL_OPT, LOCAL_VAL);

	    if(LOCAL_OPT[0] > GLOBAL_OPT)     //The case when Local Value is greater than current global.
	    {
	    	first_double_index = is_array_int(LOCAL_VAL);
	    	if (first_double_index == n)
	    	{
	    		GLOBAL_OPT = LOCAL_OPT[0];
	    		A_range[2*n] = GLOBAL_OPT;
	    		for (int i=0; i<n; i++)
	    		{
	    			GLOBAL_VAL[i] = LOCAL_VAL[i];
    				A_range[2*n + 1 + i] = GLOBAL_VAL[i];
	    		}
	    		A_range[0] = -1;  //Disable further computation
	    	}
	    	else
	    	{
	    		int val = LOCAL_VAL[first_double_index];
	    		A_range[first_double_index] = val;  //Integer Part only
	    	}

	    }
	    else
	    	A_range[0] = -1;	//Disable further computation
	}


    struct { 
        double val; 
        int   rank; 
    } G_OPT_Receiver;

    struct { 
        double val; 
        int   rank; 
    } G_OPT_Comparer;

    G_OPT_Comparer.val = GLOBAL_OPT;
    MPI_Comm_rank(depth_comm, &G_OPT_Comparer.rank);;  
    
    //MPI ALL REDUCE with comm for same depth
	MPI_Allreduce( &G_OPT_Comparer, &G_OPT_Receiver, 1, MPI_DOUBLE_INT, MPI_MAXLOC, depth_comm);
	printf("%f %d \n", G_OPT_Receiver.val, G_OPT_Receiver.rank);
	MPI_Bcast( GLOBAL_VAL, n, MPI_DOUBLE, G_OPT_Receiver.rank, depth_comm);

	GLOBAL_OPT = G_OPT_Receiver.val;
	A_range[2*n] = GLOBAL_OPT;
	for (int i=0; i<n; i++)
		A_range[2*n + 1 + i] = GLOBAL_VAL[i];

	if (2*pid+1 <= total_proc)
		{
			if(A_range[0] == -1)
			{
				MPI_Send( A_range, 2*n + n+1, MPI_DOUBLE, 2*pid, 0, MPI_COMM_WORLD);
    			MPI_Send( A_range, 2*n + n+1, MPI_DOUBLE, 2*pid+1, 0, MPI_COMM_WORLD);
    		}
    		else
    		{
    			MPI_Send( A_range, 2*n + n+1, MPI_DOUBLE, 2*pid, 0, MPI_COMM_WORLD);
    			A_range[first_double_index + n] = A_range[first_double_index] + 1;
    			A_range[first_double_index] = 0;
    			MPI_Send( A_range, 2*n + n+1, MPI_DOUBLE, 2*pid+1, 0, MPI_COMM_WORLD);
    		}
    	}
    	
    }

    MPI_Comm_free(&depth_comm);
    MPI_Finalize();
}