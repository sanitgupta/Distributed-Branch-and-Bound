#include <iostream>
#include <vector>
#include <string>
#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"


double z_opt = 2000000000; 
//int n,m;
//int x[65535];

void printMatrix(std::vector<std::vector<float> > &);
void printVector(std::vector<float>);
std::vector<float> solve(std::vector<std::vector<float> >, std::vector<float>, std::vector<float>);
//bool is_integer(float);
//float calculate_optimum(std::vector<operations_research::MPVariable*>, std::vector<float> , int);

namespace operations_research {
	bool is_integer(float value){
		if(floorf(value) - value < 0.00000001){
			return true;
		}
		return false;
	}

	float calculate_optimum(std::vector<MPVariable*> x, std::vector<float> &c, int n){
		float z = 0;
		for (int i = 0; i < n; i++){
			z += c[i]*(x[i]->solution_value());
		}
		return z;
	}

	std::vector<MPVariable*> RunTest(MPSolver::OptimizationProblemType optimization_problem_type, int n, int m, std::vector<std::vector<float> > &A, std::vector<float> &b, std::vector<float> &c) {
	  	MPSolver solver("LinearExample", optimization_problem_type);

		const double infinity = solver.infinity();

	    // Create the variables x and y.
	    std::vector<MPVariable*> variables{0};
	    variables.pop_back();
	    for(int i=0;i<n;i++){
	    	variables.push_back(solver.MakeNumVar(-infinity, infinity, "variables["+std::to_string(i)+"]"));
	    }

		MPObjective* const objective = solver.MutableObjective();

	    for(int i = 0; i < n; i++){    	
	    	objective->SetCoefficient(variables[i], c[i]);
	    }
	    objective->SetMinimization();

		std::vector<MPConstraint*> constraints(m);
	    for(int i=0;i<m;i++){
	        constraints[i] = solver.MakeRowConstraint(-infinity, b[i]);
	  		for(int j=0;j<n;j++){
	  			constraints[i]->SetCoefficient(variables[j], A[i][j]);
	  		}
	  	}

	// Call the solver and display the results.
	    solver.Solve();

 // 	for(int i = 0; i < n; i++){
 //    	printf("x%d = %.1f\n", i, variables[i]->solution_value());
	// }
		std::cout << calculate_optimum(variables,c,n) << std::endl; 
		return variables;
	}

  	std::vector<MPVariable*> RunExample(int n, int m, std::vector<std::vector<float> > &A, std::vector<float> &b, std::vector<float> &c) {
	    return RunTest(MPSolver::GLOP_LINEAR_PROGRAMMING,n,m,A,b,c);
  	}
  	void solve(int n, int m, std::vector<std::vector<float> > &A, std::vector<float> &b, std::vector<float> &c) {
	 std::vector<MPVariable*> solution_temp= operations_research::RunExample(n,m,A,b,c);
	 std::vector<std::vector<float> > A_temp_less = A;
	 std::vector<std::vector<float> > A_temp_more = A;
	 std::vector<float> b_temp_less = b;
	 std::vector<float> b_temp_more = b;

	 int i=0;
	 for(i=0; i<n;i++){
	 	if(!is_integer(solution_temp[i]->solution_value())){
	 		break;
	 	}
	 }
	 if(i==n){
	 	float z_temp = operations_research::calculate_optimum(solution_temp,c,n);
	 	if(z_opt > z_temp){
	 		z_opt= z_temp;
	 	}
	 	return;
	}

	std::vector<float> v(n);
 	v[i] = 1;
 	A_temp_less.push_back(v);
 	b_temp_less.push_back(floorf(solution_temp[i]->solution_value()));
 	solve(n,m+1,A_temp_less,b_temp_less,c);	

 	v[i] = -1;
 	A_temp_more.push_back(v);
 	b_temp_more.push_back(-ceilf(solution_temp[i]->solution_value()));
 	solve(n,m+1,A_temp_more,b_temp_more,c);	
 	return;
}
}


int main(){
	int n,m;
	std::cout << "Enter of number of variables";
	std::cin >> n;
	std::cout << "Enter of number of constraints";
	std::cin >> m;
	std::vector<std::vector<float> > A(m, std::vector<float>(n));
	std::vector<float> b(m);
	std::vector<float> c(n);

	std::cout << "Enter cost function of format c1 enter c2...cn enter";
	for (size_t i = 0; i < n; i++){
		std::cin >> c[i];
	}

	std::cout<< "Enter constraints in format a11 enter a12 enter a13....an enter b1 enter repeat for m constraints";

	for (size_t i = 0; i < m; i++){
		for (size_t j = 0; j< n; j++){
			std::cin >> A[i][j];
		}
		std::cin >> b[i];
	}

	// printMatrix(A);
	// printVector(b);
	// printVector(c);
	std::cout << std::endl;
	operations_research::solve(n,m,A,b,c);
	std::cout << z_opt << std::endl;

}

void printMatrix(std::vector<std::vector<float> > &A) {
	for (auto vec1 : A) {
		for (auto elem : vec1) {
			std::cout << elem << " ";
		}
		std::cout << std::endl;
	}
}

void printVector(std::vector<float>  A) {
	for (auto e : A) {
			std::cout << e << " ";
	}
	std::cout << std::endl;
}
