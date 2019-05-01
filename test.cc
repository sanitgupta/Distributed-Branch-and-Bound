#include <iostream>
#include <vector>
#include <string>
#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"


double z_opt = 2000000000; 
int n,m;
int x[65535];

void printMatrix(std::vector<std::vector<float> > &);
void printVector(std::vector<float>);
std::vector<float> solve(std::vector<std::vector<float> >, std::vector<float>, std::vector<float>);

namespace operations_research {
void RunTest(MPSolver::OptimizationProblemType optimization_problem_type, int n, int m, std::vector<std::vector<float> > &A, std::vector<float> &b, std::vector<float> &c) {
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

 	for(int i = 0; i < n; i++){
    	printf("x = %.1f\n", variables[i]->solution_value());
	}
    
}

  void RunExample(int n, int m, std::vector<std::vector<float> > &A, std::vector<float> &b, std::vector<float> &c) {
    RunTest(MPSolver::GLOP_LINEAR_PROGRAMMING,n,m,A,b,c);
  }
}

int main(){
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
	operations_research::RunExample(n,m,A,b,c);
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
