#include "utils_4.hpp"

using namespace std;
using namespace arma;

int main (){ 

	int N =6;
	double h = 1./N+1; 
	double a = -1./(pow(h, 2));
	double d = 2./(pow(h, 2));
	int k;
	int l;

	mat A = create_tridiagonal(N, a, d);
	mat R = mat(N, N, fill::eye);
	double maxvalue = max_offdiag_symmetric(A, k, l);
	