#include "utils.hpp"

using namespace std;
using namespace arma;

int main (){ 

	mat A = mat(4, 4, fill::eye);
	A(1,2) = -0.7;
	A(2,1) = -0.7;
	A(0,3) = 0.5;
	A(3,0) = 0.5;

	int k,l = 0;

	double maxvalue = max_offdiag_symmetric(A, k, l);


	cout << "Max value: "<<maxvalue << endl;
	cout << "In row-column position (counting from 0): "<< k <<" "<< l << endl;

	return 0;

}
