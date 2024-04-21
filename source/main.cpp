#include "matrix.h"

// example

int main () {
	m::Matrix A, B;
	A.resize(2, 3);
	B.resize(3, 2);
	A.enter();
	B.enter();
	A.multiply(B);
	A.print();
	return 0;
}
