#include "TestingTheNumericalDiffusion.h"

void test(int N, double alpha, std::string convection) {
	std::cout << "N=" << N << " convection=" << convection << std::endl;
	TestingTheNumericalDiffusion TTND(1.0, 1.0, N, N, alpha);
	Data data;
	TTND.solve(data);
	data.save("N=" + std::to_string(N), convection);
}

void test(double alpha, std::string convection) {
	test(20, alpha, convection);
	test(40, alpha, convection);
	test(80, alpha, convection);
}

void test() {
	test(1.0, "UDS");
	test(0.05, "blend");
}

int main() {
	test();
	return 0;
}
