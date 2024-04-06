#pragma once
#include "Data.h"

class TestingTheNumericalDiffusion {
	double u_x;
	double u_y;
	int Nx;
	int Ny;
	double dx;
	double dy;
	double alpha;
	void initialize(Data& data, std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P, std::vector<std::vector<double>>& A, std::vector<double>& b) {
		data.x.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			data.x[i] = (i + 0.5) * dx;
		}
		data.y.resize(Ny);
		for (int i = 0; i < Ny; i++) {
			data.y[i] = (i + 0.5) * dy;
		}
		data.phi.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			data.phi[i].resize(Ny);
		}
		A_W.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_W[i].resize(Ny);
		}
		A_S.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_S[i].resize(Ny);
		}
		A_P.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_P[i].resize(Ny);
		}
		A_N.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_N[i].resize(Ny);
		}
		A_E.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_E[i].resize(Ny);
		}
		Q_P.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			Q_P[i].resize(Ny);
		}
		A.resize(Nx * Ny);
		for (int i = 0; i < Nx * Ny; i++) {
			A[i].resize(Nx * Ny);
		}
		b.resize(Nx * Ny);
	}
	void discretize(Data& data, std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P, std::vector<std::vector<double>>& A, std::vector<double>& b) {
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				if (i == Nx - 1) {
					A_P[i][j] += u_x * dy;
				}
				else {
					A_P[i][j] += (1 + alpha) * u_x * dy / 2;
					A_E[i][j] += (1 - alpha) * u_x * dy / 2;
				}
				if (i == 0) {
					Q_P[i][j] += u_x * dy * (data.y[j] < 0.1 ? 0.0 : 1.0);
				}
				else {
					A_W[i][j] += -(1 + alpha) * u_x * dy / 2;
					A_P[i][j] += -(1 - alpha) * u_x * dy / 2;
				}
				if (j == Ny - 1) {
					A_P[i][j] += u_y * dx;
				}
				else {
					A_P[i][j] += (1 + alpha) * u_y * dx / 2;
					A_N[i][j] += (1 - alpha) * u_y * dx / 2;
				}
				if (j > 0) {
					A_S[i][j] += -(1 + alpha) * u_y * dx / 2;
					A_P[i][j] += -(1 - alpha) * u_y * dx / 2;
				}
			}
		}
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				if (i > 0) {
					A[i * Ny + j][(i - 1) * Ny + j] = A_W[i][j];
				}
				if (j > 0) {
					A[i * Ny + j][i * Ny + (j - 1)] = A_S[i][j];
				}
				A[i * Ny + j][i * Ny + j] = A_P[i][j];
				if (j < Ny - 1) {
					A[i * Ny + j][i * Ny + (j + 1)] = A_N[i][j];
				}
				if (i < Nx - 1) {
					A[i * Ny + j][(i + 1) * Ny + j] = A_E[i][j];
				}
				b[i * Ny + j] = Q_P[i][j];
			}
		}
	}
	void GEPP(std::vector<std::vector<double>>& A, std::vector<double>& b) {
		for (int i = 0; i < A.size() - 1; i++) {
			int index = i;
			for (int j = i + 1; j < A.size(); j++) {
				index = abs(A[index][i]) < abs(A[j][i]) ? j : index;
			}
			for (int j = i; j < A.size(); j++) {
				std::swap(A[i][j], A[index][j]);
			}
			std::swap(b[i], b[index]);
			for (int j = i + 1; j < A.size(); j++) {
				for (int k = i + 1; k < A.size(); k++) {
					A[j][k] -= A[j][i] / A[i][i] * A[i][k];
				}
				b[j] -= A[j][i] / A[i][i] * b[i];
			}
		}
		for (int i = A.size() - 1; i >= 0; i--) {
			for (int j = i + 1; j < A.size(); j++) {
				b[i] -= A[i][j] * b[j];
			}
			b[i] /= A[i][i];
		}
	}
	void save(Data& data, std::vector<double>& b) {
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				data.phi[i][j] = b[i * Ny + j];
			}
		}
	}
public:
	TestingTheNumericalDiffusion(double u_x, double u_y, int Nx, int Ny, double alpha) :u_x(u_x), u_y(u_y), Nx(Nx), Ny(Ny), dx(1.0 / Nx), dy(1.0 / Ny), alpha(alpha) {}
	void solve(Data& data) {
		std::vector<std::vector<double>> A_W, A_S, A_P, A_N, A_E, Q_P, A;
		std::vector<double> b;
		initialize(data, A_W, A_S, A_P, A_N, A_E, Q_P, A, b);
		discretize(data, A_W, A_S, A_P, A_N, A_E, Q_P, A, b);
		GEPP(A, b);
		save(data, b);
	}
};
