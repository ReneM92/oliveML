#ifndef SOMA_H
#define SOMA_H

#include <iostream>
#include <fstream>

using namespace std;

class soma{
	private:
		double length, d, lArea, C_M, R_A, C_m, R_a;
		double m, h, n, nk, iNa, iKdr, iK, iLeak, v;
	public:
		double inf_m(double v);
		double tau_h(double v);
		double inf_h(double v);
		double I_na_s(double v, double m, double h);
		double tau_n(double v);
		double inf_n(double v);
		double I_kdr(double v, double n);
		double alpha_nk(double v);
		double beta_nk(double v);
		double tau_nk(double v);
		double inf_nk(double v);
		double I_k(double v, double nk);
		soma(double vInit);

		void update(double I_in, double dt);
		void writeToFile(ofstream& file);
		void writeTau(ofstream& file);
		void writeIdensity(ofstream& file);

		void writeInfM(double vStart, double vEnd, double vStep);
		void writeTauH(double vStart, double vEnd, double vStep);
		void writeInfH(double vStart, double vEnd, double vStep);
		void writeTauN(double vStart, double vEnd, double vStep);
		void writeInfN(double vStart, double vEnd, double vStep);
		void writeTauNk(double vStart, double vEnd, double vStep);
		void writeInfNk(double vStart, double vEnd, double vStep);

		void writeTauAndInfNa_s(double vStart, double vEnd, double vStep);
		void writeTauAndInfKdr(double vStart, double vEnd, double vStep);
		void writeTauAndInfK(double vStart, double vEnd, double vStep);
		void writeAllTauAndInf(double vStart, double vEnd, double vStep);

		double getM();
		double getH();
		double getN();
		double getNK();
		double getV();
		double getRa();

	};

#endif