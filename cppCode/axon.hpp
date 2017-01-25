#ifndef AXON_H
#define AXON_H

#include <iostream>
#include <fstream>

using namespace std;

class axon{
	private:
		double length, d1 ,r1, r2, lArea, R_A, C_M, R_a, C_m;
		double iNa, iK, iLeak, m, h, nk, v;
	public:
		axon(double vInit);
		double inf_m(double v);
		double tau_h(double v);
		double inf_h(double v);
		double I_na_a(double v, double m, double h);
		double alpha_nk(double v);
		double beta_nk(double v);
		double tau_nk(double v);
		double inf_nk(double v);
		double I_k(double v, double nk);

		void update(double I_in, double dt);
		void writeToFile(ofstream& file);
		void writeTau(ofstream& file);
		void writeIdensity(ofstream& file);

		void writeInfM(double vStart, double vEnd, double vStep);
		void writeTauH(double vStart, double vEnd, double vStep);
		void writeInfH(double vStart, double vEnd, double vStep);
		void writeTauNk(double vStart, double vEnd, double vStep);
		void writeInfNk(double vStart, double vEnd, double vStep);

		void writeTauAndInfNa_a(double vStart, double vEnd, double vStep);
		void writeTauAndInfK(double vStart, double vEnd, double vStep);
		void writeAllTauAndInf(double vStart, double vEnd, double vStep);

		double getM();
		double getH();
		double getNk();
		double getV();
		double getRa();

};

#endif