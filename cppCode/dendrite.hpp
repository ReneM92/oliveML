#ifndef DENDRITE_H
#define DENDRITE_H

#include "caModel.hpp"
#include <iostream>
#include <fstream>

using namespace std;

class dendrite{
	private:
		double length, d, lArea, R_A, C_M, R_a, C_m;
		double k, l, r, z, n;
		double iCal, iCah, iKca, iH, iLeak, v; 
		caConcentrationModel caModel;

	public:
		dendrite(double vInit);
		double tau_k();
		double inf_k(double v);
		double tau_l(double v);
		double inf_l(double v);
		double I_cal(double v, double k, double l);
		double alpha_r(double v);
		double beta_r(double v);
		double tau_r(double v);
		double inf_r(double v);
		double I_cah(double v, double r);
		double alpha_z(double caNorm);
		double beta_z();
		double tau_z(double caNorm);
		double inf_z(double caNorm);
		double I_kca(double v, double z);
		double tau_n(double v);
		double inf_n(double v);
		double I_h(double v);

		void update(double I_in, double dt);
		void writeToFile(ofstream& file);
		void writeTau(ofstream& file);
		void writeIdensity(ofstream& file);

		void writeTauK(double vStart, double vEnd, double vStep);
		void writeInfK(double vStart, double vEnd, double vStep);
		void writeTauL(double vStart, double vEnd, double vStep);
		void writeInfL(double vStart, double vEnd, double vStep);
		void writeTauR(double vStart, double vEnd, double vStep);
		void writeInfR(double vStart, double vEnd, double vStep);
		void writeTauN(double vStart, double vEnd, double vStep);
		void writeInfN(double vStart, double vEnd, double vStep);

		void writeTauAndInfCal(double vStart, double vEnd, double vStep);
		void writeTauAndInfCah(double vStart, double vEnd, double vStep);
		void writeTauAndInfH(double vStart, double vEnd, double vStep);
		void writeAllTauAndInf(double vStart, double vEnd, double vStep);

		double getK();
		double getL();
		double getR();
		double getZ();
		double getN();
		double getCaConc();
		double getV();
		double getRa();

};

#endif