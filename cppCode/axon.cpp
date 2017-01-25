#include <math.h>
#include "axon.hpp"
#include "helper.hpp"
/**
v in Volt
**/
/**
I_na_a
**/

/**
inf_m = 1 / (1 + exp(0 - (v - -30 * pow(10,-3) / (5.5 * pow(10, -3))))
**/
double axon::inf_m(double v){
	return HHSigmoidVariable(v, -30 * pow(10, -3), 1, 5.5 * pow(10, -3));
}

// m changes instantly so no dm/dt

/**
tau_h = pow(10, -3) * 1.5 * exp((v - -40 * pow(10 , -3))/ (-33 * pow(10, -3))
**/
double axon::tau_h(double v){
	return pow(10, -3) * expTime(v, 1.5, -40 * pow(10,-3), -33 * pow(10, -3));
}

/**
inf_h = 1 / (1+ exp(0 - (v - -60*pow(10,-3)) / -5.8 * pow(10, -3)));
**/
double axon::inf_h(double v){
	return HHSigmoidVariable(v, -60 * pow(10, -3), 1, -5.8 * pow(10, -3));
}

double axon::I_na_a(double v, double m, double h){
	double condDensity = 240 * pow(10, -3);
	double erev = 55 * pow(10, -3);
	return condDensity * m * m * m * h * (erev -v) * lArea;
}

/**
x = (v - -25 * pow(10, -3)) / 10 * pow(10, -3);
	if(x == 0){
		alpha_nk = 10 * pow(10, -3);
	}else{
		alpha_nk = 10 * pow(10, -3) * x / (1 - exp(0 - x));
	}
**/
double axon::alpha_nk(double v){
	return hhExpLinearRate(v, -25 * pow(10,-3), 1.3 * pow(10,3), 10 * pow(10,-3));
}

/**
beta_nk = 1.69 * pow(10, 3) * exp((v- -35 * pow(10, -3))/ (-80 * pow(10, -3))) 
**/
double axon::beta_nk(double v){
	return hhExpRate(v, -35 * pow(10,-3), 1.69 * pow(10,3), -80 * pow(10,-3));
}

double axon::tau_nk(double v){
	return 1 / (alpha_nk(v) + beta_nk(v));
}

double axon::inf_nk(double v){
	return alpha_nk(v) / (alpha_nk(v) + beta_nk(v));
}

double axon::I_k(double v, double nk){
	double condDensity = 240 * pow(10, -3);
	double erev = -75 * pow(10, -3);
	return condDensity * nk * nk * nk * nk * (erev -v) * lArea;
}

double axon::getM(){
	return m;
}

double axon::getH(){
	return h;
}

double axon::getNk(){
	return nk;
}

double axon::getV(){
	return v;
}

double axon::getRa(){
	return R_a;
}

void axon::update(double iIn, double dt){
	double dvdt = (iIn - iNa - iK - iLeak) /C_m;
	double dhdt = (inf_h(v) - h) / tau_h(v);
	double dnkdt = (inf_nk(v) - nk) / tau_nk(v);

	m = inf_m(v);
	v = fwd_euler(v, dvdt, dt);
	h = fwd_euler(h, dvdt, dt);
	nk = fwd_euler(nk, dnkdt, dt);

	iNa = I_na_a(v, m, h);
	iK = I_k(v, nk);
	iLeak = I_leak(v, lArea);
}

axon::axon(double vInit){
	length = 0.0003; //(cm) 0.0003 cm = 3 micron
	d1 = 0.004; // (cm) 0.004 cm = 40 micrion diameter taken from proximal
	r1 = d1/2;
	r2 = 0.0011; // cm diameter taken from distal, r = d/ 2 = 0.0022 /2 
	R_A = 3;	// ohm * cm
	C_M = 0.000001; // (F) 0.000001 F / cm^2= 1 uF /cm^2
	lArea = PI * (r1 + r2) * sqrt((r1 - r2) * (r1 - r2) + length * length); // (cm^2)
	R_a = (4 * length * R_A) / (PI * d1 * d1); // (cm * ohm * cm / (cm * cm) = ohm)
	C_m = lArea * C_M; // (cm^2 * F/ cm^2 = F)
	v = vInit;
	m = inf_m(v);
	h = inf_h(v);
	nk = inf_nk(v);
	iNa = I_na_a(v, m, h);
	iK = I_k(v, nk);
	iLeak = I_leak(v, lArea);
}

void axon::writeToFile(ofstream& file){
	file << m << "\t";
	file << h << "\t";
	file << nk << "\n"; // wrong
}

void axon::writeTau(ofstream& file){
	file << tau_h(v) << "\t";
	file << tau_nk(v) << "\n";
}

void axon::writeIdensity(ofstream& file){
	file << iNa / lArea << "\t";
	file << iK / lArea << "\t";
	file << iLeak / lArea << "\n";
}

void axon::writeInfM(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/axon_inf_m.dat");
	file << "v\tinf_m\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_m(vTemp) << "\n";
	}
	file.close();

}

void axon::writeTauH(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/axon_tau_h.dat");
	file << "v\ttau_h\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_h(vTemp) << "\n";
	}
	file.close();
}

void axon::writeInfH(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/axon_inf_h.dat");
	file << "v\tinf_h\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_h(vTemp) << "\n";
	}
	file.close();
}

void axon::writeTauAndInfNa_a(double vStart, double vEnd, double vStep){
	writeInfM(vStart, vEnd, vStep);
	writeTauH(vStart, vEnd, vStep);
	writeInfH(vStart, vEnd, vStep);
}

void axon::writeTauNk(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/axon_tau_nk.dat");
	file << "v\ttau_nk\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_nk(vTemp) << "\n";
	}
	file.close();
}

void axon::writeInfNk(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/axon_inf_nk.dat");
	file << "v\tinf_nk\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_nk(vTemp) << "\n";
	}
	file.close();
}

void axon::writeTauAndInfK(double vStart, double vEnd, double vStep){
	writeTauNk(vStart, vEnd, vStep);
	writeInfNk(vStart, vEnd, vStep);
}

void axon::writeAllTauAndInf(double vStart, double vEnd, double vStep){
	writeTauAndInfNa_a(vStart, vEnd, vStep);
	writeTauAndInfK(vStart, vEnd, vStep);
}

