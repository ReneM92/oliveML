#include <math.h>
#include "soma.hpp"
#include "helper.hpp"
/**
v in Volt
**/
/**
I_na_s
**/

/**
inf_m = 1 / (1 + exp(0 - (v - -30 * pow(10,-3) / (5.5 * pow(10, -3))))
**/
double soma::inf_m(double v){
	return HHSigmoidVariable(v, -30 * pow(10,-3), 1, 5.5 * pow(10,-3));
}

// m changes instantly so no dm/dt

/**
tau_h = pow(10, -3) * 1.5 * exp((v - -40 * pow(10 , -3))/ (-33 * pow(10, -3))
**/
double soma::tau_h(double v){
	return pow(10, -3) * expTime(v, 3, -40 * pow(10,-3), -33 * pow(10,-3));
}

/**
inf_h = 1 / (1+ exp(0 - (v - -60*pow(10,-3)) / -5.8 * pow(10, -3)));
**/
double soma::inf_h(double v){
	return HHSigmoidVariable(v, -70 * pow(10,-3), 1, -5.8 * pow(10,-3));
}

double soma::I_na_s(double v, double m, double h){
	double condDensity = 150 * pow(10,-3);
	double erev = 55 * pow(10,-3);
	return condDensity * m * m * m * h * (erev -v) * lArea;
}

double soma::tau_n(double v){
	double t;
	v = v * pow(10, 3);
	t = (5 + (47 * exp(-1 * (-50 - v) / 900)));
	return t * pow(10, -3);
}

/**
inf_n = 1 / (1+ exp(0 - (v - -3 * pow(10, -3)) / 10 * pow(10, -3)))
**/
double soma::inf_n(double v){
	return HHSigmoidVariable(v, -3 * pow(10,-3), 1, 10 * pow(10, -3));
}

double soma::I_kdr(double v, double n){
	double condDensity = 9 * pow(10, -3);
	double erev = -75 * pow(10, -3);
	return condDensity * n * n * n * n * (erev -v) * lArea;
}

/**
k
**/
double soma::alpha_nk(double v){
	return hhExpLinearRate(v, -25 * pow(10, -3), 1.3 * pow(10 ,3), 10 * pow(10,-3));
}

double soma::beta_nk(double v){
	return hhExpRate(v, -35 * pow(10, -3), 1.69 * pow(10, 3), -80 * pow(10, -3));
}

double soma::tau_nk(double v){
	return 1 / (alpha_nk(v) + beta_nk(v));
}

double soma::inf_nk(double v){
	return alpha_nk(v) / (alpha_nk(v) + beta_nk(v));
}

double soma::I_k(double v, double nk){
	double condDensity = 5 * pow(10, -3);
	double erev = -75 * pow(10, -3);
	return condDensity * nk * nk * nk * nk * (erev -v) * lArea;
}

double soma::getM(){
	return m;
}

double soma::getH(){
	return h;
}

double soma::getN(){
	return n;
}
		
double soma::getNK(){
	return nk;
}
double soma::getV(){
	return v;
}

double soma::getRa(){
	return R_a;
}

void soma::update(double iIn, double dt){
	double dvdt = (iIn - iNa - iKdr - iK - iLeak)/ C_m;
	double dhdt = (inf_h(v) - h) / tau_h(v);
	double dndt = (inf_n(v) - n) / tau_h(v);
	double dnkdt = (inf_nk(v) - nk)/ tau_nk(v); 

	m = inf_m(v);
	v = fwd_euler(v, dvdt, dt);
	h = fwd_euler(h, dhdt, dt);
	n = fwd_euler(n, dndt, dt);
	dnkdt = fwd_euler(nk, dnkdt, dt);

	iNa = I_na_s(v, m, h);
	iKdr = I_kdr(v, n);
	iK = I_k(v, nk);
	iLeak = I_leak(v, lArea);
}

void soma::writeToFile(ofstream& file){
	file << m << "\t";
	file << h << "\t"; // wrong
	file << n << "\t";
	file << nk << "\n"; // wrong
}

void soma::writeTau(ofstream& file){
	file << tau_h(v) << "\t";
	file << tau_n(v) << "\t";
	file << tau_nk(v) << "\n";
}

void soma::writeIdensity(ofstream& file){
	file << iNa / lArea << "\t";
	file << iKdr / lArea << "\t";
	file << iK / lArea << "\t";
	file << iLeak / lArea << "\n";
}

soma::soma(double vInit){
	v = vInit;
	length = 0.002; // (cm) 0.002 cm = 20 microns 
	d = 0.004; //(cm) 0.004 cm = 40 microns
	R_A = 0.1; // ohm * cm
	C_M = 0.000001; // (F) 0.000001 F / cm^2= 1 uF /cm^2
	lArea = PI * d * length; // (cm^2)
	R_a = (4 * length * R_A) / (PI * d * d); //(cm * ohm * cm / (cm * cm) = ohm)
	C_m = lArea * C_M; // // (cm^2 * F/ cm^2 = F)
	m = inf_m(v);
	h = inf_h(v);
	n = inf_n(v);
	nk = inf_nk(v);
	iNa = I_na_s(v, m, h);
	iKdr = I_kdr(v, n);
	iK = I_k(v, nk);
	iLeak = I_leak(v, lArea);

	/*
	printf("init Soma\n");
	printf("v: %.6f\n",v);
	printf("l: %.6f\n",l);
	printf("d: %.6f\n",d);
	printf("R_A: %.6f\n", R_A);
	printf("C_M: %.6f\n", C_M);
	printf("lArea %.6f\n", lArea);
	printf("R_a %.6f\n", R_a);
	printf("C_m %.6f\n", C_m);
	printf("m %.6f\n", m);
	printf("h %.6f\n", h);
	printf("n %.6f\n", n);
	printf("nk %.6f\n", nk);
	printf("iNa %.6f\n", iNa);
	printf("iKdr %.6f\n", iKdr);
	printf("iK %.6f\n", iK);
	printf("iLeak %.6f\n", iLeak);
	*/
}

void soma::writeAllTauAndInf(double vStart, double vEnd, double vStep){
	writeTauAndInfNa_s(vStart, vEnd, vStep);
	writeTauAndInfKdr(vStart, vEnd, vStep);
	writeTauAndInfK(vStart, vEnd, vStep);
}

void soma::writeInfM(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/soma_inf_m.dat");
	file << "v\tinf_m\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_m(vTemp) << "\n";
	}
	file.close();

}

void soma::writeTauH(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/soma_tau_h.dat");
	file << "v\ttau_h\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_h(vTemp) << "\n";
	}
	file.close();
}

void soma::writeInfH(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/soma_inf_h.dat");
	file << "v\tinf_h\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_h(vTemp) << "\n";
	}
	file.close();
}

void soma::writeTauAndInfNa_s(double vStart, double vEnd, double vStep){
	writeInfM(vStart, vEnd, vStep);
	writeTauH(vStart, vEnd, vStep);
	writeInfH(vStart, vEnd, vStep);
}

void soma::writeTauN(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/soma_tau_n.dat");
	file << "v\ttau_n\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_n(vTemp) << "\n";
	}
	file.close();
}

void soma::writeInfN(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/soma_inf_n.dat");
	file << "v\tinf_n\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_n(vTemp) << "\n";
	}
	file.close();
}

void soma::writeTauAndInfKdr(double vStart, double vEnd, double vStep){
	writeTauN(vStart, vEnd, vStep);
	writeInfN(vStart, vEnd, vStep);
}

void soma::writeTauNk(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/soma_tau_nk.dat");
	file << "v\ttau_nk\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_nk(vTemp) << "\n";
	}
	file.close();
}

void soma::writeInfNk(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/soma_inf_nk.dat");
	file << "v\tinf_nk\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_nk(vTemp) << "\n";
	}
	file.close();
}

void soma::writeTauAndInfK(double vStart, double vEnd, double vStep){
	writeTauNk(vStart, vEnd, vStep);
	writeInfNk(vStart, vEnd, vStep);
}
