#include <math.h>
#include "dendrite.hpp"
#include "helper.hpp"

/**
v in Volt
**/
/**
cal 
**/
double dendrite::tau_k(){
	return 1 * pow(10, -3);
}

/**
inf_k = 1 / (1+ exp(0 - (v - -61 * pow(10, -3)) / 4.2 * pow(10, -3)));
**/
double dendrite::inf_k(double v){
	return HHSigmoidVariable(v, -61 * pow(10,-3), 1, 4.2 * pow(10,-3));
}

/*
on Start: k = inf_k(vInit)
dk/dt = (inf_k - k) / tau_k
*/

/**
v = v * pow(10,3); change v to same magnitude as other variables
t = ((20 * exp((v + 160) / 30) / (1 + exp((v + 84) / 7.3))) + 35) * pow(10, -3);
**/
double dendrite::tau_l(double v){
	double t;
	v = v * pow(10,3);
	t = ((20 * exp((v + 160) / 30) / (1 + exp((v + 84) / 7.3))) + 35);
	return t * pow(10,-3); 
}

/**
inf_l = 1 / (1+ exp(0 - (v - - 85.5 * pow(10, -3) / (-8.5 * pow(10, -3))));
**/
double dendrite::inf_l(double v){
	return HHSigmoidVariable(v, -85.5 * pow(10,-3), 1, -8.5 * pow(10,-3));
}

/*
on Start: l = inf_l(vInit)
dl/dt = (inf_l - l) / tau_l
*/
double dendrite::I_cal(double v, double k, double l){
	double condDensity = 5.5 * pow(10,-3);
	double erev = 120 * pow(10,-3);
	return condDensity * k * k * k * l * (erev -v) * lArea;
}

/**
cah
**/
/**
alpha_r = 1.7 * pow(10, 3) / (1 + exp(0 - (v - 5 * pow(10, -3))/(13.9 * pow(10,-3))))
**/
double dendrite::alpha_r(double v){
	return hhSigmoidRate(v , 5 * pow(10, -3), 1.7 * pow(10,3), 13.9 * pow(10,-3));
	//return hhSigmoidRate(v , 5, 1.7, 13.9);
}

/**
x = (v - - 8.5 * pow(10, -3)) / (-5 * pow(10, -3));
	if(x == 0){
		beta_r = 0.1 * pow(10, 3);
	}else{
		beta_r = 0.1 * pow(10, 3) * x / (1 - exp(0 - x));
	}
**/
double dendrite::beta_r(double v){
	return hhExpLinearRate(v, -8.5 * pow(10,-3), 0.1 * pow(10,3), -5 * pow(10,-3));
}

double dendrite::tau_r(double v){
	return 5 / (alpha_r(v) + beta_r(v));
}

double dendrite::inf_r(double v){
	return alpha_r(v) / (alpha_r(v) + beta_r(v));
}

/*
on Start: r = inf_r(vInit)
dr/dt = (inf_r - r) / tau_r
*/

double dendrite::I_cah(double v, double r){
	double condDensity = 4.5 * pow(10,-3);
	double erev = 120 * pow(10, -3);
	//cout << "gDensity: " << condDensity * r * r << "\n";
	return condDensity * r * r * (erev -v) * lArea;
}

/**
kca
**/

double dendrite::alpha_z(double caNorm){
	return min(0.00002 * caNorm, 0.01);
}

double dendrite::beta_z(){
	return 0.015;
}

double dendrite::tau_z(double caNorm){
	return 1 * pow(10, -3)/ (alpha_z(caNorm) + beta_z());
}

double dendrite::inf_z(double caNorm){
	return alpha_z(caNorm) / (alpha_z(caNorm) + beta_z());
}

/*
on Start: z = inf_z(vInit)
dz/dt = (inf_z - z) / tau_z
*/

double dendrite::I_kca(double v, double z){
	double condDensity = 45 * pow(10,-3);
	double erev = -75 * pow(10, -3);
	return condDensity * z * (erev -v) * lArea;
}

double dendrite::tau_n(double v){
	double t;
	v = v * pow(10,3);
	t = 1 /(exp(-0.086 * v - 14.6) + exp(0.070 * v - 1.87));
	return t * pow(10, -3); 
}

/**
inf_n = 1 / (1+ exp(0 - (v - -80 * pow(10, -3)) / -4 * pow(10, -3)))
**/
double dendrite::inf_n(double v){
	return HHSigmoidVariable(v, -80 * pow(10,-3), 1, -4 * pow(10,-3));
}

/*
on Start: n = inf_n(vInit)
dn/dt = (inf_n - n) / tau_n
*/

double dendrite::I_h(double v){
	double condDensity = 0.12 * pow(10,-3);
	double erev = -43 * pow(10,-3);
	return condDensity * n * (erev -v) * lArea;
}

double dendrite::getK(){
	return k;
}

double dendrite::getL(){
	return l;
}

double dendrite::getR(){
	return r;
}

double dendrite::getZ(){
	return z;
}

double dendrite::getN(){
	return n;
}

double dendrite::getCaConc(){
	return caModel.getConcentration();
}

double dendrite::getV(){
	return v;
}

double dendrite::getRa(){
	return R_a;
}

void dendrite::update(double iIn, double dt){
	//printf("start update dendrite\n");
	//cout << "vOld: " << v << "\n";
	//cout << iIn << " " << iCal << " " << iCah << " " << iKca << " " << iH << " " << iLeak << "\n";
	double dvdt = (iIn - iCal - iCah - iKca - iH - iLeak)/C_m;
	double dkdt = (inf_k(v) - k) / tau_k();
	double dldt = (inf_l(v) - l) / tau_l(v);
	double drdt = (inf_r(v) - r) / tau_r(v);
	double dzdt = (inf_z(getCaConc()) - z) / tau_z(v);
	double dndt = (inf_n(v) - n) / tau_n(n);

	//cout << "dvdt: " << dvdt << "\n";
	//printf("calculated derivatis dendrite\n");
	//cout << "dt: " << dt << "\n";
	//cout << "fwd_euler: " << fwd_euler(v, dvdt, dt) << "\n";
	v = fwd_euler(v, dvdt, dt);
	k = fwd_euler(k, dkdt, dt);
	l = fwd_euler(l, dldt, dt);
	z = fwd_euler(z, dzdt, dt);
	n = fwd_euler(n, dndt, dt);
	caModel.update(iKca, lArea, dt);
	//cout << "vNew: " << v << "\n";
	iCal = I_cal(v, k , l);
	iCah = I_cah(v, r);
	iKca = I_kca(v, z);
	iH = I_h(v);
	iLeak = I_leak(v, lArea);

}

void dendrite::writeToFile(ofstream& file){
	file << k << "\t";
	file << l << "\t";
	file << r << "\t"; // wrong
	file << z << "\t";
	file << n << "\n";
	//file << alpha_r(v) << "\t";
	//file << beta_r(v) << "\t";

}

void dendrite::writeTau(ofstream& file){
	file << tau_k() << "\t";
	file << tau_l(v) << "\t";
	file << tau_r(v) << "\t";
	file << tau_z(3.7152) << "\t";
	file << tau_n(v) << "\n";
}

void dendrite::writeIdensity(ofstream& file){
	file << iCal  << "\t";
	file << iCah  << "\t";
	file << iKca << "\t";
	file << iH  << "\t";
	file << iLeak << "\n";
}

dendrite::dendrite(double vInit): 
	caModel(13.3333333 * pow(10, -3), 0.000003, 3.7152 * pow(10, -3), 3 * pow(10,-6)){
	length = 0.1; // (cm) 0.1 cm = 1000 micron
	d = 0.0001; // (cm) 0.0001 cm = 1 micrion
	R_A = 3; // (ohm * cm)
	C_M = 0.000001; // (F) 0.000001 F / cm^2= 1 uF /cm^2
	v = vInit;
	lArea = PI * d * length; // (cm^2)
	R_a = (4 * length * R_A) / (PI * d * d); // (cm * ohm * cm / (cm * cm) = ohm) 
	//printf("PI %f\n", PI);
	//printf("d: %f\n", d);
	//printf("4 * l * R_A: %f\n", (4 * length * R_A));
	C_m = lArea * C_M; // (cm^2 * F/ cm^2 = F)

	//printf("dendrite R_a %f\n",R_a);
	// <species segmentGroup="dendrite_group" id="ca" ion="ca" concentrationModel="ca_conc" initialConcentration="3.7152 mM" initialExtConcentration="3.0mM"/>
	// caModel = caConcentrationModel(13.3333333, 0.000003, 3.7152, 3.0);
	k = inf_k(v);
	l = inf_l(v);
	r = inf_r(v);
	z = inf_z(getCaConc());
	n = inf_n(v);
	iCal = I_cal(v, k , l);
	iCah = I_cah(v, r);
	iKca = I_kca(v, z);
	iH = I_h(v);
	iLeak = I_leak(v, lArea);
}

void dendrite::writeInfK(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_inf_k.dat");
	file << "v\tinf_k\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_k(vTemp) << "\n";
	}
	file.close();

}

void dendrite::writeTauK(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_tau_k.dat");
	file << "v\ttau_k\n";
	int nSteps = (vEnd - vStart) / vStep;
	for(int i = 0; i < nSteps ; i++){
		file << vStart + i * vStep << "\t" << tau_k() << "\n";
	}
	file.close();

}

void dendrite::writeInfL(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_inf_l.dat");
	file << "v\tinf_l\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_l(vTemp) << "\n";
	}
	file.close();
}

void dendrite::writeTauL(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_tau_l.dat");
	file << "v\ttau_l\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_l(vTemp) << "\n";
	}
	file.close();
}

void dendrite::writeTauAndInfCal(double vStart, double vEnd, double vStep){
	writeTauK(vStart, vEnd, vStep);
	writeInfK(vStart, vEnd, vStep);
	writeTauL(vStart, vEnd, vStep);
	writeInfL(vStart, vEnd, vStep);
}

void dendrite::writeTauR(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_tau_r.dat");
	file << "v\ttau_r\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_r(vTemp) << "\n";
	}
	file.close();
}

void dendrite::writeInfR(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_inf_r.dat");
	file << "v\tinf_r\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << inf_r(vTemp) << "\n";
	}
	file.close();
}

void dendrite::writeTauAndInfCah(double vStart, double vEnd, double vStep){
	writeTauR(vStart, vEnd, vStep);
	writeInfR(vStart, vEnd, vStep);
}

void dendrite::writeTauN(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_tau_n.dat");
	file << "v\ttau_n\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_n(vTemp) << "\n";
	}
	file.close();
}

void dendrite::writeInfN(double vStart, double vEnd, double vStep){
	ofstream file;
	file.setf(ios::scientific, ios::floatfield);
	file.precision(6);
	file.open("channel_summary/dend_inf_n.dat");
	file << "v\tinf_n\n";
	int nSteps = (vEnd - vStart) / vStep;
	double vTemp;
	for(int i = 0; i < nSteps ; i++){
		vTemp = vStart + i * vStep;
		file << vTemp << "\t" << tau_n(vTemp) << "\n";
	}
	file.close();
}

void dendrite::writeTauAndInfH(double vStart, double vEnd, double vStep){
	writeTauN(vStart, vEnd, vStep);
	writeInfN(vStart, vEnd, vStep);
}

void dendrite::writeAllTauAndInf(double vStart, double vEnd, double vStep){
	writeTauAndInfCal(vStart, vEnd, vStep);
	writeTauAndInfCah(vStart, vEnd, vStep);
	writeTauAndInfH(vStart, vEnd, vStep);
}


 