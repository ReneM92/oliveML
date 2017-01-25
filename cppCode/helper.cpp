#include <iostream>
#include <math.h>
#include "helper.hpp"
using namespace std;

double expTime(double v, double tau, double midpoint, double scale){
	double t;
	t = tau * exp((v - midpoint)/scale);
	return t;
}

double cal_tau(double v){
	double t;
	t = ((20 * exp((v + 160) / 30) / (1 + exp((v + 84) / 7.3))) + 35);
	return t;
}

double h_tau(double v){ // TIME_SCALE = 1, Tau_Ih in h.channel.nml
	double t;
	t = 1 /(exp(-0.086 * v - 14.6) + exp(0.070 * v - 1.87));
	return t;
}

double kdr_tau(double v){
	double t;
	t = (5 + (47 * exp(-1 * (50 - v) / 900)));
	return t;
}

double hhExpRate(double v, double midpoint, double rate, double scale){
	double r;
	r = rate * exp((v-midpoint)/scale);
	return r;
}

double hhExpLinearRate(double v, double midpoint, double rate, double scale){
	double x, r;
	x = (v - midpoint) / scale;
	if(x == 0){
		r = rate;
	}else{
		r = rate * x / (1 - exp(0 - x));
	}
	return r;
}

double hhSigmoidRate(double v, double midpoint, double rate, double scale){
	double r;
	r = rate / (1 + exp(0 - (v - midpoint)/scale));
	return r;
}

double HHSigmoidVariable(double v, double midpoint, double rate, double scale){
	double x;
	x = rate / (1+ exp(0 - (v - midpoint) / scale));
	return x;
}

double I_leak(double v, double lA){
	double lArea = lA;
	double condDensity = 0.015 * pow(10,-3);
	double erev = -10 * pow(10, -3);
	return condDensity * (erev -v) * lArea;
}

double I_pulse(double t, double delay, double duration, double amplitude){ // 
	return (t>delay) * amplitude - (t>(delay + duration)) * amplitude; 
}

double fwd_euler(double x, double dxdt, double dt){
	return x + dxdt * dt;
}