#ifndef HELPER_H
#define HELPER_H
 
const double PI = 3.141592653589793;

double expTime(double v, double tau, double midpoint, double scale);

double cal_tau(double v);

double h_tau(double v);

double kdr_tau(double v);

double hhExpRate(double v, double midpoint, double rate, double scale);

double hhExpLinearRate(double v, double midpoint, double rate, double scale);

double hhSigmoidRate(double v, double midpoint, double rate, double scale);

double HHSigmoidVariable(double v, double midpoint, double rate, double scale);

double I_leak(double v, double lA);

double I_pulse(double t, double delay, double duration, double amplitude);

double fwd_euler(double x, double dxdt, double dt);
 
// This is the end of the header guard
#endif