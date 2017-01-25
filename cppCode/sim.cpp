#include <iostream>
#include <fstream>
#include <math.h>
#include "cell.hpp"
#include "helper.hpp"
//#include <fenv.h>
//#include "exception.cpp"

using namespace std;


int main(){
	double simLength = 2000 * pow(10, -3); // s
	double simStep = 0.025 * pow(10, -3); // s
	/** variables I_clamp **/
	double delay = 500 * pow(10, -3);
	double duration = 200 * pow(10, -3);
	double amplitude = -0.3 * pow(10,-9);


	double vInit = -60 * pow(10,-3);
	cell cell1 = cell(simStep, vInit);

	
	//feenableexcept(FE_INVALID | FE_OVERFLOW);

	double iClamp = 0;

	ofstream file, dendrite_file, axon_file, soma_file, dendrite_tauf, axon_tauf, soma_tauf;
	ofstream axon_if, soma_if, dendrite_if;
	dendrite_if.setf(ios::scientific, ios::floatfield);
	axon_if.setf(ios::scientific, ios::floatfield);
	soma_if.setf(ios::scientific, ios::floatfield);
	dendrite_tauf.setf(ios::scientific, ios::floatfield);
	axon_tauf.setf(ios::scientific, ios::floatfield);
	soma_tauf.setf(ios::scientific, ios::floatfield);
	dendrite_file.setf(ios::scientific, ios::floatfield);
	axon_file.setf(ios::scientific, ios::floatfield);
	soma_file.setf(ios::scientific, ios::floatfield);
	file.setf(ios::scientific, ios::floatfield); // set fixed floating format
  	file.precision(6); // for fixed format, two decimal places
  	dendrite_file.precision(6);
  	axon_file.precision(6);
  	soma_file.precision(6);
  	dendrite_tauf.precision(6);
  	axon_tauf.precision(6);
  	soma_tauf.precision(6);
  	dendrite_if.precision(6);
  	axon_if.precision(6);
  	soma_if.precision(6);
  	dendrite_file.open("test_dendrite.dat");
  	axon_file.open("test_axon.dat");
  	soma_file.open("test_soma.dat");
	file.open("test.dat");
	dendrite_tauf.open("tau_dendrite.dat");
	axon_tauf.open("tau_axon.dat");
	soma_tauf.open("tau_soma.dat");
	dendrite_if.open("i_dendrite.dat");
	axon_if.open("i_axon.dat");
	soma_if.open("i_soma.dat");
	file << "i\tvDendrite\tvAxon\tvSoma\tcaDendrite\n";

	cout.setf(ios::scientific, ios::floatfield); // set fixed floating format
  	cout.precision(6); // for fixed format, two decimal places
	//cell1.writeToFile(0, file);
	cell1.writeAllToFile(0, file, dendrite_file, axon_file, soma_file);
	cell1.writeTaus(dendrite_tauf, axon_tauf, soma_tauf);
	cell1.writeIdensity(dendrite_if, axon_if, soma_if);
	/** file to write output values **/
	for(int i = 1; i * simStep <= simLength; i++){
		iClamp = I_pulse(i * simStep, delay, duration, amplitude);
		cell1.update(iClamp);
		cell1.writeToFile(i, file);
		cell1.writeAllToFile(i, file, dendrite_file, axon_file, soma_file);
		cell1.writeIdensity(dendrite_if, axon_if, soma_if);
	}

	file.close();
	dendrite_file.close();
	axon_file.close();
	soma_file.close();
	dendrite_tauf.close();
	axon_tauf.close();
	soma_tauf.close();
	dendrite_if.close();
	axon_if.close();
	soma_if.close();
	/*
	double vStart = -100 * pow(10, -3);
	double vEnd = 100 * pow(10, -3);
	double vStep = 1 * pow(10 , -3);

	cell1.writeAllTauAndInf(vStart, vEnd, vStep);
	*/

	return 0;
}