#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <fstream>
#include "dendrite.hpp"
#include "axon.hpp"
#include "soma.hpp"

using namespace std;

class cell{
	private:
		dendrite dendrite;
		axon axon;
		soma soma;
		double dt;
	public:
		cell(double dt, double vInit);

		void update(double iIn);
		void writeToFile(int i, ofstream& file);
		void writeAllToFile(int i, ofstream& file, ofstream& dendrite_file, ofstream& axon_file, ofstream& soma_file);
		void writeTaus(ofstream& dendrite_file, ofstream& axon_file, ofstream& soma_file);
		void writeIdensity(ofstream& dendrite_file, ofstream& axon_file, ofstream& soma_file);
		void writeAllTauAndInf(double vStart, double vEnd, double vStep);
		/*
		dendrite getDendrite();
		axon getAxon();
		soma getSoma();
		*/
};

#endif