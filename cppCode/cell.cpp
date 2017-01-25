#include "cell.hpp"


cell::cell(double dtStep, double vInit):
	dendrite(vInit), soma(vInit), axon(vInit){
	dt = dtStep;

}

void cell::update(double iIn){
	//cout << "dtCell: " << dt << "\n";
	//printf("start update cell\n");
	double iInDendrite = iIn;
	//printf("iInDendrite: %.2f\n", iInDendrite);
	//printf("vDendrite: %.2f\n", dendrite.getV());
	//printf("vSoma%.2f\n", soma.getV());
	//printf("dendrite.Ra %f\n", dendrite.getRa());
	double iInSoma = (dendrite.getV() - soma.getV()) / dendrite.getRa();
	//printf("iInSoma: %.2f\n", iInSoma);
	double iInAxon = (soma.getV() - axon.getV()) / soma.getRa();
	//printf("iInAxon: %.2f\n", iInSoma);
	dendrite.update(iInDendrite, dt);
	soma.update(iInSoma, dt);
	axon.update(iInAxon, dt);
}

void cell::writeToFile(int i, ofstream& file){
	file << i << "\t";
	file << dendrite.getV() << "\t";
	file << axon.getV() << "\t";
	file << soma.getV() << "\t";
	file << dendrite.getCaConc() << "\n";
}

void cell::writeAllToFile(int i, ofstream& file, ofstream& dendrite_file, ofstream& axon_file, ofstream& soma_file){
	writeToFile(i, file);
	dendrite.writeToFile(dendrite_file);
	axon.writeToFile(axon_file);
	soma.writeToFile(soma_file);
}

void cell::writeTaus(ofstream& dendrite_file, ofstream& axon_file, ofstream& soma_file){
	dendrite.writeTau(dendrite_file);
	axon.writeTau(axon_file);
	soma.writeTau(soma_file);
}

void cell::writeIdensity(ofstream& dendrite_file, ofstream& axon_file, ofstream& soma_file){
	dendrite.writeIdensity(dendrite_file);
	axon.writeIdensity(axon_file);
	soma.writeIdensity(soma_file);
}

void cell::writeAllTauAndInf(double vStart, double vEnd, double vStep){
	dendrite.writeAllTauAndInf(vStart, vEnd, vStep);
	axon.writeAllTauAndInf(vStart, vEnd, vStep);
	soma.writeAllTauAndInf(vStart, vEnd, vStep);
}

/*
dendrite cell::getDendrite(){
	return dendrite;
}

axon cell::getAxon(){
	return axon;
}

soma cell::getSoma(){
	return soma;
}
*/