#include "caModel.hpp" 
#include "helper.hpp"

double caConcentrationModel::getConcentration(){
	return concentration;
}

void caConcentrationModel::update(double iCa, double lA, double dt){
	double dconcdt = (iCa/lA) * rho - (concentration / decayConstant);
	concentration = fwd_euler(concentration, dconcdt, dt);
	if(concentration < 0){
		concentration = 0;
	}	
}

caConcentrationModel::caConcentrationModel(double dc, double r, double init_conc, double init_extConc){
	decayConstant = dc;
	rho = r;
	concentration = init_conc;
	extConcentration = init_extConc;
}