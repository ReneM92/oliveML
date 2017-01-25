#ifndef CAMODEL_H
#define CAMODEL_H

class caConcentrationModel{
	private:
		double decayConstant, rho, concentration, extConcentration; // restingConc = 0;
	public:
		caConcentrationModel(double dc, double r, double init_conc, double init_extConc);
		double getConcentration();
		void update(double iCa, double lA, double dt);
};

#endif