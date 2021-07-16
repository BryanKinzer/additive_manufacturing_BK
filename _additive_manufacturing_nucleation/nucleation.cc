#include <iostream>
using namespace std;
double cum_prob {0.0};
// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const
{
	//Supersaturation factor
    double ssf;
    
    if (p[0] > 170){
	//was ssf = 1.23
	ssf = 1.23;
    }else{
	ssf = 1.23;
    }
  
	// Calculate the nucleation rate
	double J=k1*exp(-k2/(std::max(ssf,1.0e-6)))*exp(-tau/(this->currentTime));
	double retProb=1.0-exp(-J*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*dV);
	cum_prob = cum_prob + retProb;
	//cout << "p[0]: " << p[0] << endl;
    //p[0] > 9.5e-5 allows for heterogeneous nucleation close to the mold wall 9.5e-5/1e-4 = 0.95 of domain
    if (p[0] > 9.5e-5){
	//right now retProb is hard coded. Can make this smarter in the future, using classic nucleation theory
	retProb = 0.00002;
    }else{
	//retProb = 0.0;
    }
    return retProb;
}
