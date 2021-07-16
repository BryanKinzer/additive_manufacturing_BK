// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){

    for (unsigned int var_index=0; var_index<6; var_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(var_index));

        set_variable_name				(var_index,var_name);
    	set_variable_type				(var_index,SCALAR);
    	set_variable_equation_type		(var_index,EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(var_index, "n0, n1, n2, n3, n4, n5, z");
        set_dependencies_gradient_term_RHS(var_index, "grad(n0), grad(n1), grad(n2), grad(n3), grad(n4), grad(n5)");
	set_allowed_to_nucleate			(var_index, true);
	set_need_value_nucleation		(var_index, true);
    }

		set_variable_name				(6,"z");
		set_variable_type				(6,SCALAR);
		set_variable_equation_type		(6,EXPLICIT_TIME_DEPENDENT);

    	set_dependencies_value_term_RHS(6, "n0, n1, n2, n3, n4, n5, z");
    	set_dependencies_gradient_term_RHS(6, "grad(z)");
	set_allowed_to_nucleate			(6, false);
        set_need_value_nucleation		(6, false);

		set_variable_name				(7,"T");
		set_variable_type				(7,SCALAR);
		set_variable_equation_type		(7,EXPLICIT_TIME_DEPENDENT);

    	set_dependencies_value_term_RHS(7, "T");
    	set_dependencies_gradient_term_RHS(7, "grad(T)");
	set_allowed_to_nucleate			(7, false);
        set_need_value_nucleation		(7, false);

	// set_variable_name				(8,"phi");
	// set_variable_type				(8,SCALAR);
	// set_variable_equation_type		(8,AUXILIARY);

    // set_dependencies_value_term_RHS(8, "T");
    // set_dependencies_gradient_term_RHS(8, "");

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

dealii::VectorizedArray<double> fnV = constV(0.0);
scalarvalueType ni, nj;
scalargradType nix;


scalarvalueType z = variable_list.get_scalar_value(6);
scalargradType zx = variable_list.get_scalar_gradient(6);

scalarvalueType T = variable_list.get_scalar_value(7);
scalargradType Tx = variable_list.get_scalar_gradient(7);

// In this application, create temporary variables for the residual terms. We cannot
// call 'set_scalar_value_residual_term' and 'set_scalar_gradient_residual_term'in the
// for loop below because those functions write over the scalar value and scalar gradient
// internal variables in 'variable_list' (for performance reasons). Therefore, we wait
// to set the residual terms until all the residuals have been calculated.

std::vector<scalarvalueType> value_terms;
value_terms.resize(userInputs.number_of_variables);
std::vector<scalargradType> gradient_terms;
gradient_terms.resize(userInputs.number_of_variables);

for (unsigned int i=0; i<userInputs.number_of_variables-2; i++){

	ni = variable_list.get_scalar_value(i);
	nix = variable_list.get_scalar_gradient(i);
	fnV = - ni + ni*ni*ni; +(2.0*ni*(1.0-z)*(1.0-z));
	for (unsigned int j=0; j<userInputs.number_of_variables-2; j++){
		if (i != j){
			nj = variable_list.get_scalar_value(j);
			fnV += constV(2.0*gamma) * ni * nj*nj;
		}
	}
	value_terms[i] = ni-constV(userInputs.dtValue*Lg*mg)*fnV;
	gradient_terms[i] = constV(-userInputs.dtValue*Kg*Lg)*nix;
}

// --- Submitting the terms for the governing equations ---

// ---------------------------------------------------------------------
// Nucleation section to seed nuclei and freeze the Allen-Cahn mobility
// ---------------------------------------------------------------------
//probably a way to reduce all source_term doubles into single variable
dealii::VectorizedArray<double> source_term0 = constV(0.0);
dealii::VectorizedArray<double> source_term1 = constV(0.0);
dealii::VectorizedArray<double> source_term2 = constV(0.0);
dealii::VectorizedArray<double> source_term3 = constV(0.0);
dealii::VectorizedArray<double> source_term4 = constV(0.0);
dealii::VectorizedArray<double> source_term5 = constV(0.0);
dealii::VectorizedArray<double> lambda = constV(1.0);
seedNucleus(q_point_loc,source_term0,source_term1,source_term2,source_term3,source_term4,source_term5,lambda);
// ---------------------------------------------------------------------


for (unsigned int i=0; i<userInputs.number_of_variables-2; i++){
	if (i == 0) {
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]+source_term0);
	}else if (i == 1){
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]+source_term1);
	}else if (i == 2){
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]+source_term2);
	}else if (i == 3){
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]+source_term3);
	}else if (i == 4){
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]+source_term4);
	}else if (i == 5){
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]+source_term5);
	}else{
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]);
	}
	variable_list.set_scalar_gradient_term_RHS(i,gradient_terms[i]);
}





//scalarvalueType phi=constV(0.5*(1.0-std::tanh(theta*((/*Temperature*/Tliquidus)-1.0))));
scalarvalueType phi = constV(0.0);

for (unsigned i=0; i< T.n_array_elements;i++){
	phi[i] = (0.5)*((1.0)-std::tanh(theta*((T[i]/(Tliquidus)-1.0))));
}

// phi=std::min(phi[0],(1.0));
// phi=std::max(phi[0],(0.0));

//value terms
dealii::VectorizedArray<double> fnz_p = constV(0.0);
dealii::VectorizedArray<double> fnz_g = constV(0.0);
dealii::VectorizedArray<double> fnz_g_sum = constV(0.0);

fnz_p=2.0*z*(1.0-phi)-2.0*phi*(1.0-z);
fnz_g=(2.0*(1.0-z));
for (unsigned int i=0; i<userInputs.number_of_variables-2; i++){
	ni = variable_list.get_scalar_value(i);
    fnz_g_sum+=ni*ni;
    }


scalarvalueType value_z = z-(constV(userInputs.dtValue*Lp*mp)*fnz_p)+(constV(userInputs.dtValue*Lp*mg)*fnz_g*fnz_g_sum);
scalargradType grad_z = constV(-userInputs.dtValue*Lp*Kp)*zx;

variable_list.set_scalar_value_term_RHS(6,value_z+source_term0 + source_term1 + source_term2 + source_term3 + source_term4 + source_term5);
variable_list.set_scalar_gradient_term_RHS(6,grad_z);

// std::tanh(z[0]);

////////Dummy stuff for now////////



scalarvalueType eq_T=T;
scalargradType eqx_T=(constV(-D*userInputs.dtValue)*Tx);

variable_list.set_scalar_value_term_RHS(7,eq_T);
variable_list.set_scalar_gradient_term_RHS(7,eqx_T);

////////Dummy stuff for now////////



// dealii::VectorizedArray<double> test = constV(0.0);

// test = (0.5*(1.0-std::tanh(theta*((T[0]/Tliquidus)-1.0))));
// constV(0.5*(1.0-std::tanh(theta*((/*Temperature*/Tliquidus)-1.0))));

}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim,int degree>
void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
	dealii::VectorizedArray<double> & source_term0,
	dealii::VectorizedArray<double> & source_term1,
	dealii::VectorizedArray<double> & source_term2,
	dealii::VectorizedArray<double> & source_term3,
	dealii::VectorizedArray<double> & source_term4,
	dealii::VectorizedArray<double> & source_term5,
	dealii::VectorizedArray<double> & lambda) const {

		for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){
			if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){
				// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
				dealii::VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex), q_point_loc, thisNucleus->orderParameterIndex);

				for (unsigned i=0; i<lambda.n_array_elements;i++){
					if (weighted_dist[i] <= 1.0){
						lambda[i] = 0.0;

						// Seed a nucleus if it was added to the list of nuclei this time step
						if (thisNucleus->seedingTimestep == this->currentIncrement){
							// Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term
							dealii::Point<dim,double> q_point_loc_element;
							for (unsigned int j=0; j<dim; j++){
								q_point_loc_element(j) = q_point_loc(j)[i];
							}
							double r = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_semiaxes(thisNucleus->orderParameterIndex), q_point_loc_element, thisNucleus->orderParameterIndex);

							double avg_semiaxis = 0.0;
							for (unsigned int j=0; j<dim; j++){
								avg_semiaxis += thisNucleus->semiaxes[j];
							}
							avg_semiaxis /= dim;

							if (thisNucleus->orderParameterIndex == 0){
								//was interface_coeff instead of Lp
								source_term0[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 1){
								source_term1[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 2){
								source_term2[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 3){
								source_term3[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 4){
								source_term4[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 5){
								source_term5[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else {
								//source_terms[2][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							//source_term[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));

						}
					}
				}
			}
		}
	}



// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
