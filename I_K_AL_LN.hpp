//CREATED BY SHREYA ON 23 MARCH 2018

#ifndef INCLUDED_I_K_AL_LN_HPP
#define INCLUDED_I_K_AL_LN_HPP

#include "insilico/core/engine.hpp"

namespace insilico {

class I_K_LN_ {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, int index) {
    double gk = 10.0, ek = -95.0, vtr_k = -50, Cels = 22.0;

    //variables
    int v_index = engine::neuron_index(index, "v");
    int nk_index = engine::neuron_index(index, "nk");

    double v = variables[v_index];
    double nk = variables[nk_index];

    double v2 = v - vtr_k; //This term scales the voltage value

    double alpha_nk =  0.02*(15.0 - v2)/(exp((15.0 - v2)/5.0) - 1.0);
    double beta_nk = 0.5*exp((10.0-v2) / 40.0);
    double phi = pow(3,((Cels-36)/10));

    double tau_nk = 1/(alpha_nk + beta_nk) / phi;
    double nk_inf = alpha_nk/(alpha_nk+ beta_nk);

    // ODE set
    dxdt[nk_index] = -(nk - nk_inf)/tau_nk;

    // Current
    engine::neuron_value(index, "I_K_LN", (gk * pow(nk,4) * (v - ek)));

    } // function current
}; // class I_K_LN_

} // insilico

#endif

/*
G_Kv = 10;///////////////////////                                                                                                    
    Vtr = -50;
    S2 = 0.02;
    v2 = v - Vtr;
    Phi = pow(3,((Cels-36)/10));

    Alpha3 = 0.032*(15 - v2)/(exp((15 - v2)/5) - 1);
    Beta3 = 0.5*exp((10 - v2)/40);
    n0 = Alpha3/(Alpha3 + Beta3);     }
  void calc(double n, double &fn, double v, double x);
};

double IKv_LN::E_K = -95, IKv_LN::Cels = 22;

void IKv_LN::calc(double n, double &fn, double v, double x){

  v2 = v - Vtr;
  iKv = G_Kv* n*n*n*n*(v - E_K);
  Alpha3 = S2*(15 - v2)/(exp((15 - v2)/5) - 1);
  Beta3 = 0.5*exp((10 - v2)/40);
  tau_n = 1/(Alpha3 + Beta3) / Phi;
  n_inf = Alpha3/(Alpha3 + Beta3);

  fn  = -(n - n_inf)/tau_n;


*/
