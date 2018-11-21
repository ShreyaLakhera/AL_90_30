// CREATED BY SHREYA SOMETIME IN FEBRUARY OR MARCH 2018

#ifndef INCLUDED_I_NA_AL_PN_HPP
#define INCLUDED_I_NA_AL_PN_HPP

#include "insilico/core/engine.hpp"

namespace insilico {

class I_NA {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, int index) {
    double gna = 100, ena = 50, vtr_na = -50, Cels = 22;
    
    int v_index = engine::neuron_index(index, "v");
    int m_index = engine::neuron_index(index, "m");
    int h_index = engine::neuron_index(index, "h");

    double v = variables[v_index];
    double m = variables[m_index];
    double h = variables[h_index];

    double v2 = v - vtr_na;
    double phi = pow(3,((Cels-36)/10));
    
    double alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
    double beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);

    double alpha2 = 0.128*exp((17 - v2)/18);
    double beta2 = 4/(exp((40 - v2)/5) + 1);

    double tau_m = (1/(alpha1 + beta1)) / phi;
    double m_inf = alpha1/(alpha1 + beta1);

    double tau_h = (1/(alpha2 + beta2)) / phi;
    double h_inf = alpha2/(alpha2 + beta2);

    dxdt[m_index] = -(m - m_inf)/tau_m;
    dxdt[h_index] = -(h - h_inf)/tau_h;

    double ina = gna * pow(m,3) * h * (v - ena);

    engine::neuron_value(index, "I_NA", ina);

    } // function current
}; // class I_NA

} // insilico

#endif

    /*
    double m_inf = 1/(1+exp(-(v+20)/6.5));
    double tau_m = 1.5;
    double h_inf = 1/(1+exp((v+25)/12));
    double tau_h = (0.3*exp((v-40)/13)) + (0.002*exp(-(v-60)/29));

    std::cout << "m_inf : " << m_inf << std::endl;
    std::cout << "h_inf : " << h_inf << std::endl;
    std::cout << "tau_m : " << tau_m << std::endl;
    std::cout << "tau_h : " << tau_h << std::endl;
    */


    
    /* ODE set
    dxdt[m_index]= (-1/tau_m) * (m - m_inf);
    dxdt[h_index]= (-1/tau_h) * (h - h_inf);
    */    
    // variables

//test equations:

    /*
    double alpha_m = (2.5 - 0.1 * v) / (exp(2.5-0.1 * v) - 1.0);
    double beta_m  = 4.0 * exp(-v / 18.0);
    double alpha_h = 0.07 * exp(-v / 20.0);
    double beta_h  = 1.0 / (exp(3 - 0.1 * v) + 1);

    dxdt[m_index] = (alpha_m * (1-m) - beta_m * m);
    dxdt[h_index] = (alpha_h * (1-h) - beta_h * h);
    */

//Equations from NCellNetwork...

    /*
    G_K = 10;///////////////////////                                                                                                                                                                                                                                               
    G_Na = 100;/////////////////////                                                                                                                                                                                                                                               
    Vtr = -50;
    VtrK = -50;
    S1 = 0.32;
    S2 = 0.02;
    v2 = v - Vtr;
    v2K = v - VtrK;
    Phi = pow(3,((Cels-36)/10));
    Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
    Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
    m0 = Alpha1/(Alpha1 + Beta1);

    Alpha2 = 0.128*exp((17 - v2)/18);
    Beta2 = 4/(exp((40 - v2)/5) + 1);
    h0 = Alpha2/(Alpha2 + Beta2);

    Alpha3 = 0.02*(15 - v2)/(exp((15 - v2)/5) - 1);
    Beta3 = 0.5*exp((10 - v2)/40);
    n0 = Alpha3/(Alpha3 + Beta3);     }
  void calc(double m, double h, double n, double &fm, double &fh, double &fn,
            double v, double x);
};

double INaK::E_K = -95, INaK::E_Na = 50, INaK::Cels = 22;

void INaK::calc(double m, double h, double n, double &fm, double &fh, double &fn,
                   double v, double x){
  v2 = v - Vtr;
  v2K = v - VtrK;
  iNa = G_Na*m*m*m*h*(v - E_Na);
  Alpha1 = S1*(13 - v2)/(exp((13 - v2)/4) - 1);
  Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
  tau_m = 1/(Alpha1 + Beta1) / Phi;
  m_inf = Alpha1/(Alpha1 + Beta1);

  Alpha2 = 0.128*exp((17 - v2)/18);
  Beta2 = 4/(exp((40 - v2)/5) + 1);
  tau_h = 1/(Alpha2 + Beta2) / Phi;
  h_inf = Alpha2/(Alpha2 + Beta2);

  fm = -(m - m_inf)/tau_m;
  fh = -(h - h_inf)/tau_h;

  iK = G_K* n*n*n*n*(v - E_K);
  Alpha3 = S2*(15 - v2K)/(exp((15 - v2K)/5) - 1);
  Beta3 = 0.5*exp((10 - v2K)/40);
  tau_n = 1/(Alpha3 + Beta3) / Phi;
  n_inf = Alpha3/(Alpha3 + Beta3);

  fn  = -(n - n_inf)/tau_n;
  */
