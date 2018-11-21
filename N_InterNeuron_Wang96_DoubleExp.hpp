#ifndef INCLUDED_N_INTERNEURON_WANG96_HPP
#define INCLUDED_N_INTERNEURON_WANG96_HPP

#include <cmath>
#include <iostream>
#include "insilico/core/engine.hpp"
#include "currents/I_Na_InterNeuron_Wang96.hpp"
#include "currents/I_Leak_InterNeuron_Wang96.hpp"
#include "currents/I_K_InterNeuron_Wang96.hpp"
#include "currents/I_PeriodicPulse.hpp"
#include "currents/I_OscillatoryDrive.hpp"

#include "../../noise.hpp"

namespace insilico {

class N_InterNeuron_Wang96 : public Neuron {
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t, unsigned index) {

    //std::cout << " index is I neuron " << index << std::endl ;

    int v_index = engine::neuron_index(index, "v");

    double v = variables[v_index];
    // note the spike
    
    double I_syn;// = engine::neuron_value(index,"I_syn");
    
    I_syn = 0;
    std::vector<unsigned int> s_indices = engine::get_pre_neuron_indices(index, "g1");
    std::vector<double> esyn_values = engine::get_pre_neuron_values(index, "esyn");
    std::vector<double> gsyn_values = engine::get_pre_neuron_values(index, "gsyn");

    for(std::vector<int>::size_type iterator = 0; iterator < s_indices.size() ; ++iterator) {
      I_syn = I_syn + gsyn_values[iterator]*variables[s_indices[iterator]] * (v - esyn_values[iterator]);
    }

    //std::vector<int>::size_type iteratr   = s_indices.size() - 1;
    //I_syn = I_syn + gsyn_values[iteratr]*variables[s_indices[iteratr]] * ( (v - esyn_values[iteratr]) / (1.0 + 0.15*(exp(-0.08*(v - esyn_values[iteratr])))));

    double spike_threshold = engine::neuron_value(index, "spike_threshold");
    double is_spiked = engine::neuron_value(index, "is_spiked");
    double last_spiked = engine::neuron_value(index, "last_spiked");
    double spike_duration = engine::neuron_value(index, "spike_duration");
    
    if (v >= spike_threshold) {
      
      if (is_spiked == 0){
        last_spiked = t;
        is_spiked = 1.0;
        engine::neuron_value(index, "is_spiked", is_spiked);
        
      }
      
      if (is_spiked == 1) {
        if ((t - last_spiked) >= spike_duration){
          last_spiked = t;
        }
        
      }
    }

    engine::neuron_value(index, "last_spiked", last_spiked);
    engine::neuron_value(index,"I_syn",I_syn);

    // ODE set
    I_Na_InterNeuron_Wang96::current(variables, dxdt, t, index);
    I_K_InterNeuron_Wang96::current(variables, dxdt, t, index);
    I_Leak_InterNeuron_Wang96::current(variables, dxdt, t, index);
    I_PeriodicPulse::current(variables, dxdt, t, index);
    I_OscillatoryDrive::current(variables,dxdt,t,index);
    //I_PulseForOsc::current(variables,dxdt,t,index);

    double I_Na = engine::neuron_value(index, "I_Na_InterNeuron_Wang96");
    double I_K = engine::neuron_value(index, "I_K_InterNeuron_Wang96");
    double I_Leak = engine::neuron_value(index, "I_Leak_InterNeuron_Wang96");
    double I_Ext = engine::neuron_value(index, "I_PeriodicPulse");
    double I_DC = engine::neuron_value(index, "I_DC");
    double I_OscillatoryDrive = engine::neuron_value(index, "I_OscillatoryDrive");
    //double I_PulseForOsc = 

    double current = 0;
    current = noise::inject(index);
    double ampNoise = engine::neuron_value(index, "ampNoise");
    double I_Noise = ampNoise*current*(v+0.0);

    engine::neuron_value(index,"I_Noise",I_Noise);
    
    dxdt[v_index] = I_Ext + I_DC - I_Na - I_K - I_Leak - I_syn - I_Noise - I_OscillatoryDrive;
  } // function ode_set
}; // class N_InterNeuron_Wang96

} // insilico

#endif










