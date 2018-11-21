#ifndef INCLUDED_N_AL_PN_HPP
#define INCLUDED_N_AL_PN_HPP

#include "insilico/core/engine.hpp"
#include "insilico/core/injector.hpp"

#include "I_NA_AL_PN.hpp"
#include "I_K_AL_PN.hpp"
#include "I_A_AL_PN.hpp"
#include "I_LEAK_AL_PN.hpp"
#include "I_KLEAK_AL_PN.hpp"
#include "S_nACH.hpp"
#include "S_2GABAA.hpp"

#include <random>

namespace insilico {

class PN : public Neuron {
 public:
  void ode_set(state_type &variables, state_type &dxdt,const double t, unsigned index) {

    std::vector<unsigned> g1_indices;
    std::vector<double> gsyn_values;
    std::vector<double> esyn_values;

    std::vector<double> I_IExt;

    unsigned v_index = engine::neuron_index(index, "v");
    double v = variables[v_index];

    // Adding the external injected current (I_IExt)
    I_IExt = injector::external_current(index, t);

    double current = 0;
    if(!I_IExt.empty()) {
      current = I_IExt[0];
    }

    // incoming synaptic currents
    double I_syn = 0;
    g1_indices = engine::get_pre_neuron_indices(index, "g1");
    gsyn_values = engine::get_pre_neuron_values(index, "gsyn");
    esyn_values = engine::get_pre_neuron_values(index, "esyn");
    
    for(unsigned iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_syn += variables[g1_indices[iterator]] * gsyn_values[iterator] * (v - esyn_values[iterator]);
    }

    //std::cout << "g1_indices_size" << index <<" : " << g1_indices.size()<< std::endl;
    
    engine::neuron_value(index, "I_syn", I_syn);
    
    // note the spike
    double last_spiked = engine::neuron_value(index, "last_spike");
    double spike_duration = 5;
    double v_thresh = 0;

    // associated delay for next spikes
    if((v > v_thresh) && (t - last_spiked) > spike_duration){
      engine::neuron_value(index, "last_spike", t);
    }

    // ODE set   
    I_NA::current(variables, dxdt, t, index);
    I_K::current(variables, dxdt, t, index);
    I_A::current(variables, dxdt, t, index);
    I_LEAK::current(variables, dxdt, t, index);
    I_KLEAK::current(variables, dxdt, t, index);

    double I_Na = engine::neuron_value(index, "I_NA");
    double I_K = engine::neuron_value(index, "I_K");
    double I_A = engine::neuron_value(index, "I_A");
    double I_l = engine::neuron_value(index, "I_LEAK");
    double I_kl = engine::neuron_value(index, "I_KLEAK");
    double I_ext = engine::neuron_value(index, "I_ext");   	// The external current provided to all neurons
    
    dxdt[v_index] = I_ext + current - I_Na - I_K - I_A - I_l - I_kl - I_syn;

  } // function ode_set
}; // class N_AL_PN_HH

} // insilico

#endif

//#include "insilico/core/injector.hpp"
//#include "I_A_ANTENNAL_LOBE_PN_HH.hpp"
//#include "S_nACH.hpp"
//#include "I_KLEAK_ANTENNAL_LOBE_PN_HH.hpp"
//#include "insilico/neuron/helper/spike_list.hpp"
    
    /* Check (temporary)
    std::cout << "Neuron Index : " << index << std::endl;
    std::cout << "V : " << v << std::endl;
    std::cout << "I_Na : " << I_Na << std::endl;
    std::cout << "I_K : " << I_K << std::endl;
    std::cout << "I_A : " << I_A << std::endl;
    std::cout << "I_Leak : " << I_l << std::endl;
    std::cout << "I_KLeak : " << I_kl << std::endl;
    //std::cout << "I_Syn : " << I_Syn << std::endl;
    std::cout << "I_ext : " << I_ext << std::endl;
    std::cout << "Injected current : " << current << std::endl;
    */

    //engine::neuron_value(index, "ik", I_K);
    //double I_A = engine::neuron_value(index, "I_A_ANTENNAL_LOBE_PN_HH");
    //double I_kl = engine::neuron_value(index, "I_KLEAK_ANTENNAL_LOBE_PN_HH");  
    //ODE set
    //I_A_ANTENNAL_LOBE_PN_HH_::current(variables, dxdt, t, index);
    //I_KLEAK_ANTENNAL_LOBE_PN_HH_::current(variables, dxdt, t, index);

    /*
    double gsyn = 0.5; 
    
    double cm = 1.43 * (pow(10, -4));

    // incoming synaptic currents
    double I_Syn = 0;
    g1_indices = engine::get_pre_neuron_indices(index, "g1");
    esyn_values = engine::get_pre_neuron_values(index, "esyn");

    for(unsigned iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_Syn += variables[g1_indices[iterator]] * gsyn * (v - esyn_values[iterator]);
    }
    

    //Adding the external injected current (I_IExt)
    I_IExt = injector::external_current(index, t);
    double current = 0;
    if(!I_IExt.empty()) {
      current = I_IExt[0];
    }
    */
    //std::vector<unsigned> g1_indices;
    //std::vector<double> esyn_values;
    //std::vector<double> I_IExt;
    //double v = variables[v_index];

/*
    double last_spiked = engine::neuron_value(index, "last_spike");
    double spike_duration = engine::neuron_value(index, "spike_duration");
    double v_thresh = engine::neuron_value(index, "thresh");
*/
