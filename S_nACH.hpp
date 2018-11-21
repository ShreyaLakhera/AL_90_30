/*
 synapse/S_nACH

 Copyright (C) 2014 Pranav Kulkarni, Collins Assisi Lab, IISER, Pune <pranavcode@gmail.com>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//EDITED BY SHREYA 

#ifndef INCLUDED_S_nACH_HPP
#define INCLUDED_S_nACH_HPP

#include "insilico/core/engine.hpp"

namespace insilico {
class S_nACH: public Synapse {
 public:
  void ode_set(state_type& variables, state_type& dxdt, const double t, unsigned index) {
    
    unsigned g1_index = engine::synapse_index(index, "g1");
    double g1 = variables[g1_index];

    //unsigned post_neuron_index = engine::synapse_value(index, "post");
    unsigned pre_neuron_index = engine::synapse_value(index, "pre");
    double t_0 = engine::neuron_value(pre_neuron_index, "last_spike");

    // Constants
    double alpha = 10.0;
    double beta = 0.2;
    double A = 0.5;
    double t_max = 0.3;

    // Variables (T)
    double h1;
    if ((t_0 + t_max - t) < 0) {
    	h1 = 0;
    }
    else {
	h1 = 1;
    }

    double h2;
    if ((t - t_0) < 0) {
    	h2 = 0;
    }
    else {
    	h2 = 1;
    }

    double T = A * h1 * h2;  

    // ODE set
    dxdt[g1_index] = (alpha * T * (1.0 - g1) - beta * g1); 
    engine::synapse_value(index, "g1", g1);

  } // function ode_set
    
}; // class S_nACH

} // insilico

#endif

    //unsigned pre_neuron_index = engine::synapse_value(index, "pre"); //obtain pre-neuron index from ssets.
    //double v_index = engine::neuron_index(pre_neuron_index,"v"); //obtain the voltage of the pre-neuron
    //double v = variables[v_index];
    //double vthresh = -20.0;
    //double sig  = 1.5;
