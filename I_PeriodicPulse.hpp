/*
 current/I_PeriodicPulse.hpp - External current to be injected into the neuron

 Copyright (C) 2016 Collins Assisi Lab, IISER, Pune

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

/*
 Brief:
A pulse like current to be injected into the neuron. The baseline, maximum amplitude, start, end time
and the rise and fall time can be independently specified in the neuron.isf file for the respective
neuron.
*/

#ifndef INCLUDED_I_PERIODICPULSE_HPP
#define INCLUDED_I_PERIODICPULSE_HPP

#include "insilico/core/engine.hpp"

namespace insilico {

class I_PeriodicPulse {
 public:
  static void current(state_type &variables, state_type &dxdt, const double t, int index) {

    double PulseDuration = engine::neuron_value(index, "PulseDuration");      //duration of each cycle
    double PulseStart = engine::neuron_value(index, "PulseStart");            //Start of the pulse within the cycle
    double PulseEnd = engine::neuron_value(index, "PulseEnd");                //End of Pulse
    double tau_rise = engine::neuron_value(index, "tau_rise");                //Exponential rise time
    double tau_fall = engine::neuron_value(index, "tau_fall");                //Exponential fall time
    double PulseMax = engine::neuron_value(index, "PulseMax");                //Pulse Saturation
    double PulseMin = engine::neuron_value(index, "PulseMin");                //Pulse baseline

    double T = fmod(t,PulseDuration); //careful- check on the error of fmod operator!! 
    
    double P = PulseMin;
    if (T < PulseStart){
      P = PulseMin;
    }     
    else if ((T > PulseStart) && (T < PulseEnd)){
      P = PulseMax + (PulseMin - PulseMax)*exp(-(T - PulseStart)/tau_rise);
    }
    else if (T >= PulseEnd){
      P = PulseMin + (PulseMax - PulseMin)*exp(-(T - PulseEnd)/tau_fall);
    }

    // Current
    engine::neuron_value(index, "I_PeriodicPulse", P);

  } // function current
}; // class I_PeriodicPulse

} // insilico

#endif
