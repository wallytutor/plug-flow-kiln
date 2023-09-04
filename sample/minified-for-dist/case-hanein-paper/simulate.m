%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup and run a rotary kiln simulation.
%
% Author: Walter Dal'Maz Silva
% Date  : April 11th 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clear functions; close all force; clc;

warning("off");
addpath("../../environment/distribution-0.1.0-octave-8.2.0-w64");
addpath("../../environment/casadi-3.6.3-windows64-octave7.3.0");

import casadi.*;
import kiln.*;

case_name = "simulate.json";
save_as = "results.json";
[data, model] = run_case(case_name, save_as, plot_results=true, devel=false);

model.plot_heat_fluxes();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%