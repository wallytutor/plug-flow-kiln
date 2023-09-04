%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup and run a rotary kiln simulation.
%
% Author: Walter Dal'Maz Silva
% Date  : April 11th 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clear functions; close all force; clc;

warning("off");
addpath("../../src/octave");
addpath("../../dist/casadi-3.6.3-windows64-octave7.3.0");

import casadi.*;
import kiln.*;

case_name = "cases/simulate-barr.json";
save_as = "results/results-barr.json";
[data, model] = run_case(case_name, save_as, plot_results=false, devel=false);

model.plot_heat_fluxes();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%