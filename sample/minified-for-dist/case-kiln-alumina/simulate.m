%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup and run a rotary kiln simulation.
%
% Author: Walter Dal'Maz Silva
% Date  : April 11th 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clear functions; close all force;

warning("off");
addpath("../../environment/kilnerys-0.1.0-octave-8.2.0-w64");
addpath("../../environment/casadi-3.6.3-windows64-octave7.3.0");

import casadi.*;
import kiln.*;

name = "simulate.json";
data = load_json(name);
opts = struct("ipopt", data.ipopt);

model = RotaryKilnModel(data, devel=false);
model.solve_system(opts);
model.plot_results();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
