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

% data = load_json("simulate-source-ebu.json");
% data = load_json("simulate-source-mak.json");
data = load_json("simulate-equilibrated.json");
opts = struct("ipopt", data.ipopt);

model = RotaryKilnModel(data, devel=true);
model.solve_system(opts);
model.plot_results();
model.dump_results("results.json", debug=true, extras=true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%