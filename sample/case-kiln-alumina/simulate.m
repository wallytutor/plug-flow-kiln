%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup and run a rotary kiln simulation.
%
% Author: Walter Dal'Maz Silva
% Date  : April 11th 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clear functions; close all force;

warning("off");
addpath("../../src/octave");
addpath("../../dist/casadi-3.6.3-windows64-octave7.3.0");

import casadi.*;
import kiln.*;

name = "simulate.json";
data = load_json(strcat("cases/", name));
opts = struct("ipopt", data.ipopt);

% Add path for finind ProjectData.m.
addpath(data.project_data);

% XXX: check of material properties.
% T = linspace(298.0, 2327.0);
% bed = BedMixture(data=data.materials.bed_data);
% cp = bed.specific_heat_solid_mole(T);

model = RotaryKilnModel(data, devel=false);
model.solve_system(opts);
model.plot_results();
model.dump_results(strcat("results/", name), debug=true, extras=true);

% hmg = model.gas.enthalpy_mass(model.T_g, model.Y_g);
% cpg = model.gas.specific_heat_mass(model.T_g, model.Y_g);
% qdot_i = model.mdot_g(1) * cpg(1) * model.T_g(1);
% qdot_o = model.mdot_g(1) * cpg(end) * model.T_g(end);

% qdot_b = sum(model.bal_b);
% qdot_g = sum(model.bal_g);
% qdot_c = sum(model.bal_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%