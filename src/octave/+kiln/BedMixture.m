classdef BedMixture < Thermodata
    methods (Access = public)
        function [obj] = BedMixture(data={})
            if isempty(data)
                data = obj.load_default_data("mixture_bed.json");
            endif

            tc = data.step_temperature / 1000;
            a_lo = data.shomate_poly_lo;
            a_hi = data.shomate_poly_hi;

            k_poly = data.thermal_conductivity_poly;
            k_min = data.thermal_conductivity_min;
            k_max = data.thermal_conductivity_max;

            c_lo = @(T) obj.shomate_specific_heat(T, a_lo);
            c_hi = @(T) obj.shomate_specific_heat(T, a_hi);
            h_lo = @(T) obj.shomate_enthalpy(T, a_lo);
            h_hi = @(T) obj.shomate_enthalpy(T, a_hi);

            cp_mole = @(T) obj.evaluate(T/1000, tc, c_lo, c_hi);
            cp_mass = @(T) cp_mole(T) ./ data.molecular_mass;

            h_mole = @(T) obj.evaluate(T/1000, tc, h_lo, h_hi);
            h_mass = @(T) h_mole(T) ./ data.molecular_mass;

            k = @(T) clip(polyval(k_poly, T), k_min, k_max);

            obj.repose_angle = data.repose_angle * pi / 180;
            obj.density_mass = data.density_mass;
            obj.solid_packing = data.solid_packing;
            obj.particle_diameter = data.particle_diameter;
            obj.molecular_mass = data.molecular_mass;
            obj.specific_heat_solid_mole = cp_mole;
            obj.enthalpy_mole = h_mole;
            obj.specific_heat_solid_mass = cp_mass;
            obj.enthalpy_mass = h_mass;
            obj.thermal_conductivity = k;
        endfunction

        function [result] = n_species(self)
            % Access to number of species in mechanism [-].
            result = 2;
        endfunction

        function [result] = species_names(self)
            % Access to list of species names.
            result = {"moisture", "silica"};
        endfunction

        function [val] = specific_heat_mass(self, T, Y)
            cp_w = self.WATER_SPECIFIC_HEAT;
            cp_b = self.specific_heat_solid_mass(T);
            val = Y(:, 1) .* cp_w + Y(:, 2) .* cp_b;
            val = cp_b;
        endfunction

        function [val] = evaporation_rate(self, T, q, Tb=368.0, k=2.0)
            % Hypotethical evaporation rate [kg/s].
            scale = q ./ self.WATER_LATENT_HEAT_EVAP;
            val = scale ./ (1.0 + exp(-k .* (T - Tb)));
        endfunction

        function [val] = mdot_evap(self, z, T, Y, q, mtol=1.0e-06)
            % Mass production rates of species [kg/s].
            h = (1 + sign(Y(:, 1) - mtol)) ./ 2;
            val = h .* self.evaporation_rate(T, q);
        endfunction
    
        function [val] = sdotk(self, z, T, Y, q)
            % Surface mass production rates of species [kg/s].
            mdot = self.mdot_evap(z, T, Y, q);
            val = mdot * self.SPECIES_COEFS;
        endfunction
    endmethods

    properties (Constant, Access = public, Static = true)
        % Species stoichiometric coefficients.
        SPECIES_COEFS = [-1.0, 0.0];
    endproperties

    properties (SetAccess = private, GetAccess = public)
        % Material repose angle [rad].
        repose_angle;

        % Material specific mass [kg/m³].
        density_mass;

        % Bed solid packing [-].
        solid_packing;

        % Particle diameter [m].
        particle_diameter;

        % Material mean molar mass [kg/mol].
        molecular_mass;

        % Molar specific heat function [J/(mol.K)].
        specific_heat_solid_mole;

        % Molar enthalpy function [J/mol].
        enthalpy_mole;

        % Mass specific heat function [J/(kg.K)].
        specific_heat_solid_mass;

        % Mass enthalpy function [J/kg].
        enthalpy_mass;

        % Thermal conductivity function [W/(m².K)].
        thermal_conductivity;
    end
end
