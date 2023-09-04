classdef GasMixture < Thermodata
    methods (Access = public)
        function obj = GasMixture(data={}, kin_name="EBU")
            % Constructs an instance of gas mixture.
            if isempty(data)
                data = obj.load_default_data("mixture_gas.json");
            endif

            obj.species = {};
            obj.species_names = {};
            obj.mw = [];

            for k=1:length(data.species)
                s = data.species(k);

                tc = s.step_temperature;
                a_lo = s.nasa7_poly_lo';
                a_hi = s.nasa7_poly_hi';

                c_lo = @(T) obj.poly_nasa7_specific_heat(T, a_lo);
                c_hi = @(T) obj.poly_nasa7_specific_heat(T, a_hi);

                h_lo  = @(T) obj.poly_nasa7_enthapy(T, a_lo);
                h_hi  = @(T) obj.poly_nasa7_enthapy(T, a_hi);

                h  = @(T) obj.evaluate(T, tc, h_lo, h_hi);
                cp = @(T) obj.evaluate(T, tc, c_lo, c_hi);

                obj.species{k} = struct("cp_mole", cp, "enthalpy_mole", h);
                obj.species_names{k} = s.name;
                obj.mw(k) = s.molecular_mass;
            endfor

            if strcmp(kin_name, "MAK")
                obj.wdot = @obj.wdot_mak;
            elseif strcmp(kin_name, "EBU")
                obj.wdot = @obj.wdot_ebu;
            else
                error("Invalid kinetics %s", kin_name);
            endif
        endfunction

        function [result] = n_species(self)
            % Access to number of species in mechanism [-].
            result = numel(self.mw);
        endfunction

        function [result] = molecular_masses(self)
            % Access to array of molecular masses [kg/mol].
            result = self.mw;
        endfunction

        function [mw] = mean_molecular_mass_xfrac(self, X)
            % Mixture mean molecular mass from mole fractions [kg/mol].
            mw = X * self.mw';
        endfunction

        function [mw] = mean_molecular_mass_yfrac(self, Y)
            % Mixture mean molecular mass from mass fractions [kg/mol].
            mw = (sum(Y' ./ self.mw').^(-1))';
        endfunction

        function [Y] = mole_to_mass_fraction(self, X)
            % Convert mole to mass fractions.
            Y = X .* self.mw ./ self.mean_molecular_mass_xfrac(X);
        endfunction

        function [X] = mass_to_mole_fraction(self, Y)
            % Convert mass to mole fractions.
            X = ((Y .* self.mean_molecular_mass_yfrac(Y))' ./ self.mw')';
        endfunction

        function [rho] = density_mass(self, T, P, Y)
            % Mixture specific mass [kg/m³].
            m = self.mean_molecular_mass_yfrac(Y);
            rho = P .* m ./ (Thermodata.GAS_CONSTANT .* T);
        endfunction

        function [cp] = specific_heat_mass(self, T, Y)
            % Mixture mass-averaged specific heat [J/(kg.K)].
            cp = 0.0;
            for k=1:self.n_species
                cp = cp + self.species{k}.cp_mole(T) .* Y(:, k) ./ self.mw(k);
            endfor
        endfunction

        function [h] = enthalpy_mass(self, T, Y)
            % Mixture mass-averaged enthalpy [J/kg].
            h = sum((Y .* self.enthalpies_mass(T))')';
        endfunction

        function [hs] = enthalpies_mass(self, T)
            % Matrix of species enthalpies [J/kg].
            hs = [];
            for k=1:self.n_species
                hs = horzcat(hs, self.species{k}.enthalpy_mole(T) ./ self.mw(k));
            endfor
        endfunction

        function [mu] = viscosity(self, T, Y)
            % Gas molecular viscosity [Pa.s].
            mu = 1.0e-05 * (0.1672 * sqrt(T) - 1.058);
        endfunction

        function [k] = thermal_conductivity(self, T, Y)
            % Gas thermal conductivity [W/(m³.K)].
            k = 1.581e-17;
            k = T .* k - 9.463e-14;
            k = T .* k + 2.202e-10;
            k = T .* k - 2.377e-07;
            k = T .* k + 1.709e-04;
            k = T .* k - 7.494e-03;
        endfunction

        function [wdot] = wdot_mak(self, z, T, Y, L)
            % Mass action kinetics methane combustion rate [kg/(m³.s)].
            k0 = 1.1e+07;
            Ea = 83680.0;

            X = self.mass_to_mole_fraction(Y);
            C = (X * self.PRESSURE ./ (self.GAS_CONSTANT .* T));
            k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
            rt = k .* C(:, 1) .* C(:, 2).^0.5;

            wdot = rt * (self.mw .* self.SPECIES_COEFS);
        endfunction

        function [wdot] = wdot_ebu(self, z, T, Y, L)
            % Eddy break-up kinetics methane combustion rate [kg/(m³.s)].
            cr = 4.000e+00;
            bo = 4.375e+00;
            k0 = 1.600e+10;
            Ea = 1.081e+05;

            k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
            rho = self.density_mass(T, self.PRESSURE, Y);
            yf = Y(:, 1);
            yo = Y(:, 2);

            % TODO implement this in ProjectData
            ke = z ./ L;

            R_ebu = (rho.^1) .* cr .* ke .* min(yf, yo ./ bo);
            R_arr = (rho.^2) .* yf .* yo .* k;

            rt = min(R_ebu, R_arr) / self.mw(1);

            wdot = rt * (self.mw .* self.SPECIES_COEFS);
        endfunction

        function hdot = heat_release_rate(self, h, mdotk)
            % Heat release rate [W/m³].
            hdot = sum((mdotk .* h)')';
        endfunction

        function [res] = equilibrate_mak(self, x)
            % Equilibrate solution state.
            error("function not finished");
            x = x';
            T = x(1);
            Y = x(2:end);

            h = self.enthalpies_mass(T);
            wdot = self.wdot_mak(0, T, Y, 0);
            hdot = self.heat_release_rate(h, wdot);

            res = vertcat(hdot, wdot(1:end-1)', sum(Y)-1);
        endfunction
    endmethods

    properties (Constant, Access = public, Static = true)
        % Species stoichiometric coefficients.
        SPECIES_COEFS = [-1.0, -2.0, 1.0, 2.0, 0.0, 0.0];

        % Operating pressure [Pa].
        PRESSURE = 101325.0;
    endproperties

    properties (SetAccess = private, GetAccess = public)
        % Array of mixture species.
        species;

        % Array of mixture species names.
        species_names;

        % Array of species  molecular weight.
        mw;

        % Function handle to kinetic rates.
        wdot;
    endproperties
endclassdef
