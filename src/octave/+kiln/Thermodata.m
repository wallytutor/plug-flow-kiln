classdef Thermodata
    properties (Constant, Static = true)
        % Ideal gas constant [J/(mol.K)].
        GAS_CONSTANT = 8.31446261815324;

        % Stefan-Boltzmann constant W/(m².K⁴).
        SIGMA = 5.6703744191844314e-08;

        % Water specific heat [J/(kg.K)].
        WATER_SPECIFIC_HEAT = 4186.0;

        % Water latent heat of evaporation [J/kg].
        WATER_LATENT_HEAT_EVAP = 2260000.0;
    end

    methods (Access = public, Static = true)
        function [data] = load_default_data(name)
            % Load default database for this class.
            file_path = fileparts(mfilename("fullpath"));
            data_file = strcat(file_path, "/data/", name);
            data = load_json(data_file);
        endfunction

        function [p] = poly_nasa7_specific_heat(T, c)
            % Molar specific heat from NASA7 polynomial [J/(mol.K)].
            p = c(1)+T.*(c(2)+T.*(c(3)+T.*(c(4)+c(5).*T)));
            p = Thermodata.GAS_CONSTANT .* p;
        endfunction

        function [p] = poly_nasa7_enthapy(T, c)
            % Molar enthalpy from NASA7 polynomial [J/mol].
            d = c(1:5) ./ linspace(1, 5, 5);
            p = d(1)+T.*(d(2)+T.*(d(3)+T.*(d(4)+d(5).*T)))+c(6)./T;
            p = Thermodata.GAS_CONSTANT .* T .* p;
        end

        function [p] = shomate_specific_heat(T, c)
            % Molar specific heat with Shomate equation [J/(mol.K)].
            p = T.*(c(2)+T.*(c(3)+c(4).*T))+c(5)./T.^2+c(1);
        endfunction

        function [p] = shomate_enthalpy(T, c)
            % Molar enthalpy with Shomate equation [J/mol].
            p = T.*(c(1)+T.*(c(2)/2+T.*(c(3)/3+c(4)/4.*T)))-c(5)./T+c(6)-c(8);
        endfunction

        function [p] = shomate_entropy(T, c)
            % Entropy with Shomate equation [J/K].
            p = c(1).*log(T)+T.*(c(2)+T.*(c(3)/2+c(4)/3.*T))+c(5)./(2.*T.^2)+c(7);
        endfunction

        function [p] = evaluate(T, Tc, f_lo, f_hi)
            % Evaluation of a step-wise function.
            h = (1 + sign(T - Tc)) ./ 2;
            p = (1 - h) .* f_lo(T) + h .* f_hi(T);
        endfunction
    endmethods
endclassdef
