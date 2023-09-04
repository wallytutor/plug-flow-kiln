classdef RotaryKilnModel < handle
    methods (Access = public)
        function [obj] = RotaryKilnModel(data, devel=false)
            % Experimental features flag.
            obj.devel = devel;

            % If true, do not evaluate kinetic equations.
            obj.equilibrated = data.hypothesis.gas_equilibrated;

            % Dimensionless gas film thickness for h_cwb.
            obj.gas_film_thickness = data.hypothesis.gas_film_thickness;

            % Initialize materials/models inside rotary kiln.
            obj.gas = GasMixture(data.materials.gas_data,...
                                 data.hypothesis.gas_kinetics);
            obj.bed = BedMixture(data.materials.bed_data);
            obj.radcal = obj.get_radcal();

            % Setup numerical section.
            obj.nz = data.numerical.number_of_slices;

            % Get shortcuts.
            method = data.numerical.discretization_method;
            geom = data.geometry;
            htp = data.heat_transfer_parameters;

            % Note: do not change order, interdependencies!
            obj.build_main_geometry(geom, method);
            obj.build_operation(data.operation, data.initial_guess);
            obj.build_heat_transfer(geom, htp)
            obj.build_bed_geometry();
            obj.build_heat_geometry();
        endfunction

        function build_main_geometry(self, geom, method)
            % Setup main geometric features of kiln.
            self.slope = geom.kiln_slope * pi / 180;
            self.L = geom.kiln_inner_length;
            self.R = geom.kiln_inner_diameter / 2;
            self.dam_height = geom.kiln_dam_height;
            
            % Handle space discretization.
            switch method
                case "linear"
                    self.z = linspace(0.0, self.L, self.nz)';
                otherwise
                    error("Unknown discretization method %s", method);
            end
        endfunction

        function build_operation(self, operation, guess)
            % Rotation rate is stored internally in [rev/s].
            self.rot_rate = operation.rotation_rate / 60.0;

            % Gas flow rate corresponds to a sum of primary/secondary.
            % TODO currently composition is provided for the sum, but it
            % should be averaged here instead!
            flow_rate = operation.gas_flow_rate + operation.gas_leak_rate;
            feed_rate = operation.bed_feed_rate / 3600.0;

            % Set vectorized flow rates.
            self.mdot_g = flow_rate * ones(size(self.z));
            self.mdot_b = feed_rate * ones(size(self.z));
            
            % Store feed rate internally for Kramer's equation.
            % TODO implement variable flow in bed profile!
            self.feed_rate = feed_rate;
            self.feed_humidity = operation.bed_feed_humidity;
            self.Y_b = zeros(self.nz, self.bed.n_species);
            self.Y_b(:, 1) = self.feed_humidity;
            self.Y_b(:, 2) = 1 - self.feed_humidity;

            % Set initial composition profile.
            x0 = operation.gas_feed_composition';
            self.X_g = x0 .* ones(self.nz, length(x0));
            self.Y_g = self.gas.mole_to_mass_fraction(self.X_g);

            % Set linearized initial guess of temperatures.
            tg0 = operation.gas_feed_temperature;
            tb0 = operation.bed_feed_temperature;
            tg1 = guess.gas_final_temperature;
            tb1 = guess.bed_final_temperature;
            self.T_b = linspace(tb1, tb0, self.nz)';
            self.T_g = linspace(tg0, tg1, self.nz)';
        endfunction

        function build_heat_transfer(self, geom, htp)
            % Setup heat transfer parameters and initial guess.
            % Compute different radii and set conductivities.
            % Get configuration geometrical features.
            t_coat = geom.thickness_coating;
            t_refr = geom.thickness_refractory;
            t_shell = geom.thickness_shell;

            % Get configuration thermal conductivities.
            k_coat = htp.conductivity_coating;
            k_refr = htp.conductivity_refractory;
            k_shell = htp.conductivity_shell;

            % Get external HTC.
            h_env = htp.environment_htc;
            T_env = htp.environment_temperature;

            % Make a single ones array.
            ones_z = ones(size(self.z));

            % Set heat transfer parameters
            self.eps_bed = htp.emissivity_bed;
            self.eps_ref = htp.emissivity_refractory;
            self.eps_env = htp.emissivity_shell;
            
            thickness_coat = self.get_function(               ...
                cond  = (t_coat < 0),                         ...
                fun_a = @ProjectData.thickness_coating,       ...
                fun_b = @(z) t_coat .* ones_z                 ...
            );
            thickness_refr = self.get_function(               ...
                cond  = (t_refr < 0),                         ...
                fun_a = @ProjectData.thickness_refractory,    ...
                fun_b = @(z) t_refr .* ones_z                 ...
            );
            thickness_shell = self.get_function(              ...
                cond  = (t_shell < 0),                        ...
                fun_a = @ProjectData.thickness_shell,         ...
                fun_b = @(z) t_shell .* ones_z                ...
            );

            k_coat = self.get_function(                       ...
                cond  = (k_coat < 0),                         ...
                fun_a = @ProjectData.conductivity_coating,    ...
                fun_b = @(z, T) k_coat                        ...
            );
            k_refr = self.get_function(                       ...
                cond  = (k_refr < 0),                         ...
                fun_a = @ProjectData.conductivity_refractory, ...
                fun_b = @(z, T) k_refr                        ...
            );
            k_shell = self.get_function(                      ...
                cond  = (k_shell < 0),                        ...
                fun_a = @ProjectData.conductivity_shell,      ...
                fun_b = @(z, T) k_shell                       ...
            );

            h_env = self.get_function(                        ...
                cond  = (h_env < 0),                          ...
                fun_a = @ProjectData.environment_htc,         ...
                fun_b = @(z) h_env * ones_z                   ...
            );
            T_env = self.get_function(                        ...
                cond  = (T_env < 0),                          ...
                fun_a = @ProjectData.environment_temperature, ...
                fun_b = @(z) T_env * ones_z                   ...
            );

            self.R_cv = self.R - thickness_coat(self.z);
            self.R_cr = self.R .* ones_z;
            self.R_rs = self.R_cr + thickness_refr(self.z);
            self.R_sh = self.R_rs + thickness_shell(self.z);

            self.k_coat = @(T) k_coat(self.z, T);
            self.k_refr = @(T) k_refr(self.z, T);
            self.k_shell = @(T) k_shell(self.z, T);

            self.h_env = h_env(self.z);
            self.T_env = T_env(self.z);

            % Initial guesses for shell properties.
            self.T_wg = 0.9 * self.T_g;
            self.T_cr = 0.8 * self.T_g;
            self.T_rs = 1.1 * self.T_env;
            self.T_sh = 1.0 * self.T_env;
        endfunction

        function build_bed_geometry(self)
            % Initialize bed geometry of rotary kiln.
            % TODO: consider material loss in the future!
            % TODO: odeset, values below are quite arbitrary.
            Rdot = -diff(self.R_cv) ./ diff(self.z);
            alpha = atan(cat(1, Rdot, Rdot(end)));

            alpha = @(x) interp1(self.z, alpha, x, "linear");
            radius = @(x) interp1(self.z, self.R_cv, x, "linear");

            [x, h] = ode15s(                                     ...
                @(z, h) self.kramers_model(z, h, alpha, radius), ...
                [0, self.z(end)],                                ...
                self.dam_height,                                              ...
                odeset('MaxStep', 0.1)                           ...
            );

            h = interp1(x, h, self.z, "linear");

            phi = 2 * acos(1 - h ./ self.R_cv);
            lgb = 2 .* self.R_cv .* sin(phi / 2);

            Ab = (1/2)*(phi .* self.R_cv.^2 - lgb .* (self.R_cv - h));
            Ag = pi * self.R_cv.^2 - Ab;

            eta_loc = (phi - sin(phi)) / (2 * pi);
            eta_bar = trapz(self.z, eta_loc) / self.L;

            % Kiln equivalent diameter used in Tschengs HTC self.
            fn = @(n, b) (3 - n) * pi - b / n + sin(b / n);
            de = 0.5 .* (2 .* self.R_cv) .* fn(1, phi) ./ fn(2, phi);

            self.bed_height = h;
            self.bed_cord_length = lgb;
            self.bed_cross_area = Ab;
            self.gas_cross_area = Ag;
            self.central_angle = phi;
            self.diameter_eff = de;
            self.local_loading = eta_loc;
            self.mean_loading = eta_bar;

            % TODO: here we assume equal spacing in cell_length
            % and that should not be the case because geometric
            % grids are more appropriate in general!
            self.cell_length = self.z(2) - self.z(1);
        endfunction

        function build_heat_geometry(self)
            % Initialize kiln heat transfer geometrical features.

            % Areas for pair gas-wall: since the bed covers the central
            % angle, theta=2pi-phi and areas are computed as follows.
            self.P_cgw = (2 * pi - self.central_angle) .* self.R_cv;
            self.P_rgw = (2 * pi - self.central_angle) .* self.R_cv;
            self.A_cgw = self.P_cgw .* self.cell_length;
            self.A_rgw = self.A_cgw;

            % Areas for pair gas-bed: this area is as trapezoidal
            % section because in fact bed is an inclined plane. Here,
            % to avoid useless complications it is computed as a
            % rectangle of exposed bed.
            self.P_cgb = self.bed_cord_length;
            self.P_rgb = self.bed_cord_length;
            self.A_cgb = self.P_cgb .* self.cell_length;
            self.A_rgb = self.A_cgb;

            % Areas for pair wall-bed: in this case there are two
            % different areas because radiation comes through the
            % exposed surface and conduction from contact with wall.
            % XXX: should'n the CWB be based on emitting surface?
            self.P_cwb = self.central_angle .* self.R_cv;
            self.P_rwb = self.bed_cord_length;
            self.A_cwb = self.P_cwb .* self.cell_length;
            self.A_rwb = self.P_rwb .* self.cell_length;

            % External shell area.
            self.P_env = 2 * pi .* self.R_sh;
            self.A_env = self.P_env .* self.cell_length;

            % View factor for RWB is the ratio of receiving bed
            % surface to emitting walls (interface gas-wall).
            self.omega = self.P_rwb ./ self.P_rgw;

            % Gorog's optical beam correlation used by Hanein (2016).
            D = 2 .* self.R_cv;
            self.beam = 0.95 .* D .* (1 - self.bed_height ./ D);

            % Effective area for radiative heat transfer.
            a1 = (1.0 - self.eps_ref) ./ (self.eps_ref .* self.A_rgw);
            a2 = (1.0 ./ (1.0 .* self.A_rwb));
            a3 = (1.0 - self.eps_bed) ./ (self.eps_bed .* self.A_rwb);
            self.A_rwb_eff = 1.0 ./ (a1 + a2 + a3);
        endfunction

        function unpack_parameters_walls(self, x, finished=false)
            % Unpack inputs/outputs for wall constraints.
            nz = self.nz;
            self.T_wg = x(1 + 0 * nz:1 * nz);
            self.T_cr = x(1 + 1 * nz:2 * nz);
            self.T_rs = x(1 + 2 * nz:3 * nz);
            self.T_sh = x(1 + 3 * nz:4 * nz);

            if (finished)
                self.residuals = self.lhs_walls(x);
            endif
        endfunction

        function unpack_parameters_internal(self, x, finished=false)
            % Unpack inputs/outputs for internal constraints.
            nz = self.nz;
            ng = self.gas.n_species();
            nb = self.bed.n_species();

            self.mdot_g = x(1 + 0 * nz:1 * nz);
            self.mdot_b = x(1 + 1 * nz:2 * nz);
            self.T_g    = x(1 + 2 * nz:3 * nz);
            self.T_b    = x(1 + 3 * nz:4 * nz);

            i = 1 + 4 * nz;
            j = i + ng * nz - 1;
            self.Y_g    = reshape(x(i:j), ng, nz)';

            i = 1 + j;
            j = i + nb * nz - 1;
            self.Y_b    = reshape(x(i:j), nb, nz)';

            self.X_g    = self.gas.mass_to_mole_fraction(self.Y_g);

            if (finished)
                self.residuals = self.lhs_internal(x);
            endif
        endfunction

        function unpack_parameters_system(self, x, finished=false)
            % Unpack inputs/outputs for system constraints.
            nz = self.nz;
            ng = self.gas.n_species();
            nb = self.bed.n_species();

            self.T_wg   = x(1 + 0 * nz:1 * nz);
            self.T_cr   = x(1 + 1 * nz:2 * nz);
            self.T_rs   = x(1 + 2 * nz:3 * nz);
            self.T_sh   = x(1 + 3 * nz:4 * nz);
            self.mdot_g = x(1 + 4 * nz:5 * nz);
            self.mdot_b = x(1 + 5 * nz:6 * nz);
            self.T_g    = x(1 + 6 * nz:7 * nz);
            self.T_b    = x(1 + 7 * nz:8 * nz);

            i = 1 + 8 * nz;
            j = i + ng * nz - 1;
            self.Y_g    = reshape(x(i:j), ng, nz)';

            i = 1 + j;
            j = i + nb * nz - 1;
            self.Y_b    = reshape(x(i:j), nb, nz)';

            self.X_g    = self.gas.mass_to_mole_fraction(self.Y_g);

            if (finished)
                self.residuals = self.lhs_system(x);
                x = abs(self.residuals);

                % Full system solved.
                self.res_T_wg   = x(1 + 0 * nz:1 * nz);
                self.res_T_cr   = x(1 + 1 * nz:2 * nz);
                self.res_T_rs   = x(1 + 2 * nz:3 * nz);
                self.res_T_sh   = x(1 + 3 * nz:4 * nz);

                % One extremity (depending on phase) ignored.
                self.res_mdot_g = x(1 + 4 * nz - 0:5 * nz - 1);
                self.res_mdot_b = x(1 + 5 * nz - 1:6 * nz - 2);
                self.res_T_g    = x(1 + 6 * nz - 2:7 * nz - 3);
                self.res_T_b    = x(1 + 7 * nz - 3:8 * nz - 4);

                i = 1 + 8 * nz - 4;
                j = i + ng * (nz - 1) - 1;
                self.res_Y_g    = reshape(x(i:j), ng, nz - 1)';

                i = 1 + j;
                j = i + nb * (nz - 1) - 1;
                self.res_Y_b    = reshape(x(i:j), nb, nz - 1)';
            endif
        endfunction

        function update_htc(self)
            % Update different heat transfer coefficients in kiln.
            P_g = self.gas.PRESSURE;

            % Evaluate gas and bed properties.
            rho_g = self.gas.density_mass(self.T_g, P_g, self.Y_g);
            mu_g = self.gas.viscosity(self.T_g, self.Y_g);
            k_g = self.gas.thermal_conductivity(self.T_g, self.Y_g);
            k_s = self.bed.thermal_conductivity(self.T_b);

            % Evaluate effective conductivity and associated diffusivity.
            k_b = self.effective_thermal_conductivity(k_g, k_s);

            % Compute gas speed and kiln angular speed.
            u = self.mdot_g ./ (rho_g .* self.gas_cross_area);
            w = 2 * pi * self.rot_rate;

            % Evaluate Re numbers only once.
            de = self.diameter_eff;
            re_d = rho_g .* u .* de ./ mu_g;
            re_w = rho_g .* w .* de .* de ./ mu_g;

            % Compute HTC for all modes/couples.
            self.h_cgb = self.htc_cgb_tscheng(k_g, re_d, re_w);
            self.h_cgw = self.htc_cgw_tscheng(k_g, re_d, re_w);
            self.h_cwb = self.htc_cwb_hanein(k_g, k_b);
        endfunction

        function update_radiative_properties(self)
            % Update gas radiative properties.
            % TODO train radcal with mass fractions to avoid conversion!
            self.X_g = self.gas.mass_to_mole_fraction(self.Y_g);

            X_co2 = self.X_g(:, 3);
            X_h2o = self.X_g(:, 4);

            T_g = self.T_g;
            T_w = self.A_rgb .* self.T_b + self.A_rgw .* self.T_wg;
            T_w = T_w ./ (self.A_rgb + self.A_rgw);

            SMALL = 1.0e-10;

            pg = X_co2 + X_h2o;
            pgl = pg .* self.beam;
            xco2 = X_co2 ./ (pg + SMALL);

            M = horzcat(T_w, T_g, pgl, xco2);
            E = self.radcal(M);

            self.eps_g = E(:, 1);
            self.abs_g = E(:, 2);
        endfunction

        function compute_external_exchanges(self)
            % Compute wall terms heat exchanges fluxes.
            self.q_env   = self.fn_q_env(self.T_sh);
            self.q_shell = self.fn_q_shell(self.T_rs, self.T_sh);
            self.q_refr  = self.fn_q_refr(self.T_cr,  self.T_rs);
            self.q_coat  = self.fn_q_coat(self.T_wg,  self.T_cr);
        endfunction

        function compute_internal_exchanges(self)
            % Compute internal heat exchanges fluxes.
            self.update_htc();
            self.update_radiative_properties();
            self.q_cgw = self.fn_q_cgw(self.T_g,  self.T_wg);
            self.q_rgw = self.fn_q_rgw(self.T_g,  self.T_wg);
            self.q_cgb = self.fn_q_cgb(self.T_g,  self.T_b);
            self.q_rgb = self.fn_q_rgb(self.T_g,  self.T_b);
            self.q_cwb = self.fn_q_cwb(self.T_wg, self.T_b);
            self.q_rwb = self.fn_q_rwb(self.T_wg, self.T_b);
        endfunction

        function heat_balances(self)
            % Compute energy balance on a phase basis.
            self.bal_g = self.q_cgw + self.q_cgb + self.q_rgw + self.q_rgb;
            self.bal_b = self.q_cwb + self.q_cgb + self.q_rwb + self.q_rgb;
            self.bal_c = self.q_cgw + self.q_rgw - self.q_rwb - self.q_cwb;
        endfunction

        function rhs_internal(self)
            % Compute RHS of internal ODE system.
            A_g = self.gas_cross_area;
            A_b = self.bed_cross_area;
            P_g = self.bed_cord_length;
            P_b = self.bed_cord_length;

            % Get gas rates.
            if not(self.equilibrated)
                wdotk_g = self.gas.wdot(self.z, self.T_g, self.Y_g, self.L);
                sdotk_g = zeros(size(wdotk_g));
            else
                wdotk_g = zeros(size(self.Y_g));
                sdotk_g = zeros(size(wdotk_g));
            endif

            % Recover energy balances [qdotv is in W/m].
            qdotv_g = -1 * self.bal_g ./ self.cell_length;
            qdotv_b = +1 * self.bal_b ./ self.cell_length;

            % Get bed rates.
            % Heat supplied to each cell for evaporation computation
            % must be the total, thus we use `bal_b`.
            %
            % The values returned by `sdotk` below are given in [kg/s],
            % but the model expect rates given in [kg/(m².s)] for later
            % computations. Thus divide the values by the exposed area
            % of each cell. Notice that the report does not currently
            % formulates the problem this way but in the course of the
            % development it seemed easier to interpret.
            %
            % TODO make report consistent with this text.
            if (self.devel)
                wdotk_b = 0;
                sdotk_b = self.bed.sdotk(self.z, self.T_b, self.Y_b, self.bal_b);
                sdotk_b = sdotk_b ./ (P_b);
                % sdotk_b = sdotk_b ./ (self.cell_length .* P_b);
                % sdotk_b = sdotk_b ./ (self.cell_length);
            else
                wdotk_b = zeros(size(self.Y_b));
                sdotk_b = zeros(size(wdotk_b));
            endif

            % TODO give steam away to the gas.
            % 4 = H2O index for evaporation towards gas.
            % sdotk_g(:, 4) = -1 * wdotk_b;

            % Compute specific heats.
            cp_g = self.gas.specific_heat_mass(self.T_g, self.Y_g);
            cp_b = self.bed.specific_heat_mass(self.T_b, self.Y_b);

            % Gas equations (2).
            h_g = self.gas.enthalpies_mass(self.T_g);
            hdotv_g = -A_g .* self.gas.heat_release_rate(h_g, wdotk_g);
            hdots_g = -P_g .* self.gas.heat_release_rate(h_g, sdotk_g);

            % Gas equations (3).
            Ydotv_g = A_g .* wdotk_g;
            Ydots_g = P_g .* sdotk_g;

            % Bed equations (2).
            % h_b = [self.bed.WATER_SPECIFIC_HEAT * self.T_b,
            %        self.bed.enthalpy_mass(self.T_b)
            %        ];
            hdotv_b = -A_b .* 0;
            hdots_b = -P_b .* 0;
            % hdots_b = -P_b .* sum(sdotk_b')' * self.bed.WATER_LATENT_HEAT_EVAP;

            % Bed equations (3).
            Ydotv_b = A_b .* wdotk_b;
            Ydots_b = P_b .* sdotk_b;

            % Assembly terms into space derivatives.
            sdot_g = P_g .* sum(sdotk_g')';
            sdot_b = P_b .* sum(sdotk_b')';
            Tdot_g = (hdotv_g + hdots_g + qdotv_g) ./ (self.mdot_g .* cp_g);
            Tdot_b = (hdotv_b + hdots_b + qdotv_b) ./ (self.mdot_b .* cp_b);
            Ydot_g = (Ydotv_g + Ydots_g - self.Y_g .* sdot_g) ./ self.mdot_g;
            Ydot_b = (Ydotv_b + Ydots_b - self.Y_b .* sdot_b) ./ self.mdot_b;

            % Select equations to solve: because we are handling counter flows,
            % gas B.C. is located on the first index (which should not be solved
            % for) and bed B.C. on the last index, also not solved. To assembly
            % the proper optimization problem to solve steady state one should
            % take care when filtering the derivatives. Since derivatives are
            % *plug-like*, forward propagation is used, e.g. xdot(1) is used
            % with the finite difference x(2) - x(1) = dz * xdot(1) for gas and
            % and analogous reversed for the bed.
            self.sdot_g = sdot_g(1:end-1);
            self.Tdot_g = Tdot_g(1:end-1);
            self.Ydot_g = Ydot_g(1:end-1, :);
            self.sdot_b = sdot_b(2:end);
            self.Tdot_b = Tdot_b(2:end);
            self.Ydot_b = Ydot_b(2:end, :);

            self.Ydot_g = flatten(self.Ydot_g);
            self.Ydot_b = flatten(self.Ydot_b);
        endfunction

        function rhs_finite_differences(self)
            % Compute forward gradient by finite differences.
            %
            % Notice here that *forward* derivatives of bed are reversed (flux
            % from right to left), as illystrated by the samples below:
            % e.g. xdot(1) = x(1) - x(2), xdot(n-1) = x(n-1) - x(n).
            self.dm_g = self.mdot_g(2:end)   - self.mdot_g(1:end-1);
            self.dT_g = self.T_g(2:end)      - self.T_g(1:end-1);
            self.dY_g = self.Y_g(2:end, :)   - self.Y_g(1:end-1, :);
            self.dm_b = self.mdot_b(1:end-1) - self.mdot_b(2:end);
            self.dT_b = self.T_b(1:end-1)    - self.T_b(2:end);
            self.dY_b = self.Y_b(1:end-1, :) - self.Y_b(2:end, :);

            self.dY_g = flatten(self.dY_g);
            self.dY_b = flatten(self.dY_b);
        endfunction

        function [hdot] = kramers_model(self, z, h, alpha, radius)
            rho = self.bed.density_mass();
            repose = self.bed.repose_angle();
            vdot = self.feed_rate / rho;

            R = radius(z);
            alpha = alpha(z);

            phi = (3/4) * vdot / (pi * self.rot_rate * R^3);
            terml = tan(self.slope + alpha) / sin(repose);
            termr = phi * ((2 - h / R) * h / R)^(-3/2);

            hdot = -tan(repose) * (terml - termr);
        endfunction

        function [k_eff] = effective_thermal_conductivity(self, k_g, k_s)
            % Maxwell effective medium theory approximation.
            % NOTE: at first I had understood wrong on Hanein's paper. In
            % this function `phi` is actually the solid packing!
            phi = self.bed.solid_packing;
            f_sum = 2 * k_g + k_s;
            f_dif = k_s - k_g;
            num = f_sum + 2 * phi .* f_dif;
            den = f_sum - 1 * phi .* f_dif;
            k_eff = (num ./ den) .* k_g;
        endfunction

        function [htc] = htc_cgb_tscheng(self, k_g, re_d, re_w)
            n = self.local_loading;
            de = self.diameter_eff;
            a = [+0.46, +0.535, +0.104, -0.341];
            nu = a(1) .* re_d.^a(2) .* re_w.^a(3) .* n.^a(4);
            htc = (k_g ./ de) .* nu;
        endfunction

        function [htc] = htc_cgw_tscheng(self, k_g, re_d, re_w)
            n = self.local_loading;
            de = self.diameter_eff;
            a = [+1.54, +0.575, -0.292, +0.000];
            nu = a(1) .* re_d.^a(2) .* re_w.^a(3) .* n.^a(4);
            htc = (k_g ./ de) .* nu;
        endfunction

        function [htc] = htc_cwb_hanein(self, k_g, k_b, w)
            % Note: the equation is wrong in Hanein's paper (unit of
            % rotation rate, it is not inversed but actually 1/htc,
            % definition of theta not clear,...). Checking the original
            % Li (2005) allowed to find the expected values!
            chi = self.gas_film_thickness;
            d_p = self.bed.particle_diameter;
            rho_b = self.bed.density_mass;
            cp_b = self.bed.specific_heat_mass(self.T_b, self.Y_b);
            theta = self.central_angle / 2;
            n = 60 .* self.rot_rate;

            term1 = (chi * d_p) ./ k_g;
            term2 = 2 .* k_b .* rho_b .* cp_b .* n ./ theta;

            htc = 1 ./ (term1 + 0.5 ./ sqrt(term2));
        endfunction

        function [qdot] = fn_q_coat(self, T_wg, T_cr)
            % Heat flux accross internal coating [W].
            km = self.k_coat((T_wg + T_cr) / 2);
            den = 2 * pi .* self.cell_length .* km .* (T_wg - T_cr);
            qdot = den ./ log(self.R_cr ./ self.R_cv);
        endfunction

        function [qdot] = fn_q_refr(self, T_cr, T_rs)
            % Heat flux accross refractory [W].
            km = self.k_refr((T_cr + T_rs) / 2);
            den = 2 * pi .* self.cell_length .* km .* (T_cr - T_rs);
            qdot = den ./ log(self.R_rs ./ self.R_cr);
        endfunction

        function [qdot] = fn_q_shell(self, T_rs, T_sh)
            % Heat flux accross shell [W].
            km = self.k_shell((T_rs + T_sh) / 2);
            den = 2 * pi .* self.cell_length .* km .* (T_rs - T_sh);
            qdot = den ./ log(self.R_sh ./ self.R_rs);
        endfunction

        function [qdot] = fn_q_env(self, T_sh)
            % Heat flux towards environment [W].
            con = self.h_env .* (T_sh - self.T_env);
            rad = self.eps_env .* (T_sh.^4 - self.T_env.^4);
            qdot = self.A_env .* (con + Thermodata.SIGMA .* rad);
        endfunction

        function [qdot] = fn_q_cgw(self, T_g, T_w)
            % Convection from gas to wall Eq. (18).
            qdot = self.h_cgw .* self.A_cgw .* (T_g - T_w);
        endfunction

        function [qdot] = fn_q_rgw(self, T_g, T_w)
            % Radiation from gas to wall Eq. (20).
            E = (1.0 + self.eps_ref) ./ 2.0;
            A = self.A_rgw;
            eu = self.eps_g;
            au = self.abs_g;
            qdot = Thermodata.SIGMA .* E .* A .* (eu.*T_g.^4 - au.*T_w.^4);
        endfunction

        function [qdot] = fn_q_cgb(self, T_g, T_b)
            % Convection from gas to bed Eq. (18).
            qdot = self.h_cgb .* self.A_cgb .* (T_g - T_b);
        endfunction

        function [qdot] = fn_q_rgb(self, T_g, T_b)
            % Radiation from gas to bed Eq. (20).
            E = (1.0 + self.eps_bed) / 2.0;
            A = self.A_rgb;
            eu = self.eps_g;
            au = self.abs_g;
            qdot = Thermodata.SIGMA .* E .* A .* (eu.*T_g.^4 - au.*T_b.^4);
        endfunction

        function [qdot] = fn_q_cwb(self, T_w, T_b)
            % Conduction(-like) from wall to bed Eq. (22).
            qdot = self.h_cwb .* self.A_cwb .* (T_w - T_b);
        endfunction

        function [qdot] = fn_q_rwb(self, T_w, T_b)
            % Radiation from wall to bed Eq. (21).
            qdot = Thermodata.SIGMA .* self.A_rwb_eff .* (T_w.^4 - T_b.^4);
        endfunction

        function [lhs] = lhs_walls(self, x)
            % Left-hand side of walls system model.
            self.unpack_parameters_walls(x);
            self.compute_external_exchanges();
            self.compute_internal_exchanges();
            self.heat_balances();

            lhs = vertcat(                                      ...
                self.q_coat  - self.bal_c,                      ...
                self.q_refr  - self.q_coat,                     ...
                self.q_shell - self.q_refr,                     ...
                self.q_env   - self.q_shell                     ...
            );
        endfunction

        function [lhs] = lhs_internal(self, x)
            % Left-hand side of internal system model.
            self.unpack_parameters_internal(x);
            self.compute_internal_exchanges();
            self.heat_balances();
            self.rhs_internal();
            self.rhs_finite_differences();

            lhs = vertcat(                                      ...
                self.dm_g    - self.cell_length .* self.sdot_g, ...
                self.dm_b    - self.cell_length .* self.sdot_b, ...
                self.dT_g    - self.cell_length .* self.Tdot_g, ...
                self.dT_b    - self.cell_length .* self.Tdot_b, ...
                self.dY_g    - self.cell_length .* self.Ydot_g, ...
                self.dY_b    - self.cell_length .* self.Ydot_b  ...
            );
        endfunction

        function [lhs] = lhs_system(self, x)
            % Left-hand side of coupled system model.
            self.unpack_parameters_system(x);
            self.compute_external_exchanges();
            self.compute_internal_exchanges();
            self.heat_balances();
            self.rhs_internal();
            self.rhs_finite_differences();

            lhs = vertcat(                                      ...
                self.q_coat  - self.bal_c,                      ...
                self.q_refr  - self.q_coat,                     ...
                self.q_shell - self.q_refr,                     ...
                self.q_env   - self.q_shell,                    ...
                self.dm_g    - self.cell_length .* self.sdot_g, ...
                self.dm_b    - self.cell_length .* self.sdot_b, ...
                self.dT_g    - self.cell_length .* self.Tdot_g, ...
                self.dT_b    - self.cell_length .* self.Tdot_b, ...
                self.dY_g    - self.cell_length .* self.Ydot_g, ...
                self.dY_b    - self.cell_length .* self.Ydot_b  ...
            );
        endfunction

        function solve_system(self, opts)
            tic();

            % Initialize arrays of bounds.
            % TODO parametrize mass boundaries.
            mdot_max = self.mdot_g + self.mdot_b;
            mdot_min = 0.000 * ones(self.nz, 1);
            mdot_max = 5.000 * ones(self.nz, 1);
            T_min = 200.0 .* ones(self.nz, 1);
            T_max = 3000.0 .* ones(self.nz, 1);
            Y_min_g = zeros(self.nz*self.gas.n_species, 1);
            Y_max_g = ones(self.nz*self.gas.n_species, 1);
            Y_min_b = zeros(self.nz*self.bed.n_species, 1);
            Y_max_b = ones(self.nz*self.bed.n_species, 1);

            % Copy arrays for specific modifications.
            mdot_min_g = mdot_min;
            mdot_min_b = mdot_min;
            mdot_max_g = mdot_max;
            mdot_max_b = mdot_max;
            T_min_g = T_min;
            T_min_b = T_min;
            T_max_g = T_max;
            T_max_b = T_max;

            % Enforce B.C. in phase-specific arrays.
            mdot_min_g(1)      = self.mdot_g(1);
            mdot_max_g(1)      = self.mdot_g(1);
            T_min_g(1)         = self.T_g(1);
            T_max_g(1)         = self.T_g(1);
            mdot_min_b(end)    = self.mdot_b(end);
            mdot_max_b(end)    = self.mdot_b(end);
            T_min_b(end)       = self.T_b(end);
            T_max_b(end)       = self.T_b(end);
            Y_min_g(1:6)       = self.Y_g(1, :);
            Y_max_g(1:6)       = self.Y_g(1, :);

            % if (self.devel)
            %     Y_min_b(end-1:end) = self.Y_b(end, :);
            %     Y_max_b(end-1:end) = self.Y_b(end, :);
            % endif

            M = [self.T_wg,          T_min,       T_max;      ...
                 self.T_cr,          T_min,       T_max;      ...
                 self.T_rs,          T_min,       T_max;      ...
                 self.T_sh,          T_min,       T_max;      ...
                 self.mdot_g,        mdot_min_g,  mdot_max_g; ...
                 self.mdot_b,        mdot_min_b,  mdot_max_b; ...
                 self.T_g,           T_min_g,     T_max_g;    ...
                 self.T_b,           T_min_b,     T_max_b;    ...
                 flatten(self.Y_g),  Y_min_g,       Y_max_g;  ...
                 flatten(self.Y_b),  Y_min_b,       Y_max_b;  ...
            ];

            results = self.solve_ipopt(
                x0   = M(:, 1),
                g    = @self.lhs_system,
                lbx  = M(:, 2),
                ubx  = M(:, 3),
                opts = opts
            );
            self.unpack_parameters_system(results, finished=true);

            printf("\nSimulation took %.2f s\n", toc());
        endfunction

        function [h] = plot_results(self)
            % Plot all results/references together.
            tc = (self.R_cr - self.R_cv) .* 1000.0;

            baseline = tan(self.slope) * self.z;

            b = self.bed_height;
            R = self.R_cv;
            r = self.R_cr;
            eta_loc = 100 * self.local_loading;
            eta_bar = 100 * self.mean_loading * ones(size(self.z));

            % Integrate heat balances.
            bal_g = -1 .* cumtrapz(self.bal_g) ./ 1000;
            bal_b = +1 .* cumtrapz(self.bal_b) ./ 1000;
            bal_c = +1 .* cumtrapz(self.bal_c) ./ 1000;
            bal_b = flip(bal_b);

            % Heat flux density [kW/m].
            q_cgw = self.q_cgw ./ (1000 .* self.cell_length);
            q_cgb = self.q_cgb ./ (1000 .* self.cell_length);
            q_cwb = self.q_cwb ./ (1000 .* self.cell_length);
            q_rgw = self.q_rgw ./ (1000 .* self.cell_length);
            q_rgb = self.q_rgb ./ (1000 .* self.cell_length);
            q_rwb = self.q_rwb ./ (1000 .* self.cell_length);

            % Compute bed residence time.
            tau_bed = self.bed.density_mass * self.bed_cross_area;
            tau_bed = tau_bed .* self.cell_length ./ self.mdot_b;
            tau_bed = cumtrapz(tau_bed) / 60.0;
            tau_bed = tau_bed(end) - tau_bed;

            % Reversed axis for matching Hanein paper.
            z = self.L - self.z;

            h = figure("units", "normalized", "outerposition",...
                       [0.05 0.07 0.9 0.9]);

            %% Column 1

            subplot(3, 3, 1);
            plot(z, tc, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Thickness [mm]");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            subplot(3, 3, 4);
            plot(z, baseline-R+b, "linewidth", 2); hold on;
            plot(z, baseline-R+0, "linewidth", 2); hold on;
            plot(z, baseline-r+0, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Diameter [m]");
            legend("Bed", "Coating", "Refractory");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            subplot(3, 3, 7);
            plot(z, eta_loc, "linewidth", 2); hold on;
            plot(z, eta_bar, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Loading [%]");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            %% Column 2

            subplot(3, 3, 2);
            plot(z, self.T_g  - 273.15, "linewidth", 2); hold on;
            plot(z, self.T_b  - 273.15, "linewidth", 2); hold on;
            plot(z, self.T_wg - 273.15, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Temperature [°C]");
            legend("Gas", "Bed", "Wall");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            subplot(3, 3, 5);
            plot(z, self.T_sh - 273.15, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Temperature [°C]");
            legend("Shell");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            subplot(3, 3, 8);
            plot(z, self.X_g(:, 1), "linewidth", 2); hold on;
            plot(z, self.X_g(:, 2), "linewidth", 2); hold on;
            plot(z, self.X_g(:, 3), "linewidth", 2); hold on;
            plot(z, self.X_g(:, 4), "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Mole fraction [-]");
            legend("CH4", "O2", "CO2", "H2O");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            %% Column 3

            subplot(3, 3, 3);
            plot(z, bal_g, "linewidth", 2); hold on;
            plot(z, bal_b, "linewidth", 2); hold on;
            plot(z, bal_c, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Heat flux [kW]");
            legend("Gas", "Bed", "Coating");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            subplot(3, 3, 6);
            plot(z, q_cgw, "linewidth", 2); hold on;
            plot(z, q_cgb, "linewidth", 2); hold on;
            plot(z, q_cwb, "linewidth", 2); hold on;
            plot(z, q_rgw, "linewidth", 2); hold on;
            plot(z, q_rgb, "linewidth", 2); hold on;
            plot(z, q_rwb, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Heat flux [kW/m]");
            legend(
                "Q^{cv}_{g-w}",
                "Q^{cv}_{g-s}", 
                "Q^{cd}_{w-s}", 
                "Q^{rd}_{g-w}", 
                "Q^{rd}_{g-s}", 
                "Q^{rd}_{w-s}"
            );
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");

            subplot(3, 3, 9);
            plot(tau_bed, self.T_b - 273.15, "linewidth", 2); hold on;
            grid();
            xlabel("Residence time [min]");
            ylabel("Temperature [°C]");
            set(gca, "GridLineStyle", ":");
            % set(gca, "xdir", "reverse");
        endfunction

        function [h] = plot_coating_thickness(self)
            % Illustrate coating profile accross kiln length.
            tc = (self.R_cr - self.R_cv) .* 1000.0;

            h = figure();

            subplot(1, 1, 1);
            plot(self.z, tc, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Thickness [mm]");
            set(gca, "GridLineStyle", ":");
        endfunction

        function [h] = plot_bed_profile(self)
            % Illustrate bed profile with kiln slope.
            baseline = tan(self.slope) * self.z;

            b = self.bed_height;
            R = self.R_cv;
            r = self.R_cr;
            eta_loc = 100 * self.local_loading;
            eta_bar = 100 * self.mean_loading * ones(size(self.z));

            h = figure();

            subplot(2, 1, 1);
            plot(self.z, baseline-R+b, "r", "linewidth", 2); hold on;
            plot(self.z, baseline-R+0, "b", "linewidth", 2); hold on;
            plot(self.z, baseline-r+0, "k", "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Diameter [m]");
            set(gca, "GridLineStyle", ":");

            subplot(2, 1, 2);
            plot(self.z, eta_loc, "linewidth", 2); hold on;
            plot(self.z, eta_bar, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Loading [%]");
            set(gca, "GridLineStyle", ":");
        endfunction

        function [h] = plot_temperature_profiles(self)
            h = figure();

            subplot(2, 1, 1);
            plot(self.z, self.T_g  - 273.15, "linewidth", 2); hold on;
            plot(self.z, self.T_b  - 273.15, "linewidth", 2); hold on;
            plot(self.z, self.T_wg - 273.15, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Temperature [°C]");
            legend("Gas", "Bed", "Wall");
            set(gca, "GridLineStyle", ":");

            subplot(2, 1, 2);
            plot(self.z, self.T_sh - 273.15, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Temperature [°C]");
            legend("Shell");
            set(gca, "GridLineStyle", ":");
        endfunction

        function [h] = plot_internal_profiles(self)
            self.X_g = self.gas.mass_to_mole_fraction(self.Y_g);

            h = figure();

            subplot(2, 1, 1);
            plot(self.z, self.T_g  - 273.15, "linewidth", 2); hold on;
            plot(self.z, self.T_b  - 273.15, "linewidth", 2); hold on;
            plot(self.z, self.T_wg - 273.15, "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Temperature [°C]");
            legend("Gas", "Bed", "Wall");
            set(gca, "GridLineStyle", ":");

            subplot(2, 1, 2);
            plot(self.z, self.X_g(:, 1), "linewidth", 2); hold on;
            plot(self.z, self.X_g(:, 2), "linewidth", 2); hold on;
            plot(self.z, self.X_g(:, 3), "linewidth", 2); hold on;
            plot(self.z, self.X_g(:, 4), "linewidth", 2); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Mole fraction [-]");
            legend("CH4", "O2", "CO2", "H2O");
            set(gca, "GridLineStyle", ":");
        endfunction

        function [h] = plot_heat_fluxes(self)
            % Heat flux density [kW/m].
            q_cgw = self.q_cgw ./ (1000 .* self.cell_length);
            q_cgb = self.q_cgb ./ (1000 .* self.cell_length);
            q_cwb = self.q_cwb ./ (1000 .* self.cell_length);
            q_rgw = self.q_rgw ./ (1000 .* self.cell_length);
            q_rgb = self.q_rgb ./ (1000 .* self.cell_length);
            q_rwb = self.q_rwb ./ (1000 .* self.cell_length);

            % Reversed axis for matching Hanein paper.
            z = self.L - self.z;

            h = figure(1, "position", [0, 0, 800, 400]);;

            subplot(1, 1, 1);
            plot(z, q_cgw, "g", "linewidth", 1); hold on;
            plot(z, q_cgb, "b", "linewidth", 1); hold on;
            plot(z, q_cwb, "r", "linewidth", 1); hold on;
            plot(z, q_rgw, "m", "linewidth", 1); hold on;
            plot(z, q_rgb, "k", "linewidth", 1); hold on;
            plot(z, q_rwb, "c", "linewidth", 1); hold on;
            grid();
            xlabel("Position [m]");
            ylabel("Heat flux [kW/m]");
            legend(
                "Q^{cv}_{g-w}",
                "Q^{cv}_{g-s}", 
                "Q^{cd}_{w-s}", 
                "Q^{rd}_{g-w}", 
                "Q^{rd}_{g-s}", 
                "Q^{rd}_{w-s}"
            );
            xlim([0.0 self.L]);
            set(gca, "GridLineStyle", ":");
        endfunction

        function dump_results(self, fname, debug=false, extras=false)
            % Dump results (with optional debug and extras) to file.
            results = struct();
            results.z = self.z;
            results.local_loading = self.local_loading;
            results.mdot_b = self.mdot_b;
            results.mdot_g = self.mdot_g;
            results.T_g = self.T_g;
            results.T_b = self.T_b;
            results.T_wg = self.T_wg;
            results.T_cr = self.T_cr;
            results.T_rs = self.T_rs;
            results.T_sh = self.T_sh;
            results.T_env = self.T_env;

            species_names = self.gas.species_names;

            for i=1:length(species_names)
                Y_name = strcat("Y_", cell2mat(species_names(i)));
                X_name = strcat("X_", cell2mat(species_names(i)));
                results = setfield(results, Y_name, self.Y_g(:, i));
                results = setfield(results, X_name, self.X_g(:, i));
            endfor

            species_names = self.bed.species_names;

            for i=1:length(species_names)
                Y_name = strcat("Y_", cell2mat(species_names(i)));
                results = setfield(results, Y_name, self.Y_b(:, i));
            endfor

            results.h_cgb = self.h_cgb;
            results.h_cgw = self.h_cgw;
            results.h_cwb = self.h_cwb;
            results.h_env = self.h_env;
            results.q_cgw = self.q_cgw;
            results.q_cgb = self.q_cgb;
            results.q_cwb = self.q_cwb;
            results.q_rgw = self.q_rgw;
            results.q_rgb = self.q_rgb;
            results.q_rwb = self.q_rwb;
            results.q_env = self.q_env;
            results.bal_b = self.bal_b;
            results.bal_c = self.bal_c;
            results.bal_g = self.bal_g;
            results.abs_g = self.abs_g;
            results.eps_g = self.eps_g;

            if (debug)
                results.A_cgb = self.A_cgb;
                results.A_cgw = self.A_cgw;
                results.A_cwb = self.A_cwb;
                results.A_env = self.A_env;
                results.A_rgb = self.A_rgb;
                results.A_rgw = self.A_rgw;
                results.A_rwb = self.A_rwb;
                results.A_rwb_eff = self.A_rwb_eff;
                results.P_cgb = self.P_cgb;
                results.P_cgw = self.P_cgw;
                results.P_cwb = self.P_cwb;
                results.P_env = self.P_env;
                results.P_rgb = self.P_rgb;
                results.P_rgw = self.P_rgw;
                results.P_rwb = self.P_rwb;
                results.R_cr = self.R_cr;
                results.R_cv = self.R_cv;
                results.R_rs = self.R_rs;
                results.R_sh = self.R_sh;
                results.bed_cord_length = self.bed_cord_length;
                results.beam = self.beam;
                results.omega = self.omega;
                results.q_coat = self.q_coat;
                results.q_refr = self.q_refr;
                results.q_shell = self.q_shell;
                results.k_coat = self.k_coat(0.5*(self.T_wg + self.T_cr));
                results.k_refr = self.k_refr(0.5*(self.T_cr + self.T_rs));
                results.k_shell = self.k_shell(0.5*(self.T_rs + self.T_sh));

                if not(all(size(results.k_coat) == size(self.z)))
                    results.k_coat = results.k_coat .* ones(size(self.z));
                endif

                if not(all(size(results.k_refr) == size(self.z)))
                    results.k_refr = results.k_refr .* ones(size(self.z));
                endif

                if not(all(size(results.k_shell) == size(self.z)))
                    results.k_shell = results.k_shell .* ones(size(self.z));
                endif
            endif

            if (extras)
                results.sdot_b = cat(1, self.sdot_b, 0);
                results.sdot_g = cat(1, 0, self.sdot_g);
      
                results.Tdot_b = cat(1, self.Tdot_b, 0);
                results.Tdot_g = cat(1, 0, self.Tdot_g);

                nb = self.bed.n_species;
                species_names = self.bed.species_names;
                Ydot_b = reshape(self.Ydot_b, self.nz-1, nb);
                Ydot_b = cat(1, Ydot_b, zeros(1, nb));

                for i=1:length(species_names)
                    Y_name = strcat("Ydot_b_", cell2mat(species_names(i)));
                    results = setfield(results, Y_name, Ydot_b(:, i));
                endfor

                ng = self.gas.n_species;
                species_names = self.gas.species_names;
                Ydot_g = reshape(self.Ydot_g, self.nz-1, ng);
                Ydot_g = cat(1, zeros(1, ng), Ydot_g);

                for i=1:length(species_names)
                    Y_name = strcat("Ydot_g_", cell2mat(species_names(i)));
                    results = setfield(results, Y_name, Ydot_g(:, i));
                endfor
            endif

            dump_json(fname, results);
        endfunction
    endmethods

    methods (Access = public, Static = true)
        function [results] = solve_ipopt(x0, g, lbx, ubx, opts)
            x = SX.sym("x", numel(x0), 1);
            nlp = struct("x", x, "f", 1, "g", g(x));
            solver = nlpsol("solver", "ipopt", nlp, opts);

            results = solver(...
                x0=x0,       ...
                p=[],        ...
                lbx=lbx,     ...
                ubx=ubx,     ...
                lbg=0.0,     ...
                ubg=0.0,     ...
                lam_x=[],    ...
                lam_g=[]     ...
            );

            results = full(results);
        endfunction

        function [x] = relu(x)
            % ReLU function compatible with CasADi differentiation.
            x = (1 + sign(x)) / 2 .* x;
        endfunction

        function [radcal] = get_radcal()
            file_path = fileparts(mfilename("fullpath"));
            neural_file = strcat(file_path, "/data/radcal.json");
            scaler_file = strcat(file_path, "/data/scaler.json");

            scaler = load_json(scaler_file);
            neural = load_json(neural_file);

            mean = scaler.mean;
            dev = sqrt(scaler.var);

            % TODO generalize for any shape!
            m1 = cell2mat(neural.M(1));
            m2 = cell2mat(neural.M(2));
            m3 = cell2mat(neural.M(3));
            m4 = cell2mat(neural.M(4));
            m5 = cell2mat(neural.M(5));

            b1 = cell2mat(neural.b(1));
            b2 = cell2mat(neural.b(2));
            b3 = cell2mat(neural.b(3));
            b4 = cell2mat(neural.b(4));
            b5 = cell2mat(neural.b(5));

            function [z] = radcal_evaluate(x)
                z = (x' - mean) ./ dev;
                z = RotaryKilnModel.relu(m1' * z + b1);
                z = RotaryKilnModel.relu(m2' * z + b2);
                z = RotaryKilnModel.relu(m3' * z + b3);
                z = RotaryKilnModel.relu(m4' * z + b4);
                z = RotaryKilnModel.relu(m5' * z + b5);
                z = max(0.0, min(1.0, z))';
            endfunction

            radcal = @radcal_evaluate;
        endfunction

        function [f] = get_function(cond, fun_a, fun_b)
            if (cond)
                f = fun_a;
            else
                f = fun_b;
            endif
        endfunction
    endmethods

    properties (SetAccess = private, GetAccess = public)
        % Enable/disable experimental features.
        devel;

        % Hypotheses.
        equilibrated;
        gas_film_thickness;

        % Number of cells in system.
        nz;

        % Kiln geometrical features.
        L;
        R;
        z;
        slope;
        dam_height;
        rot_rate;
        feed_rate;

        % Different layers radii.
        R_cv;
        R_cr;
        R_rs;
        R_sh;

        % Problem integration profiles.
        T_wg;
        T_cr;
        T_rs;
        T_sh;
        T_b;
        T_g;
        Y_g;
        Y_b;
        X_g;

        % Internal objects.
        gas;
        bed;
        radcal;

        % Mass fraction of water in feed [-].
        feed_humidity;

        % Bed solution profile data.
        bed_height;
        bed_cord_length;
        bed_cross_area;
        gas_cross_area;
        central_angle;
        diameter_eff;
        local_loading;
        mean_loading;
        cell_length;

        % Vectorized flow rates.
        mdot_g;
        mdot_b;

        % Perimeters and areas for pair GW.
        P_cgw;
        P_rgw;
        A_cgw;
        A_rgw;

        % Perimeters and areas for pair GB.
        P_cgb;
        P_rgb;
        A_cgb;
        A_rgb;

        % Perimeters and areas for pair WB.
        P_cwb;
        P_rwb;
        A_cwb;
        A_rwb;

        % Perimeter and area for shell-environment.
        P_env;
        A_env;

        % Radiative heat transfer parameters
        omega;
        beam;
        A_rwb_eff;

        % Emissivities of bed [-].
        eps_bed;

        % Emissivities of refractory [-].
        eps_ref;

        % Emissivity of environment [-].
        eps_env;

        % Shell heat transfer coefficient [W/(m².K)].
        h_env;

        % Environment temperature [K].
        T_env;

        % Layers thermal conductivities f(z, T).
        k_coat;
        k_refr;
        k_shell;

        % Heat transfer parameters.
        h_cgb;
        h_cgw;
        h_cwb;
        eps_g;
        abs_g;

        % Shell heat exchange fluxes.
        q_env;
        q_shell;
        q_refr;
        q_coat;

        % Internal heat exchange fluxes.
        q_cgw;
        q_rgw;
        q_cgb;
        q_rgb;
        q_cwb;
        q_rwb;

        % Heat balances.
        bal_g;
        bal_b;
        bal_c;

        % Internal derivates.
        sdot_g;
        Tdot_g;
        Ydot_g;
        sdot_b;
        Tdot_b;
        Ydot_b;

        % Finite differences.
        dm_g;
        dT_g;
        dY_g;
        dm_b;
        dT_b;
        dY_b;

        % Final residual
        residuals;
        res_T_wg;
        res_T_cr;
        res_T_rs;
        res_T_sh;
        res_mdot_g;
        res_mdot_b;
        res_T_g;
        res_T_b;
        res_Y_g;
        res_Y_b;
    endproperties
endclassdef
