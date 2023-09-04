% Run validation of program.

clc; clear; clear functions;

addpath("../src/octave");
import kiln.*


function [vals_ref] = load_data()
    load("validate.mat");

    vals_ref = struct(
        "phi",   data(:, 1),
        "T",     data(:, 2),
        "cp",    data(:, 3),
        "h",     data(:, 4),
        "qdot",  data(:, 5),
        "mw",    data(:, 6),
        "rho",   data(:, 7),
        "mu",    data(:, 8),
        "k",     data(:, 9),
        "X",     data(:, 10:15),
        "Y",     data(:, 16:21),
        "hm",    data(:, 22:27),
        "wdotk", data(:, 28:33)
    );
endfunction


function plot_viscosity(T, Y, mu_ref, sol)
    mu_sol = sol.viscosity(T, Y);

    h = figure();
    plot(T, mu_ref, "linewidth", 2); hold on;
    plot(T, mu_sol, "linewidth", 2); hold on;
    grid();
    xlabel("Temperature [K]");
    ylabel("Viscosity [Pa.s]");
    legend("Reference", "Computed");
    set(gca, "GridLineStyle", ":");
endfunction


function plot_conductivity_bed(bed)
    % Display of bed thermal conductivity.
    T_arr = [295.15, 779.15, 1086.15, 1290.15, 1484.15];
    k_arr = [0.07, 0.12, 0.16, 0.21, 0.30];

    % IMPORTANT: coating conductivity is set to be twice
    % the bed value arbitrarily.
    T = linspace(300, 1500, 100);
    k = bed.thermal_conductivity(T);

    h = figure();

    subplot(1,1,1);
    plot(T, k, ".", "linewidth", 4); hold on;
    plot(T_arr, k_arr, ".r", "MarkerSize", 20); hold on;
    grid();
    xlabel("Temperature [K]");
    ylabel("Conductivity [W/(m²K)]");
    set(gca, "GridLineStyle", ":");
endfunction


function plot_thermal_conductivity(T, Y, k_ref, sol)
    k_sol = sol.thermal_conductivity(T, Y);

    h = figure();
    plot(T, k_ref, "linewidth", 2); hold on;
    plot(T, k_sol, "linewidth", 2); hold on;
    grid();
    xlabel("Temperature [K]");
    ylabel("Thermal Conductivity [W/(m.K)]");
    legend("Reference", "Computed");
    set(gca, "GridLineStyle", ":");
endfunction


function validate_array(val_ref, val_sol, tol, name)
    residuals = (val_sol - val_ref) ./ (val_ref + eps);
    respected = sum(flatten(residuals) < tol);
    invalid = not(respected == numel(val_ref));

    if (invalid)
        error("Method `%s` fails!", name);
    else
        printf("\nMethod `%s` is validated!\n", name);
    endif
endfunction


function validate_gas()
    P    = 101325.0;
    tol  = 1.0e-06;
    vals_ref = load_data();

    sol = GasMixture(data={}, kin_name="MAK");

    plot_viscosity(vals_ref.T, vals_ref.Y, vals_ref.mu, sol);
    plot_thermal_conductivity(vals_ref.T, vals_ref.Y, vals_ref.k, sol);

    validate_array(
        val_ref = vals_ref.mw,
        val_sol = 1000 * sol.mean_molecular_mass_xfrac(vals_ref.X),
        tol = tol,
        name = "mean_molecular_mass_xfrac"
    );
    validate_array(
        val_ref = vals_ref.mw,
        val_sol = 1000 * sol.mean_molecular_mass_yfrac(vals_ref.Y),
        tol = tol,
        name = "mean_molecular_mass_yfrac"
    );
    validate_array(
        val_ref = vals_ref.Y,
        val_sol = sol.mole_to_mass_fraction(vals_ref.X),
        tol = tol,
        name = "mole_to_mass_fraction"
    );
    validate_array(
        val_ref = vals_ref.X,
        val_sol = sol.mass_to_mole_fraction(vals_ref.Y),
        tol = tol,
        name = "mass_to_mole_fraction"
    );
    validate_array(
        val_ref = vals_ref.rho,
        val_sol = sol.density_mass(vals_ref.T, P, vals_ref.Y),
        tol = tol,
        name = "density_mass"
    );
    validate_array(
        val_ref = vals_ref.cp,
        val_sol = sol.specific_heat_mass(vals_ref.T, vals_ref.Y),
        tol = tol,
        name = "specific_heat_mass"
    );
    validate_array(
        val_ref = vals_ref.h,
        val_sol = sol.enthalpy_mass(vals_ref.T, vals_ref.Y),
        tol = tol,
        name = "enthalpy_mass"
    );
    validate_array(
        val_ref = vals_ref.hm,
        val_sol = sol.enthalpies_mass(vals_ref.T),
        tol = tol,
        name = "enthalpies_mass"
    );
    validate_array(
        val_ref = vals_ref.wdotk,
        val_sol = sol.wdot_mak(0, vals_ref.T, vals_ref.Y, 0),
        tol = tol,
        name = "wdot_mak"
    );
    validate_array(
        val_ref = vals_ref.qdot,
        val_sol = sol.heat_release_rate(vals_ref.hm, vals_ref.wdotk),
        tol = tol,
        name = "heat_release_rate"
    );
endfunction


function validate_bed()
    bed = BedMixture(data={});
    plot_conductivity_bed(bed);
endfunction


function validate_radcal()
    X = [1401.3, 1397.8, 0.709, 0.438;
         1000.0, 1200.0, 0.500, 0.400;
         1100.0, 1150.0, 0.300, 0.500];

    radcal = RotaryKilnModel.get_radcal();
    E = radcal(X);
    eps_g = E(:, 1)'
    abs_g = E(:, 2)'
    % TODO add reference values to test against!
endfunction


function main()
    validate_gas();
    validate_bed();
    validate_radcal();
endfunction


main();

