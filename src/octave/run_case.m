function [data, model] = run_case(case_name, save_as, plot_results=false, devel=false)
    % Simple wrapper for running a standard simulation.
    data = load_json(case_name);
    opts = struct("ipopt", data.ipopt);

    model = RotaryKilnModel(data, devel=devel);
    model.solve_system(opts);
    model.dump_results(save_as, debug=true, extras=true);

    if (plot_results)
        model.plot_results();
    endif
endfunction