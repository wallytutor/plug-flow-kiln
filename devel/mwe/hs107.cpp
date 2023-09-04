// Minimal working example for testing environment build.
#include <casadi/casadi.hpp>
#include <iostream>
#include <iomanip>

using casadi::Dict;
using casadi::DM;
using casadi::Function;
using casadi::NlpBuilder;

int main()
try
{
    std::cout << "** MWE of CasADi **" << std::endl;

    NlpBuilder nl;
    nl.import_nl("hs107.nl");

    Dict opts;
    opts["expand"] = true;
    opts["ipopt.linear_solver"] = "mumps";

    Function solver = nlpsol("nlpsol", "ipopt", nl, opts);
    std::map<std::string, DM> arg;
    std::map<std::string, DM> res;

    arg["lbx"] = nl.x_lb;
    arg["ubx"] = nl.x_ub;
    arg["lbg"] = nl.g_lb;
    arg["ubg"] = nl.g_ub;
    arg["x0"]  = nl.x_init;

    res = solver(arg);

    for (auto&& s : res) {
        std::cout << std::setw(10) << s.first << ": " 
                << std::vector<double>(s.second) 
                << std::endl;
    }

    return 0;
}
catch (std::exception &err) 
{
    std::cerr << err.what() << std::endl;
}
