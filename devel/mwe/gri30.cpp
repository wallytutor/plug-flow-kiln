// Minimal working example for testing environment build.
#include <cantera/thermo/IdealGasPhase.h>
#include <iostream>
#include <iomanip>

using Cantera::newSolution;
using Cantera::OneAtm;

int main()
try
{
    std::cout << "** MWE of Cantera **" << std::endl;

    auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();

    gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");

    std::cout << gas->meanMolecularWeight() << std::endl;

    return 0;
}
catch (std::exception &err) 
{
    std::cerr << err.what() << std::endl;
}
