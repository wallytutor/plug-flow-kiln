// Minimal program to parse dumped CPU hash.
#include "get_unique_id.hpp"

int main()
{
    std::ifstream rf("rklg.dat", std::ios::out | std::ios::binary);

    if(!rf) {
        std::cout << "Cannot open file!" << std::endl;
        return 1;
    }

    LicenceInputs data;
    rf.read((char *) &data, sizeof(LicenceInputs));
    rf.close();

    if(!rf.good()) {
        std::cout << "Error occurred at reading time!" << std::endl;
        return 1;
    }

    std::cout << "Machine details:" << std::endl;
    std::cout << "Unique ID: " << data.unique_id << std::endl;
    std::cout << "Timestamp: " << data.timestamp << std::endl;

    return 0;
}