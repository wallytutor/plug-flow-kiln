// Minimal program to get CPU hash.
// https://stackoverflow.com/questions/6412042
#include "get_unique_id.hpp"

uint16_t get_cpu_hash()
{        
   int cpuinfo[4] = {0, 0, 0, 0};
   __cpuid(cpuinfo, 0);

   uint16_t hash = 0;
   uint16_t* ptr = (uint16_t*)(&cpuinfo[0]);

   for (uint32_t i = 0; i < 8; i++) {
      hash += ptr[i];
   }

   return hash;
}

int main()
{
   std::ofstream wf("rklg.dat", std::ios::out | std::ios::binary);

   if(!wf) {
      std::cout << "Cannot open file!" << std::endl;
      return 1;
   }

   LicenceInputs data;
   data.unique_id = int(get_cpu_hash());
   data.timestamp = std::time(nullptr);

   wf.write((char *) &data, sizeof(LicenceInputs));
   wf.close();

   if(!wf.good()) {
      std::cout << "Error occurred at writing time!" << std::endl;
      return 1;
   }

   return 0;
}
