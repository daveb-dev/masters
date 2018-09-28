 /*------------------------------------------------------------------
 * Brief description of this file:
 *
 * This is an example of how to use QUESO classes and algorithms in order to define and solve a statistical inverse problem (SIP).
 * The SIP consists on calibrating the magnitude 'E', Young's Modulus, a measure of the initial slope of the stress-strain curve for steel. The solution of the SIP is the posterior probability density function (PDF) of 'E'.

 *
 * The code consists of 5 files:
 * - 'main.C' (this file)
 * - 'compute.C' (the driving application code)
 * - 'compute.h'
 * - 'likelihood.C' (necessary for the SIP)
 * - 'likelihood.h'
 *-----------------------------------------------------------------*/

#include <compute.h>

int main(int argc, char* argv[])
{
  // Initialize QUESO environment
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(argv[1],"",NULL);
#endif

  // Call application
  computeE(*env);

  // Finalize QUESO environment
  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
