/*
 * This file handles the statistical inverse problem (SIP) for estimating
 *   the magnitude 'E' of Young's modulus
 *
 * The SIP definition requires a user defined likelihood function; refer to
 * files 'e_likelihood.h' and 'e_likelihood.C'. 
 */

#include <cmath>
#include <sys/time.h>

#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorFunction.h>
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/InverseGammaVectorRV.h>
#include <queso/ConcatenationSubset.h>
#include <queso/ConcatenatedVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>

#include <compute.h>
#include <likelihood.h>


void computeE(const QUESO::FullEnvironment& env) 
{
  struct timeval timevalNow;

  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "\nBeginning run of 'Steel Model' at "
              << ctime(&timevalNow.tv_sec)
              << "\n my fullRank = "         << env.fullRank()
              << "\n my subEnvironmentId = " << env.subId()
              << "\n my subRank = "          << env.subRank()
              << "\n my interRank = "        << env.inter0Rank()
               << std::endl << std::endl;
  }

  // Just examples of possible calls
  if ((env.subDisplayFile()       ) &&
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Beginning run of 'Steel Model' at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  env.fullComm().Barrier();
  env.subComm().Barrier();  // Just an example of a possible call

  //================================================================
  // Statistical inverse problem (SIP): find posterior PDF for 'E'
  //================================================================
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'SIP -> Modulus estimation' at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  std::vector<std::string> paramNamesA(4,"");
  paramNamesA [0] = "youngModulus";
  QUESO::VectorSpace<> paramSpaceA(env, "paramA_", paramNamesA.size(),\
				  &paramNamesA);
  
  QUESO::VectorSpace<> paramSpaceB(env, "paramB_",3,NULL);
  QUESO::VectorSpace<> paramSpace(env, "param_", 7,NULL);

  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinA(paramSpaceA.zeroVector());
  QUESO::GslVector paramMaxA(paramSpaceA.zeroVector());
  paramMinA[0] = 1000;
  paramMaxA[0] = 300000;
  paramMinA[1] = .0000001;
  paramMaxA[1] = .0001;
  paramMinA[2] = .00001;
  paramMaxA[2] = .001;
  paramMinA[3] = .0001;
  paramMaxA[3] = .01;
  QUESO::BoxSubset<> paramDomainA("paramA_", paramSpaceA, paramMinA,
      paramMaxA);

  QUESO::GslVector paramMinB(paramSpaceB.zeroVector());
  QUESO::GslVector paramMaxB(paramSpaceB.zeroVector());

  paramMinB[0] = .001;
  paramMaxB[0] = .1;
  paramMinB[1] = .01;
  paramMaxB[1] = .1;
  paramMinB[2] = .1;
  paramMaxB[2] = .3;
  QUESO::BoxSubset<> paramDomainB("paramB_", paramSpaceB, paramMinB,
      paramMaxB);

  QUESO::ConcatenationSubset<>paramDomain("param_",
      paramSpace,paramDomainA,paramDomainB);

  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function
  // object to be used by QUESO.
  //------------------------------------------------------
  Likelihood<> lhood("like_", paramDomain);

  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  QUESO::UniformVectorRV<> priorRvA("priorA_", paramDomainA);

  QUESO::GslVector alpha(paramSpaceB.zeroVector());
  alpha[0] = 20;
  alpha[1] = 20;
  alpha[2] = 20;

  QUESO::GslVector beta(paramSpaceB.zeroVector());
  beta[0] = 1;
  beta[1] = 1;
  beta[2] = 1;
  QUESO::InverseGammaVectorRV<QUESO::GslVector,QUESO::GslMatrix>
      priorRvB("priorB_",paramDomainB,alpha,beta);

  QUESO::ConcatenatedVectorRV<QUESO::GslVector,QUESO::GslMatrix>
      priorRv("prior_",priorRvA,priorRvB,paramDomain);

  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  // Extra prefix before the default "rv_" prefix
  QUESO::GenericVectorRV<> postRv("post_", paramSpace);

  // No extra prefix before the default "ip_" prefix
  QUESO::StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV
  //------------------------------------------------------
  
  ip.solveWithBayesMLSampling();


}
