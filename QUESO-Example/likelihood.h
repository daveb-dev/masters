#ifndef QUESO_EXAMPLE_GRAVITY_LIKELIHOOD_H
#define QUESO_EXAMPLE_GRAVITY_LIKELIHOOD_H

#include <queso/ScalarFunction.h>
#include <queso/GslMatrix.h>

template<class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:
  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain);
  virtual ~Likelihood();
  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const;
  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const;

private:
  std::vector<double> stress; // stress
  std::vector<double> strains;   // strains
  std::vector<double> stdevs; // standard deviations

};

#endif
