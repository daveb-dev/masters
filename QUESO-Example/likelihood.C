#include <cmath>
#include <iostream>
#include <fstream> // Stream class to read from files
#include <sstream> // To parse strings
#include <limits> // for NaN
#include <string>

#include <queso/GenericScalarFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <likelihood.h>


template<class V, class M>
Likelihood<V, M>::Likelihood(const char * prefix,
    const QUESO::VectorSet<V, M> & domain)
  : QUESO::BaseScalarFunction<V, M>(prefix, domain),
    stress(0),strains(0),stdevs(0)
{
    std::ifstream data ("steeldata.csv");

    // Determine the number of tests (columns)
    std::string line, cell;
    int cols = 0;
    int maxstress = 600; ///////// SET THIS FOR DIFFERENT REGIMES
    if (data.is_open())
    {
	getline(data,line);
	std::stringstream lineStream(line);
	while(getline(lineStream,cell,','))
	{
	    cols += 1;
	}
    }
    else{std::cout << "Data not found" << std::endl;}

    data.clear(); //Clear warnings
    data.seekg(0,data.beg); //Start at beginning of data
 
    // Put stress/strain values into vectors and determine std dev
    int i = 0;
    int var1;
    double var2;
    while(getline(data,line))
    {
	std::stringstream lineStream(line);
	int j = 0;
	int k = 0;
	double mean = 0;
	double stdvar = 0;
	while(getline(lineStream,cell,','))
	{
	    std::istringstream s(cell);
	    if((i % cols) == 0) //First column is stress value
	    {
		s >> var1;
		if(var1 > maxstress)
		{
		    break;
		}
		stress.push_back(var1);
	    }
	    else if(cell.empty()) //Cell is empty -> NaN
	    {
		strains.push_back(std::numeric_limits<double>::quiet_NaN());
	    }
	    else //Remaining columns are strain
	    {
		s >> var2;
		strains.push_back(var2);
		mean += var2;
		k += 1;
	    }
	    j += 1;
	    i += 1; 
      	}
	if(var1 > maxstress) //Stop for regime
	{
	    break;
	}
	if(j<cols) //Empty column at end
 	{
	    strains.push_back(std::numeric_limits<double>::quiet_NaN());
	    i += 1;
	}
	//Calculate std dev
	mean = mean/k;
	for(int l=0; l<k; l++){
	    if(isnan(strains.rbegin()[l])){continue;}
	    else{stdvar += pow(strains.rbegin()[l]-mean,2);}
	}
	stdvar = pow(stdvar/(k-1),.5);
	if(stdvar<1e-6){stdvar = 1e-6;}
	stdevs.push_back(stdvar);
    }
}

template<class V, class M>
Likelihood<V, M>::~Likelihood()
{
  // Deconstruct here
}

template<class V, class M>
double
Likelihood<V, M>::lnValue(const V & domainVector, const V * domainDirection,
    V * gradVector, M * hessianMatrix, V * hessianEffect) const
{
    
  double E = domainVector[0];
  double s1 = domainVector[1];
  double s2 = domainVector[2];
  double s3 = domainVector[3];
  double s4 = domainVector[4];
  double s5 = domainVector[5];
  double s6 = domainVector[6];
  double misfitValue = 0.0, ratio, modelStrain;
  int cols = strains.size()/stress.size();
  for (unsigned int j = 0; j < stress.size(); ++j)
  {
      modelStrain = stress[j]/E;
      for(int k = 0; k < cols; ++k)
      {
	  if(isnan(strains[k+j*cols]))
	  {
	      continue;
	  }
	  else if(stress[j] < 100)
	  {
	      ratio = (modelStrain - strains[k+j*cols])/s1;
	  }
	  else if(stress[j] < 200)
	  {
	      ratio = (modelStrain - strains[k+j*cols])/s2;
	  }
	  else if(stress[j] < 300)
	  {
	      ratio = (modelStrain - strains[k+j*cols])/s3;
	  }
	  else if(stress[j] < 400)
	  {
	      ratio = (modelStrain - strains[k+j*cols])/s4;
	  }
	  else if(stress[j] < 500)
	  {
	      ratio = (modelStrain - strains[k+j*cols])/s5;
	  }
	  else
	  {
	      ratio = (modelStrain - strains[k+j*cols])/s6;
	  }
	  misfitValue += ratio*ratio;
      }
  }
  return -0.5 * misfitValue - 0.5*(34.*log(std::pow(2.*M_PI*s1*s2,2))+
		     +33.*log(std::pow(4.*M_PI*M_PI*s3*s4*s5*s6,2)));
}

template<class V, class M>
double
Likelihood<V, M>::actualValue(const V & domainVector,
    const V * domainDirection, V * gradVector, M * hessianMatrix,
    V * hessianEffect) const
{
  return std::exp(this->lnValue(domainVector, domainDirection, gradVector,
        hessianMatrix, hessianEffect));
}

template class Likelihood<QUESO::GslVector, QUESO::GslMatrix>;
