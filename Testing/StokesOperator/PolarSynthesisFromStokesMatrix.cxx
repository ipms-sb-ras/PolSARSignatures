/*=========================================================================

  Program:   
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) IPMS SB RAS


     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.

=========================================================================*/

#include "itkMacro.h"

#include "otbSinclairToStokesMatrixFunctor.h"
#include "itkVariableLengthVector.h"
#include "otbMath.h"

int main(int itkNotUsed(argc), char * itkNotUsed(argv)[])
{
  typedef std::complex<double>                   ComplexType;
  typedef itk::VariableLengthVector<double>  OutputType;

  typedef otb::Functor::SinclairToStokesMatrixFunctor<ComplexType, ComplexType, ComplexType, ComplexType, OutputType >         FunctorType;

//  const double const_PI_90 = 2. * otb::CONST_PI_180;
  const double psiI = otb::CONST_PI/6.;
  const double khiI = otb::CONST_PI/12.;
  // Define the incident Stokes vector
  double si0 = 1.0;
  double si1 = std::cos(2.*psiI) * std::cos(2.*khiI);
  double si2 = std::sin(2.*psiI) * std::cos(2.*khiI);
  double si3 = std::sin(2.*khiI);

  const double psiR = otb::CONST_PI/8.; //otb::CONST_PI/6.;
  const double khiR = otb::CONST_PI/8.; //otb::CONST_PI/12.;
  // Define the recieved Stokes vector
  double sr0 = 1.0;
  double sr1 = std::cos(2.*psiR) * std::cos(2.*khiR);
  double sr2 = std::sin(2.*psiR) * std::cos(2.*khiR);
  double sr3 = std::sin(2.*khiR);
  
  FunctorType funct;
  OutputType m = funct.operator ()( ComplexType(1., 4.), ComplexType(2., 3.), ComplexType(3., 2.), ComplexType(4., 1.) );
  
  // Calculating polarimetric response
  long double result = sr0 * ( m[0]*si0 + m[1]*si1 + m[2]*si2 + m[3]*si3);
  result += sr1 * ( m[4]*si0 + m[5]*si1 + m[6]*si2 + m[7]*si3);
  result += sr2 * ( m[8]*si0 + m[9]*si1 + m[10]*si2 + m[11]*si3);
  result += sr3 * ( m[12]*si0 + m[13]*si1 + m[14]*si2 + m[15]*si3);
  
  //const double exactResult = 217./8.;
  const double exactResult = 19.7657508123141;
  
  std::cout.precision(15);
  std::cout<<result<<std::endl;
  std::cout<<vcl_abs(exactResult-result)<<std::endl;
  
  if( vcl_abs(exactResult-result) > 1e-10  )
  {
    std::cout<<"Test gives :"<<std::endl;
    std::cout<<result<<std::endl;
    std::cout<<"Wanted results are :"<<std::endl;
    std::cout<<exactResult<<std::endl;

    return EXIT_FAILURE;
  }
  std::cout<<"Test passed."<<std::endl;
  return EXIT_SUCCESS;
}
