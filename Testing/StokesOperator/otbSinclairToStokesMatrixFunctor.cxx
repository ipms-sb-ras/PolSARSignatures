/*=========================================================================

  Program:   
  Language:  C++
  Date:      $Date$
  Version:   $Revision$



     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  

=========================================================================*/
#include "itkMacro.h"

#include "otbSinclairToStokesMatrixFunctor.h"
#include "itkVariableLengthVector.h"

int main(int itkNotUsed(argc), char * itkNotUsed(argv)[])
{
  typedef std::complex<double>                   ComplexType;
  typedef itk::VariableLengthVector<double>  OutputType;

  typedef otb::Functor::SinclairToStokesMatrixFunctor<ComplexType, ComplexType, ComplexType, ComplexType, OutputType >         FunctorType;

  OutputType  result(16);
  FunctorType funct;
  OutputType outputFunct;

  result[0]  =  15.;
  result[1]  =   0.;
  result[2]  =  14.;
  result[3]  =  -5.;
  result[4]  =   0.;
  result[5]  =   2.;
  result[6]  =   0.;
  result[7]  =   0.;
  result[8]  =  11.;
  result[9]  =   0.;
  result[10] =  10.;
  result[11] =  -5.;
  result[12] = -10.;
  result[13] =   0.;
  result[14] = -10.;
  result[15] =   2.;

  outputFunct = funct.operator ()( ComplexType(1., 4.), ComplexType(2., 3.), ComplexType(3., 2.), ComplexType(4., 1.) );

  if( vcl_abs(result[0]-outputFunct[0]) > 1e-10 ||
      vcl_abs(result[1]-outputFunct[1]) > 1e-10 ||
      vcl_abs(result[2]-outputFunct[2]) > 1e-10 ||
      vcl_abs(result[3]-outputFunct[3]) > 1e-10 ||
      vcl_abs(result[4]-outputFunct[4]) > 1e-10 ||
      vcl_abs(result[5]-outputFunct[5]) > 1e-10 ||
      vcl_abs(result[6]-outputFunct[6]) > 1e-10 ||
      vcl_abs(result[7]-outputFunct[7]) > 1e-10 ||
      vcl_abs(result[8]-outputFunct[8]) > 1e-10 ||
      vcl_abs(result[9]-outputFunct[9]) > 1e-10 ||
      vcl_abs(result[10]-outputFunct[10]) > 1e-10 ||
      vcl_abs(result[11]-outputFunct[11]) > 1e-10 ||
      vcl_abs(result[12]-outputFunct[12]) > 1e-10 ||
      vcl_abs(result[13]-outputFunct[13]) > 1e-10 ||
      vcl_abs(result[14]-outputFunct[14]) > 1e-10 ||
      vcl_abs(result[15]-outputFunct[15]) > 1e-10    )
  {
    std::cout<<"Test gives :"<<std::endl;
    std::cout<<outputFunct<<std::endl;
    std::cout<<"Wanted results are :"<<std::endl;
    std::cout<<result<<std::endl;

    return EXIT_FAILURE;
  }
  std::cout<<"Test passed."<<std::endl;
  return EXIT_SUCCESS;
}
