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
#ifndef otbSinclairToStokesMatrixFunctor_h
#define otbSinclairToStokesMatrixFunctor_h

#include <complex>

namespace otb
{
namespace Functor
{
/** \class SinclairToStokesMatrixFunctor
 *  \brief Construct Stokes scattering operator from Sinclair matrix.
 *  Elements of the Stokes scattering operator (matrix) are extract from 
 *  Remote Sensing with Imaging Radar.  John A. Richards p 97-98, 2009
 *
 *  Output value are:
 *  - channel #0  : \f$ m_{11} = 0.25 * ( S_{HH}S_{HH}^{*} + S_{HV}S_{HV}^{*} + S_{VH}S_{VH}^{*} + S_{VV}S_{VV}^{*} ) \f$
 *  - channel #1  : \f$ m_{12} = 0.25 * ( S_{HH}S_{HH}^{*} - S_{HV}S_{HV}^{*} + S_{VH}S_{VH}^{*} - S_{VV}S_{VV}^{*} ) \f$
 *  - channel #2  : \f$ m_{13} = 0.25 * ( S_{HH}S_{HV}^{*} + S_{HV}S_{HH}^{*} + S_{VH}S_{VV}^{*} + S_{VV}S_{VH}^{*} ) \f$
 *  - channel #3  : \f$ m_{14} = 0.25j * ( S_{HH}S_{HV}^{*} - S_{HV}S_{HH}^{*} + S_{VH}S_{VV}^{*} - S_{VV}S_{VH}^{*} ) \f$
 *  - channel #4  : \f$ m_{21} = 0.25 * ( S_{HH}S_{HH}^{*} + S_{HV}S_{HV}^{*} - S_{VH}S_{VH}^{*} - S_{VV}S_{VV}^{*} ) \f$
 *  - channel #5  : \f$ m_{22} = 0.25 * ( S_{HH}S_{HH}^{*} - S_{HV}S_{HV}^{*} - S_{VH}S_{VH}^{*} + S_{VV}S_{VV}^{*} ) \f$
 *  - channel #6  : \f$ m_{23} = 0.25 * ( S_{HH}S_{HV}^{*} + S_{HV}S_{HH}^{*} - S_{VH}S_{VV}^{*} - S_{VV}S_{VH}^{*} ) \f$
 *  - channel #7  : \f$ m_{24} = 0.25j * ( S_{HH}S_{HV}^{*} - S_{HV}S_{HH}^{*} - S_{VH}S_{VV}^{*} + S_{VV}S_{VH}^{*} ) \f$
 *  - channel #8  : \f$ m_{31} = 0.25 * ( S_{HH}S_{VH}^{*} + S_{HV}S_{VV}^{*} + S_{VH}S_{HH}^{*} + S_{VV}S_{HV}^{*} ) \f$
 *  - channel #9  : \f$ m_{32} = 0.25 * ( S_{HH}S_{VH}^{*} - S_{HV}S_{VV}^{*} + S_{VH}S_{HH}^{*} - S_{VV}S_{HV}^{*} ) \f$
 *  - channel #10 : \f$ m_{33} = 0.25 * ( S_{HH}S_{VV}^{*} + S_{HV}S_{VH}^{*} + S_{VH}S_{HV}^{*} + S_{VV}S_{HH}^{*} ) \f$
 *  - channel #11 : \f$ m_{34} = 0.25j * ( S_{HH}S_{VV}^{*} - S_{HV}S_{VH}^{*} + S_{VH}S_{HV}^{*} - S_{VV}S_{HH}^{*} ) \f$
 *  - channel #12 : \f$ m_{41} = 0.25j * ( S_{HH}S_{VH}^{*} + S_{HV}S_{VV}^{*} - S_{VH}S_{HH}^{*} - S_{VV}S_{HV}^{*} ) \f$
 *  - channel #13 : \f$ m_{42} = 0.25j * ( S_{HH}S_{VH}^{*} - S_{HV}S_{VV}^{*} - S_{VH}S_{HH}^{*} + S_{VV}S_{HH}^{*} ) \f$
 *  - channel #14 : \f$ m_{43} = 0.25j * ( S_{HH}S_{VV}^{*} + S_{HV}S_{VH}^{*} - S_{VH}S_{HV}^{*} - S_{VV}S_{HH}^{*} ) \f$
 *  - channel #15 : \f$ m_{44} = 0.25 * ( -S_{HH}S_{VV}^{*} + S_{HV}S_{VH}^{*} + S_{VH}S_{HV}^{*} - S_{VV}S_{HH}^{*} ) \f$
 *
 *
 * Output is a not a complex. The output pixel has 16 channels : each element of the Stokes matrix.
 * The order of the channels corresponds to :
 * \f$  \begin{pmatrix}
 * {channel #0 }&{channel #1 }&{channel #2 }&{channel #3 } \\
 * {channel #4 }&{channel #5 }&{channel #6 }&{channel #7 } \\
 * {channel #8 }&{channel #9 }&{channel #10}&{channel #11} \\
 * {channel #12}&{channel #13}&{channel #14}&{channel #15} \\
 * \end{pmatrix}  \f$
 *
 *  \ingroup Functor
 *  \ingroup SARPolarimetry
 *
 *  \sa SinclairImageFilter
 *  \sa SinclairToCircularCovarianceMatrixFunctor
 *  \sa SinclairToCoherencyMatrixFunctor
 *  \sa SinclairToCovarianceMatrixFunctor
 *  \sa SinclairToReciprocalCircularCovarianceMatrixFunctor
 *  \sa SinclairToReciprocalCoherencyMatrixFunctor
 *  \sa SinclairToReciprocalCovarianceMatrixFunctor
 *
 *
 * \ingroup OTBPolarimetry
 */
template <class TInput1, class TInput2, class TInput3,
          class TInput4, class TOutput>
class SinclairToStokesMatrixFunctor
{
public:
  /** Some typedefs. */
  typedef typename std::complex <double>           ComplexType;
  typedef typename TOutput::ValueType              OutputValueType;
  typedef double                                   RealType;

  inline TOutput operator ()(const TInput1& Shh, const TInput2& Shv,
                             const TInput3& Svh, const TInput4& Svv)
  {
    TOutput result;

    result.SetSize(m_NumberOfComponentsPerPixel);

    const ComplexType Txx = static_cast<ComplexType>(Shh);
    const ComplexType Txy = static_cast<ComplexType>(Shv);
    const ComplexType Tyx = static_cast<ComplexType>(Svh);
    const ComplexType Tyy = static_cast<ComplexType>(Svv);

    const ComplexType conjTxx = std::conj(Txx);
    const ComplexType conjTxy = std::conj(Txy);
    const ComplexType conjTyx = std::conj(Tyx);
    const ComplexType conjTyy = std::conj(Tyy);
    
    const RealType  normTxx = std::norm(Txx);
    const RealType  normTxy = std::norm(Txy);
    const RealType  normTyx = std::norm(Tyx);
    const RealType  normTyy = std::norm(Tyy);

    result[0]  = static_cast<OutputValueType>( 0.25 * (normTxx + normTxy + normTyx + normTyy) );
    result[1]  = static_cast<OutputValueType>( 0.25 * (normTxx - normTxy + normTyx - normTyy) );
    result[2]  = static_cast<OutputValueType>( 0.5  * (Txx*conjTxy + Tyx*conjTyy).real()  );
    result[3]  = static_cast<OutputValueType>( 0.5  * (Txy*conjTxx + Tyy*conjTyx).imag()  );
    result[4]  = static_cast<OutputValueType>( 0.25 * (normTxx + normTxy - normTyx - normTyy) );
    result[5]  = static_cast<OutputValueType>( 0.25 * (normTxx - normTxy - normTyx + normTyy) );
    result[6]  = static_cast<OutputValueType>( 0.5  * (Txx*conjTxy - Tyx*conjTyy).real()  );
    result[7]  = static_cast<OutputValueType>( 0.5  * (Txy*conjTxx + Tyx*conjTyy).imag()  );
    result[8]  = static_cast<OutputValueType>( 0.5  * (Txx*conjTyx + Txy*conjTyy).real()  );
    result[9]  = static_cast<OutputValueType>( 0.5  * (Txx*conjTyx - Txy*conjTyy).real()  );
    result[10] = static_cast<OutputValueType>( 0.5  * (Txx*conjTyy + Txy*conjTyx).real()  );
    result[11] = static_cast<OutputValueType>( 0.5  * (Tyy*conjTxx + Txy*conjTyx).imag()  );
    result[12] = static_cast<OutputValueType>( 0.5  * (Tyx*conjTxx + Tyy*conjTxy).imag()  );
    result[13] = static_cast<OutputValueType>( 0.5  * (Tyx*conjTxx + Txy*conjTyy).imag()  );
    result[14] = static_cast<OutputValueType>( 0.5  * (Tyy*conjTxx + Tyx*conjTxy).imag()  );
    result[15] = static_cast<OutputValueType>( 0.5  * (Txy*conjTyx - Txx*conjTyy).real()  );

    return (result);
  }

  unsigned int GetNumberOfComponentsPerPixel()
  {
    return m_NumberOfComponentsPerPixel;
  }

  /** Constructor */
  SinclairToStokesMatrixFunctor() : m_NumberOfComponentsPerPixel(16) {}

  /** Destructor */
  virtual ~SinclairToStokesMatrixFunctor() {}

protected:


private:
    unsigned int m_NumberOfComponentsPerPixel;
};

} // namespace Functor
} // namespace otb

#endif //otbSinclairToStokesMatrixFunctor_h
