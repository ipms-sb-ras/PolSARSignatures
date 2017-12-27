/*
 * otbPolarimetricSynthesisStokesImageFunction.txx
 *
 *  Created on: Jun 7, 2016
 *      Author: dav
 */

#ifndef __otbPolarimetricSynthesisStokesImageFunction_txx
#define __otbPolarimetricSynthesisStokesImageFunction_txx

#include "otbPolarimetricSynthesisStokesImageFunction.h"
#include "itkLocalStatisticImageFunctionBase.hxx"
#include "itkConstNeighborhoodIterator.h"

#include <functional>


namespace otb {

/**
 * Force Copolar mode
 */
template <typename TInputImage, typename TCoordRep>
void
PolarimetricSynthesisStokesImageFunction<TInputImage, TCoordRep>
::ForceCoPolar() {
    SetPsiR(m_PsiI);
    SetKhiR(m_KhiI);
    SetMode(1);
}

/**
 * Force Crosspolar mode
 */
template <typename TInputImage, typename TCoordRep>
void
PolarimetricSynthesisStokesImageFunction<TInputImage, TCoordRep>
::ForceCrossPolar() {
    SetPsiR(m_PsiI + 90);
    SetKhiR(-m_KhiI);
    SetMode(2);
}

template <typename TInputImage, typename TCoordRep>
RealType
PolarimetricSynthesisStokesImageFunction<TInputImage, TCoordRep>
::EvaluateAtIndex(const IndexType & index) const
 {
	if (!this->GetInputImage())
	{
		return (NumericTraits< RealType >::max());
	}

	if (!this->IsInsideBuffer(index))
	{
		return (NumericTraits< RealType >::max());
	}

	if (!this->m_NeighborhoodPointer)
	{
		itkExceptionMacro(<< "Neighborhood is not set!");
	}

	this->m_NeighborhoodPointer->UseCentralPixelOn();

	ConstNeighborhoodIterator< TInputImage >
		It(this->m_NeighborhoodPointer->GetRadius(), this->GetInputImage(),
		this->GetInputImage()->GetBufferedRegion());

	// Set the iterator at the desired location
	It.SetLocation(index);

	// Averaging backscattering values
	typedef typename itk::Vector<ComplexType, 4>     ComplexVectorType;
	ComplexVectorType zeroVector;
	zeroVector.Fill(ComplexType(0.0,0.0));

	ComplexVectorType averagedVector = std::accumulate(
		this->m_NeighborhoodPointer->Begin(), this->m_NeighborhoodPointer->End(),
		zeroVector, std::plus<ComplexVectorType>() )
	/ static_cast<double>(this->m_NeighborhoodPointer->GetSumOfWeights());

	// Calculating polarimetric response

	double polResponse = this->CalculatePolarimetricResponse(averagedVector[0],
			averagedVector[1], averagedVector[2], averagedVector[3]);

	return polResponse;
 }



/**
 * Printself
 */
template <typename TInputImage, typename TCoordRep>
void
PolarimetricSynthesisStokesImageFunction<TInputImage, TCoordRep>
::PrintSelf(std::ostream& os, itk::Indent indent) const {
    this->Baseclass::PrintSelf(os, indent);
    os << indent << "PsiI: " << m_PsiI << std::endl;
    os << indent << "KhiI: " << m_KhiI << std::endl;
    os << indent << "PsiR: " << m_PsiR << std::endl;
    os << indent << "KhiR: " << m_KhiR << std::endl;
    os << indent << "Gain: " << m_Gain << std::endl;
}

template <typename TInputImage, typename TCoordRep>
inline double
PolarimetricSynthesisStokesImageFunction<TInputImage, TCoordRep>
::CalculatePolarimetricResponse(const InputImagePixelType Shh,
        const InputImagePixelType Shv, const InputImagePixelType Svh,
        const InputImagePixelType Svv) {

    typedef typename itk::Matrix<double, 4, 4>  StokesMatrixType;
    typedef typename itk::Vector<double, 4>     StokesVectorType;

    const ComplexType I(0.0, 1.0); // imagery unit
    const ComplexType S_hh = static_cast<ComplexType> (Shh);
    const ComplexType S_hv = static_cast<ComplexType> (Shv);
    const ComplexType S_vh = static_cast<ComplexType> (Svh);
    const ComplexType S_vv = static_cast<ComplexType> (Svv);

    const ComplexType S_hh_conj = vcl_conj(S_hh);
    const ComplexType S_hv_conj = vcl_conj(S_hv);
    const ComplexType S_vh_conj = vcl_conj(S_vh);
    const ComplexType S_vv_conj = vcl_conj(S_vv);

    const double S_hh_norm = std::norm(S_hh);
    const double S_hv_norm = std::norm(S_hv);
    const double S_vh_norm = std::norm(S_vh);
    const double S_vv_norm = std::norm(S_vv);

    // Define the Stokes scattering operator
    StokesMatrixType StokesMatrix;
    StokesMatrix[0][0] = 0.25 * (S_hh_norm + S_hv_norm + S_vh_norm + S_vv_norm);
    StokesMatrix[0][1] = 0.25 * (S_hh_norm - S_hv_norm + S_vh_norm - S_vv_norm);
    StokesMatrix[0][2] = 0.25 * std::real(S_hh * S_hv_conj + S_hv * S_hh_conj + S_vh * S_vv_conj + S_vv * S_vh_conj);
    StokesMatrix[0][3] = 0.25 * std::real(I * (S_hh * S_hv_conj - S_hv * S_hh_conj + S_vh * S_vv_conj - S_vv * S_vh_conj));
    StokesMatrix[1][0] = 0.25 * (S_hh_norm + S_hv_norm - S_vh_norm - S_vv_norm);
    StokesMatrix[1][1] = 0.25 * (S_hh_norm - S_hv_norm - S_vh_norm + S_vv_norm);
    StokesMatrix[1][2] = 0.25 * std::real(S_hh * S_hv_conj + S_hv * S_hh_conj - S_vh * S_vv_conj - S_vv * S_vh_conj);
    StokesMatrix[1][3] = 0.25 * std::real(I * (S_hh * S_hv_conj - S_hv * S_hh_conj - S_vh * S_vv_conj + S_vv * S_vh_conj));
    StokesMatrix[2][0] = 0.25 * std::real(S_hh * S_vh_conj + S_hv * S_vv_conj + S_vh * S_hh_conj + S_vv * S_hv_conj);
    StokesMatrix[2][1] = 0.25 * std::real(S_hh * S_vh_conj - S_hv * S_vv_conj + S_vh * S_hh_conj - S_vv * S_hv_conj);
    StokesMatrix[2][2] = 0.25 * std::real(S_hh * S_vv_conj + S_hv * S_vh_conj + S_vh * S_hv_conj + S_vv * S_hh_conj);
    StokesMatrix[2][3] = 0.25 * std::real(I * (S_hh * S_vv_conj - S_hv * S_vh_conj + S_vh * S_hv_conj - S_vv * S_hh_conj));
    StokesMatrix[3][0] = 0.25 * std::real(I * (S_hh * S_vh_conj + S_hv * S_vv_conj - S_vh * S_hh_conj - S_vv * S_hv_conj));
    StokesMatrix[3][1] = 0.25 * std::real(I * (S_hh * S_vh_conj - S_hv * S_vv_conj - S_vh * S_hh_conj + S_vv * S_hv_conj));
    StokesMatrix[3][2] = 0.25 * std::real(I * (S_hh * S_vv_conj + S_hv * S_vh_conj - S_vh * S_hv_conj - S_vv * S_hh_conj));
    StokesMatrix[3][3] = 0.25 * std::real(-S_hh * S_vv_conj + S_hv * S_vh_conj + S_vh * S_hv_conj - S_vv * S_hh_conj);

    const double const_PI_90 = 2. * CONST_PI_180;
    // Define the incident Stokes vector
    StokesVectorType Si;
    Si[0] = 1.0;
    Si[1] = cos(m_PsiI * const_PI_90) * cos(m_KhiI * const_PI_90);
    Si[2] = sin(m_PsiI * const_PI_90) * cos(m_KhiI * const_PI_90);
    Si[3] = sin(m_KhiI * const_PI_90);

    // Evaluate the product (Stokes scattering operator x Stokes incidence vector)
    StokesVectorType Smi = StokesMatrix * Si;

    // Define the reflected Stokes vector
    StokesVectorType Sr;
    Sr[0] = 1.0;
    Sr[1] = cos(m_PsiR * const_PI_90) * cos(m_KhiR * const_PI_90);
    Sr[2] = sin(m_PsiR * const_PI_90) * cos(m_KhiR * const_PI_90);
    Sr[3] = sin(m_KhiR * const_PI_90);

    double sigma = /* 4.* CONST_PI * */
        (Sr[0] * Smi[0] + Sr[1] * Smi[1] + Sr[2] * Smi[2] + Sr[3] * Smi[3]);

    return (sigma);
}

}  // namespace otb


#endif /* __otbPolarimetricSynthesisStokesImageFunction_txx */
