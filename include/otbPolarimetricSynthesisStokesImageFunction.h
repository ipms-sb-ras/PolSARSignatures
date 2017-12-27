/*
 * otbPolarimetricSynthesisStokesImageFunction.h
 *
 *  Created on: Jun 7, 2016
 *      Author: dav
 */

#ifndef __otbPolarimetricSynthesisStokesImageFunction_h
#define __otbPolarimetricSynthesisStokesImageFunction_h

#include "itkLocalStatisticImageFunctionBase.hxx"
#include "itkMacro.h"
#include <complex>

namespace otb {

/**
 * \class
 * \brief Calculate the mean value in the neighborhood of a pixel
 *
 * Calculate the mean pixel value over the standard 8, 26, etc. connected
 * neighborhood.  This calculation uses a ZeroFluxNeumannBoundaryCondition.
 *
 * If called with a ContinuousIndex or Point, the calculation is performed
 * at the nearest neighbor.
 *
 * This class is templated over the input image type and the
 * coordinate representation type (e.g. float or double).
 *
 * \sa VectorLocalGetisOrdImageFunction
 *
 * \ingroup ImageFunctions
 * \ingroup ITKImageFunction
 */
template<typename TInputImage, typename TCoordRep = float>
class PolarimetricSynthesisStokesImageFunction: public itk::LocalStatisticImageFunctionBase<
		TInputImage, TCoordRep> {
public:
	/** Standard class typedefs. */
	typedef PolarimetricSynthesisStokesImageFunction Self;
	typedef itk::LocalStatisticImageFunctionBase<TInputImage, TCoordRep> Baseclass;

	typedef typename Baseclass::RealType RealType;
	typedef typename Baseclass::IndexType IndexType;

	typedef SmartPointer<Self> Pointer;
	typedef SmartPointer<const Self> ConstPointer;

	typedef typename TInputImage::InternalPixelType InputImagePixelType;
	typedef typename std::complex <double>          ComplexType;

	/** Run-time type information (and related methods). */
	itkTypeMacro(PolarimetricSynthesisStokesImageFunction, LocalStatisticImageFunctionBase);

	/** Method for creation through the object factory. */
	itkNewMacro(Self)	;


	/** Set/Get PsiI */
    itkSetMacro(PsiI, double);
    itkGetMacro(PsiI, double);
    /** Set/Get KhiI */
    itkSetMacro(KhiI, double);
    itkGetMacro(KhiI, double);

    /** Set/Get PsiR */
    itkSetMacro(PsiR, double);
    itkGetMacro(PsiR, double);
    /** Set/Get KhiR */
    itkSetMacro(KhiR, double);
    itkGetMacro(KhiR, double);

    /** Set/Get Mode */
    itkSetMacro(Mode, int);
    itkGetMacro(Mode, int);

    /** Set the gain */
    itkSetMacro(Gain, double);
    itkGetMacro(Gain, double);

	/** Evalulate the function at specified index */
	virtual RealType EvaluateAtIndex(const IndexType & index) const ITK_OVERRIDE;

	void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;


    /** Force the copolar mode */
    void ForceCoPolar();
    /** Force the crosspolar mode */
    void ForceCrossPolar();


protected:
	PolarimetricSynthesisStokesImageFunction() :
			Baseclass(),
			m_PsiI(.0), m_KhiI(.0), m_PsiR(.0), m_KhiR(.0),
			m_Gain(1.0), m_Mode(0) {};

	~PolarimetricSynthesisStokesImageFunction() {}

    /** */
    inline double CalculatePolarimetricResponse(const InputImagePixelType Shh,
            const InputImagePixelType Shv, const InputImagePixelType Svh,
            const InputImagePixelType Svv);

private:
	PolarimetricSynthesisStokesImageFunction(const Self &); //purposely not implemented
	void operator=(const Self &);    //purposely not implemented

    /** Psi Incident */
    double m_PsiI;
    /** Khi Incident */
    double m_KhiI;
    /** Psi Reflected */
    double m_PsiR;
    /** Khi Reflected */
    double m_KhiR;

    /** Gain */
    double m_Gain;

    /** None = 0 , copolar = 1 , crosspolar = 2 */
    int m_Mode;

};

}  // namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbPolarimetricSynthesisStokesImageFunction.txx"
#endif

#endif /* otbPolarimetricSynthesisStokesImageFunction_h */
