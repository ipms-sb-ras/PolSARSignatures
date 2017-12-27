/*=========================================================================

Program:
Language:  C++
Date:      $Date: $
Version:   $Revision: $

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

@file    itkLocalSigmaImageFunction.hxx

*=========================================================================*/
#ifndef itkLocalSigmaImageFunction_hxx
#define itkLocalSigmaImageFunction_hxx

#include "itkLocalStatisticImageFunctionBase.hxx"
#include "vnl/vnl_math.h"

namespace itk
{
	/**
	 * \class LocalSigmaImageFunction
	 * \brief Calculate the sigma value in the neighborhood of a pixel
	 *
	 * Calculate the sigma pixel value over the given neighborhood.  
	 * This calculation uses a ZeroFluxNeumannBoundaryCondition.
	 *
	 * If called with a ContinuousIndex or Point, the calculation is performed
	 * at the nearest neighbor.
	 *
	 * This class is templated over the input image type and the
	 * coordinate representation type (e.g. float or double).
	 *
	 * \sa LocalSigmaImageFunction
	 *
	 * \ingroup ImageFunctions
	 * \ingroup ITKImageFunction
	 */
	template< typename TInputImage, typename TCoordRep = float >
	class LocalSigmaImageFunction :
		public LocalStatisticImageFunctionBase < TInputImage, TCoordRep >
	{
	public:
		/** Standard class typedefs. */
		typedef LocalSigmaImageFunction Self;
		typedef LocalStatisticImageFunctionBase< TInputImage, TCoordRep > Baseclass;

                typedef typename Baseclass::RealType RealType;
                typedef typename Baseclass::IndexType IndexType;
                

		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Run-time type information (and related methods). */
		itkTypeMacro(LocalSigmaImageFunction, LocalStatisticImageFunctionBase);

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Evalulate the function at specified index */
		virtual RealType EvaluateAtIndex(const IndexType & index) const ITK_OVERRIDE
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

			ConstNeighborhoodIterator< TInputImage >
				It(this->m_NeighborhoodPointer->GetRadius(), this->GetInputImage(),
				this->GetInputImage()->GetBufferedRegion());
			
			/** "��������" ����������� ������ */
			this->m_NeighborhoodPointer->UseCentralPixelOn();
			
			// Set the iterator at the desired location
			It.SetLocation(index);
			itk::NeighborhoodInnerProduct<TInputImage> innerProduct;
			/** Computing local mean value*/
			RealType localMean = static_cast< RealType >(
				innerProduct(It, *this->m_NeighborhoodPointer) / 
				this->m_NeighborhoodPointer->GetSumOfWeights());
			
			/** Computing local sigma value*/
			size_t itSize = It.Size();
			RealType sum = NumericTraits<RealType>::ZeroValue();
			for (size_t i = 0; i < itSize; ++i)
			{
				sum += vnl_math_sqr(It.GetPixel(i) * 
                                        (*this->m_NeighborhoodPointer)[i]);
			}
			RealType sigma = std::sqrt(sum / 
                                this->m_NeighborhoodPointer->GetSumOfWeights() - vnl_math_sqr(localMean)
                                );
			return (sigma);
		}

	protected:
		LocalSigmaImageFunction() : Baseclass(){};
		~LocalSigmaImageFunction(){};

	private:
		LocalSigmaImageFunction(const Self &); //purposely not implemented
		void operator=(const Self &);    //purposely not implemented
	};
} // end namespace itk

#endif //itkLocalSigmaImageFunction_hxx
