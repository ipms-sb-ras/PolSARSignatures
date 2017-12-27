/*=========================================================================

Program:
Language:  C++
Date:      $Date: $
Version:   $Revision: $

    Copyright (c) IPMS SB RAS

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

@file    itkLocalMeanImageFunction.hxx

*=========================================================================*/
#ifndef itkLocalMeanImageFunction_hxx
#define itkLocalMeanImageFunction_hxx

#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include "itkLocalStatisticImageFunctionBase.hxx"

namespace itk
{
	/**
	 * \class LocalMeanImageFunction
	 * \brief Calculate the mean value in the neighborhood of a pixel
	 *
	 * Calculate the mean pixel value over the given neighborhood.  
	 * This calculation uses a ZeroFluxNeumannBoundaryCondition.
	 *
	 * If called with a ContinuousIndex or Point, the calculation is performed
	 * at the nearest neighbor.
	 *
	 * This class is templated over the input image type and the
	 * coordinate representation type (e.g. float or double).
	 *
	 * \sa LocalMeanImageFunction
	 *
	 * \ingroup ImageFunctions
	 * \ingroup ITKImageFunction
	 */
	template< typename TInputImage, typename TCoordRep = float >
	class LocalMeanImageFunction :
		public LocalStatisticImageFunctionBase < TInputImage, TCoordRep >
	{
	public:
		/** Standard class typedefs. */
		typedef LocalMeanImageFunction Self;
		typedef LocalStatisticImageFunctionBase< TInputImage, TCoordRep > Baseclass;
                
                typedef typename Baseclass::RealType RealType;
                typedef typename Baseclass::IndexType IndexType;

		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Run-time type information (and related methods). */
		itkTypeMacro(LocalMeanImageFunction, LocalStatisticImageFunctionBase);

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

			itk::ConstNeighborhoodIterator< TInputImage >
				It(this->m_NeighborhoodPointer->GetRadius(), this->GetInputImage(),
				this->GetInputImage()->GetBufferedRegion());
			
			/** "Turn on" central pixel */
			this->m_NeighborhoodPointer->UseCentralPixelOn();
			
			// Set the iterator at the desired location
			It.SetLocation(index);
			itk::NeighborhoodInnerProduct<TInputImage> innerProduct;
			/** Computing local mean value*/
			RealType localMean = static_cast< RealType >(
				innerProduct(It, *(this->m_NeighborhoodPointer)) / 
				this->m_NeighborhoodPointer->GetSumOfWeights());

			return (localMean);
		}

	protected:
		LocalMeanImageFunction() : Baseclass(){};
		~LocalMeanImageFunction(){};

	private:
		LocalMeanImageFunction(const Self &); //purposely not implemented
		void operator=(const Self &);    //purposely not implemented
	};
} // end namespace itk

#endif //itkLocalMeanImageFunction_hxx
