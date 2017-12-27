/*=========================================================================

Program:
Language:  C++
Date:      $Date: $
Version:   $Revision: $

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

@file    LocalStatisticImageFunctionBase.hxx

*=========================================================================*/
#ifndef itkLocalStatisticImageFunctionBase_hxx
#define itkLocalStatisticImageFunctionBase_hxx

#include <memory>

#include "itkImageFunction.h"
#include "itkLocalStatisticNeighborhoodBase.hxx"
#include "itkNumericTraits.h"

namespace itk
{
	/**
	* \class LocalStatisticImageFunctionBase
	* \brief A base class for local statistic image functions
	*
	* \sa LocalStatisticImageFunctionBase
	*
	* \ingroup ImageFunctions
	* \ingroup ITKImageFunction
	*/
	template< typename TInputImage, typename TCoordRep = float >
	class LocalStatisticImageFunctionBase :
		public ImageFunction < TInputImage,
		typename NumericTraits< typename TInputImage::PixelType >::RealType,
		TCoordRep >
	{
	public:
		/** Standard class typedefs. */
		typedef LocalStatisticImageFunctionBase Self;
		typedef ImageFunction< TInputImage,
			typename NumericTraits< typename TInputImage::PixelType >::RealType,
			TCoordRep >                     Superclass;

		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Run-time type information (and related methods). */
		itkTypeMacro(LocalStatisticImageFunctionBase, ImageFunction);

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** InputImageType typedef support. */
		typedef TInputImage InputImageType;

		/** OutputType typdef support. */
		typedef typename Superclass::OutputType OutputType;

		/** Index typedef support. */
		typedef typename Superclass::IndexType IndexType;

		/** ContinuousIndex typedef support. */
		typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

		/** Point typedef support. */
		typedef typename Superclass::PointType PointType;

		/** Datatype used */
		typedef typename NumericTraits< typename InputImageType::PixelType >::RealType
			RealType;

		typedef LocalStatisticNeighborhoodBase<typename InputImageType::PixelType>	NeighborhoodType;
		typedef std::shared_ptr<NeighborhoodType>				NeighborhoodPointerType;

		/** Evalulate the function at specified index */
		virtual RealType EvaluateAtIndex(const IndexType & index) const ITK_OVERRIDE 
		{
			itkExceptionMacro(<< "You must override this function in instantiated class!");
			return (NumericTraits< RealType >::max()); 
		}

		/** Evaluate the function at non-integer positions */
		virtual RealType Evaluate(const PointType & point) const ITK_OVERRIDE
		{
			IndexType index;
			this->ConvertPointToNearestIndex(point, index);
			return this->EvaluateAtIndex(index);
		}

		virtual RealType EvaluateAtContinuousIndex(
		const ContinuousIndexType & cindex) const ITK_OVERRIDE
		{
			IndexType index;
			this->ConvertContinuousIndexToNearestIndex(cindex, index);
			return this->EvaluateAtIndex(index);
		}

		/**
		* Overriding SetInputImage method in order to compute image statistics
		*/
		virtual void SetInputImage(const TInputImage *ptr)
		{
			Superclass::SetInputImage(ptr);
		}

		/** Setting neighborhood. */
		template <typename TNeighborhood>
		void SetNeighborhood(const SizeValueType s)
		{
			m_NeighborhoodPointer = std::make_shared<TNeighborhood>();
			/** "Turn off" central pixel by default */
			m_NeighborhoodPointer->UseCentralPixelOff();
			m_NeighborhoodPointer->SetRadius(s);
			this->Modified();
		}
		/** Setting neighborhood. */
		void SetNeighborhood(NeighborhoodPointerType _arg)
		{
			if (this->m_NeighborhoodPointer != _arg){
				this->m_NeighborhoodPointer = _arg;
				this->m_NeighborhoodPointer->UseCentralPixelOff();
				this->Modified();
			}
		}

		itkGetConstMacro(NeighborhoodPointer, NeighborhoodPointerType);

		void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE
		{
			this->Superclass::PrintSelf(os, indent);
			os << indent << "Neighborhood: " << *m_NeighborhoodPointer << std::endl;
		}

	protected:
		LocalStatisticImageFunctionBase() :
			Superclass(),
			m_NeighborhoodPointer(nullptr){};

		virtual ~LocalStatisticImageFunctionBase(){};
		
		NeighborhoodPointerType    m_NeighborhoodPointer;

	private:
		LocalStatisticImageFunctionBase(const Self &); //purposely not implemented
		void operator=(const Self &);    //purposely not implemented
	};
} // end namespace itk

#endif //itkLocalStatisticImageFunctionBase_hxx
