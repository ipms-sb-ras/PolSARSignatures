/*=========================================================================

Program:
Language:  C++
Date:      $Date: $
Version:   $Revision: $

   Copyright (c) IPMS SB RAS

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

@file    itkLocalStatisticNeighborhoodBase.hxx

=========================================================================*/
#ifndef itkLocalStatisticNeighborhoodBase_hxx
#define itkLocalStatisticNeighborhoodBase_hxx

#include "itkNeighborhood.h"
#include "itkNeighborhoodAllocator.h"
#include "itkNumericTraits.h"

namespace itk {

	/**
	* \class LocalStatisticNeighborhoodBaseBase
	*
	* \brief A Base class for neighborhoods used in spatial statistics calculation.
	*
	*/
	template< typename TPixel, unsigned int VDimension = 2,
		typename TAllocator = NeighborhoodAllocator< TPixel > >
	class LocalStatisticNeighborhoodBase :
		public itk::Neighborhood < TPixel, VDimension, TAllocator > {
	public:
		/** Standard typedefs */
		typedef LocalStatisticNeighborhoodBase Self;
		typedef Neighborhood< TPixel, VDimension, TAllocator > Superclass;

		typedef typename Superclass::SizeType SizeType;
		typedef typename Superclass::SizeValueType SizeValueType;
		typedef typename Superclass::DimensionValueType DimensionValueType;
                typedef typename Superclass::OffsetType OffsetType;

		itkTypeMacro(LocalStatisticNeighborhoodBase, Neighborhood);

		LocalStatisticNeighborhoodBase(const Self & other) :
			Superclass(other),
			m_UseCentralPixel(other.GetUseCentralPixel()),
			m_SumOfWeights(other.GetSumOfWeights())
		{}

		LocalStatisticNeighborhoodBase() :
			Superclass(),
			m_UseCentralPixel(false),
			m_SumOfWeights(NumericTraits<TPixel>::OneValue())
		{}

		LocalStatisticNeighborhoodBase(const SizeType & r) :
			m_UseCentralPixel(false)
		{
			this->SetRadius(r);
		}

		itkGetConstMacro(UseCentralPixel, bool);
		itkGetConstMacro(SumOfWeights, TPixel);
		itkBooleanMacro(UseCentralPixel);

		virtual void SetRadius(const SizeValueType s)
		{
			SizeType k;
			for (DimensionValueType i = 0; i < VDimension; i++)
			{
				k[i] = s;
			}
			this->SetRadius(k);
		}

		virtual void SetRadius(const SizeType & r)
		{
			Superclass::SetRadius(r);
		}

		/*
		* Setting central pixel value
		*/
		void SetUseCentralPixel(bool _arg) {
			if (m_UseCentralPixel != _arg) {
				m_UseCentralPixel = _arg;
				if (m_UseCentralPixel) {
					this->operator[](this->GetCenterNeighborhoodIndex()) =
						NumericTraits<TPixel>::OneValue();
					m_SumOfWeights += NumericTraits<TPixel>::OneValue();
				}
				else {
					this->operator[](this->GetCenterNeighborhoodIndex()) =
						NumericTraits<TPixel>::ZeroValue();
					m_SumOfWeights -= NumericTraits<TPixel>::OneValue();
				}
			}
		}

		/**
		* Assignment operator
		*/
		Self & operator=(const Self & other) {
			Superclass::operator=(other);
			return *this;
		}

		/**
		* Prints some debugging information
		*/
		virtual void PrintSelf(std::ostream & os, Indent ind) const ITK_OVERRIDE{
			os << ind << "LocalStatisticNeighborhoodBase { this=" << this << "}" << std::endl;
			Superclass::PrintSelf(os, ind.GetNextIndent());
		}

	protected:
		/*
		* Using or not central pixel. False by default
		*/
		bool m_UseCentralPixel;
		/*
		* Sum of neighborhood coefficients
		*/
		TPixel m_SumOfWeights;
	};

	template< typename TPixel, unsigned int VDimension, typename TContainer >
	std::ostream & operator<<(std::ostream & os, const LocalStatisticNeighborhoodBase< TPixel, VDimension, TContainer > & neighborhood)
	{
		os << "Neighborhood: " << neighborhood.GetNameOfClass() << std::endl;
		os << "\tRadius:" << neighborhood.GetRadius() << std::endl;
		os << "\tSize:" << neighborhood.GetSize() << std::endl;
		os << "\tUseCentralPixel:" << neighborhood.GetUseCentralPixel() << std::endl;
		os << "\tSumOfWeights:" << neighborhood.GetSumOfWeights() << std::endl;
		os << "\tDataBuffer:" << neighborhood.GetBufferReference() << std::endl;
		os << "\tData:" << std::endl;
#ifdef DEBUG
		typename LocalStatisticNeighborhoodBase< TPixel, VDimension, TContainer >::SizeType k = neighborhood.GetSize();
		for (SizeValueType i = 0; i < k[0]; ++i)
		{
			os << "\t\t";
			for (SizeValueType j = 0; j < k[1]; ++j)
			{
				os << neighborhood[i*k[1] + j] << "\t";
			}
			os << std::endl;
		}
#endif
		return os;
	}

} // namespace itk

#endif //itkLocalStatisticNeighborhoodBase_hxx
