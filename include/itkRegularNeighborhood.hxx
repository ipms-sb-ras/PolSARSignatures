/*=========================================================================

Program:
Language:  C++
Date:      $Date: $
Version:   $Revision: $

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

@file    itkRegularNeighborhood.hxx

=========================================================================*/
#ifndef itkRegularNeighborhood_hxx
#define itkRegularNeighborhood_hxx

#include <numeric>
#include "itkLocalStatisticNeighborhoodBase.hxx"
#include "itkNeighborhoodAllocator.h"

namespace itk {

	/**
	* \class RegularNeighborhood
	*
	* \brief A Neighborhood for using in spatial statistics calculation.
	*
	*/
	template< typename TPixel, unsigned int VDimension = 2,
		typename TAllocator = NeighborhoodAllocator< TPixel > >
	class RegularNeighborhood :
		public itk::LocalStatisticNeighborhoodBase < TPixel, VDimension, TAllocator > {
	public:
		/** Standard typedefs */
		typedef RegularNeighborhood Self;
		typedef LocalStatisticNeighborhoodBase< TPixel, VDimension, TAllocator > Baseclass;

		itkTypeMacro(RegularNeighborhood, Baseclass);

		void SetRadius(const SizeValueType s)
		{
			Baseclass::SetRadius(s);
		}

		void SetRadius(const typename Baseclass::SizeType & r) {
			Baseclass::SetRadius(r);
			// Filling all with One
			std::fill(this->Begin(), this->End(), NumericTraits<TPixel>::OneValue());

			// Filling Center of Neighborhood with corresponding value
			if (!this->GetUseCentralPixel()) {
				this->operator[](this->GetCenterNeighborhoodIndex()) =
					NumericTraits<TPixel>::ZeroValue();
			}
			// Calculating sum of all neighborhood's values
			this->m_SumOfWeights = std::accumulate(this->Begin(), this->End(),
				NumericTraits<TPixel>::ZeroValue());
		}

		/**
		* Prints some debugging information
		*/
		virtual void PrintSelf(std::ostream & os, Indent i) const ITK_OVERRIDE{
			os << i << "RegularNeighborhood { this=" << this << "}" << std::endl;
			Baseclass::PrintSelf(os, i.GetNextIndent());
		}

	private:
	};
} // namespace itk

#endif //itkRegularNeighborhood_hxx
