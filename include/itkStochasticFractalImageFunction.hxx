/*=========================================================================
 
  Program:   
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.*
 
    Code for fractal dimension calculation is taken from file 
    itkStochasticFractalDimensionImageFilter.hxx
    See http://www.insight-journal.org/browse/publication/318
 * 
 *=========================================================================*/
#ifndef itkStochasticFractalImageFunction_hxx
#define itkStochasticFractalImageFunction_hxx

#include "itkLocalStatisticImageFunctionBase.hxx"
#include "itkConstNeighborhoodIterator.h"


namespace itk
{
	/**
	 * \class StochasticFractalImageFunction
	 * \brief Calculate the fractal dimension of the given neighborhood
	 *
	 * This class is templated over the input image type and the
	 * coordinate representation type (e.g. float or double).
	 *
	 * \ingroup ImageFunctions
	 * \ingroup ITKImageFunction
	 */
	template< typename TInputImage, typename TCoordRep = float >
	class StochasticFractalImageFunction :
		public LocalStatisticImageFunctionBase < TInputImage, TCoordRep >
	{
	public:
		/** Standard class typedefs. */
		typedef StochasticFractalImageFunction Self;
		typedef LocalStatisticImageFunctionBase< TInputImage, TCoordRep > Baseclass;

		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		typedef typename Baseclass::InputImageType InputImageType;
                typedef typename Baseclass::RealType RealType;
                typedef typename Baseclass::IndexType IndexType;

		/** Run-time type information (and related methods). */
		itkTypeMacro(StochasticFractalImageFunction, LocalStatisticImageFunctionBase);

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** ImageDimension constants */
		itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

		/** Evalulate the function at specified index */
		virtual RealType EvaluateAtIndex(const IndexType & index) const ITK_OVERRIDE
		{
			typedef typename InputImageType::PixelType  InputPixelType;
			typedef typename InputImageType::PointType  PointType;

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

			const InputImageType *inputImage = this->GetInputImage();

			/** Turning on central pixel */
			this->m_NeighborhoodPointer->UseCentralPixelOn();

			typename InputImageType::SpacingType spacing = inputImage->GetSpacing();
			RealType minSpacing = spacing[0];
			for (unsigned int d = 0; d < ImageDimension; d++)
			{
                            if (spacing[d] < minSpacing)
                            {
                                minSpacing = spacing[d];
                            }
			}

			std::vector< RealType > distances;
			std::vector< RealType > distancesFrequency;
			std::vector< RealType > averageAbsoluteIntensityDifference;

			ConstNeighborhoodIterator< InputImageType >
				It(this->m_NeighborhoodPointer->GetRadius(), inputImage, inputImage->GetBufferedRegion());
			// Set the iterator at the desired location
			It.SetLocation(index);

			for (unsigned int i = 0; i < It.GetNeighborhood().Size(); i++)
			{
				bool IsInBounds1;
				InputPixelType pixel1 = It.GetPixel(i, IsInBounds1);

				if (!IsInBounds1)
				{
					continue;
				}
				PointType point1;
				inputImage->TransformIndexToPhysicalPoint(It.GetIndex(i), point1);

				for (unsigned int j = 0; j < It.GetNeighborhood().Size(); j++)
				{
					if (i == j)
					{
						continue;
					}

					bool           IsInBounds2;
					InputPixelType pixel2 = It.GetPixel(j, IsInBounds2);

					if (!IsInBounds2)
					{
						continue;
					}
					PointType point2;
					inputImage->TransformIndexToPhysicalPoint(It.GetIndex(j), point2);

					const RealType distance = point1.SquaredEuclideanDistanceTo(point2);

					bool distanceFound = false;
					for (unsigned int k = 0; k < distances.size(); k++)
					{
						if (vnl_math_abs(distances[k] - distance) < 0.5 * minSpacing)
						{
							distancesFrequency[k]++;
							averageAbsoluteIntensityDifference[k] += vnl_math_abs(pixel1 - pixel2);
							distanceFound = true;
							break;
						}
					}

					if (!distanceFound)
					{
						distances.push_back(distance);
						distancesFrequency.push_back(1);
						averageAbsoluteIntensityDifference.push_back(vnl_math_abs(pixel1 - pixel2));
					}
				}
			}

			RealType sumY = 0.0;
			RealType sumX = 0.0;
			RealType sumXY = 0.0;
			RealType sumXX = 0.0;

			for (unsigned int k = 0; k < distances.size(); k++)
			{
				if (distancesFrequency[k] == 0)
				{
					continue;
				}

				averageAbsoluteIntensityDifference[k] /= static_cast<RealType>(distancesFrequency[k]);
				averageAbsoluteIntensityDifference[k] = std::log(averageAbsoluteIntensityDifference[k]);

				const RealType distance = std::log(std::sqrt(distances[k]));

				sumY += averageAbsoluteIntensityDifference[k];
				sumX += distance;
				sumXX += (distance * distance);
				sumXY += (averageAbsoluteIntensityDifference[k] * distance);
			}

			const RealType N = static_cast<RealType>(distances.size());
			const RealType slope = (N * sumXY - sumX * sumY) / (N * sumXX - sumX * sumX);

			RealType fd = static_cast<RealType>(3.0 - slope);
			return fd;
		}

		void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE
		{
			Baseclass::PrintSelf(os, indent);
		}

	protected:
		StochasticFractalImageFunction() : Baseclass(){};
		~StochasticFractalImageFunction(){}

	private:
		StochasticFractalImageFunction(const Self &); //purposely not implemented
		void operator=(const Self &);    //purposely not implemented
	};
} // end namespace itk


#endif // itkStochasticFractalImageFunction_hxx
