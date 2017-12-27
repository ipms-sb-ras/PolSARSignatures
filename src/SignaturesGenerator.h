/*=========================================================================

  Program:   SignaturesGenerator
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) IPMS SB RAS

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  

=========================================================================*/

#ifndef SignaturesGenerator_h
#define	SignaturesGenerator_h

#include "ezETAProgressBar.hpp"
#include "SignaturesCalculationHelper.h"

namespace Helper{

template <class TInpuImage, class TFilter>
void SignaturesGenerator(const typename TInpuImage::Pointer inputImage, 
        const IndexVectorType& pois,
        const IndexCalculationDescriptorVectorType& indexCalculationVector, 
        PolarizationMode polMode, double step = 3.)
{
    size_t maxXSize = inputImage->GetLargestPossibleRegion().GetSize(0);
    size_t maxYSize = inputImage->GetLargestPossibleRegion().GetSize(1);

    size_t radius = indexCalculationVector[0].second->GetNeighborhoodPointer()->GetRadius()[0];
    IndexType constrains = {{static_cast<IndexValueType>(maxXSize-radius), 
            static_cast<IndexValueType>(maxYSize-radius)}};

    RealType psiStart, psiStop, psiStep;
    RealType khiStart, khiStop, khiStep;
    psiStart = 0.;
    psiStop = 180.;
	khiStart = -45.;
	khiStop = 45.;
    psiStep = khiStep = step;

    std::cout.precision(2);
    
    ez::ezETAProgressBar pg(static_cast<size_t>((psiStop-psiStart)/psiStep) + 1);
    pg.start();
    
    typedef typename TFilter::Pointer SynthesisFilter;
    for (RealType psi = psiStart; psi <= psiStop; psi+=psiStep) 
    {
        for (RealType khi = khiStart; khi <= khiStop; khi+=khiStep) 
        {
            SynthesisFilter polarimetricSynthesisFilter = TFilter::New();
            polarimetricSynthesisFilter->SetInput(inputImage);
            polarimetricSynthesisFilter->SetPsiI(psi);
            polarimetricSynthesisFilter->SetKhiI(khi);
            switch (polMode) 
            {
                case CO_POLAR:
                    polarimetricSynthesisFilter->ForceCoPolar();
                    break;
                case CROSS_POLAR:
                    polarimetricSynthesisFilter->ForceCrossPolar();
                    break;
                default: // ANY_POLAR
                    polarimetricSynthesisFilter->SetPsiI(0.);
                    polarimetricSynthesisFilter->SetKhiI(0.);
                    polarimetricSynthesisFilter->SetPsiR(psi);
                    polarimetricSynthesisFilter->SetKhiR(khi);
                    polarimetricSynthesisFilter->SetMode(ANY_POLAR);
            }
            
            polarimetricSynthesisFilter->Update();
            RealImageType::Pointer synthesizedImage =
                    polarimetricSynthesisFilter->GetOutput();

            for (auto idxCalc : indexCalculationVector) {
                idxCalc.second->SetInputImage(synthesizedImage);
                double valSum = 0.;
                double valSquaredSum = 0.;
                size_t counter = 0;
                double val;
                for (auto poi: pois)  {
                    if (    poi[0] > constrains[0] || poi[0] < 0 
                        ||  poi[1] > constrains[1] || poi[1] < 0) 
                    {
                        std::cerr << "Neighborhood of POI [" << poi[0] << ", "
                                << poi[0] << "] exceeds the image  boundaries! This POI will be skipped.\n";
                        continue;
                    }

                    val = static_cast<double>(idxCalc.second->EvaluateAtIndex(poi));
                    if ( !std::isnan(val) ){
                        valSum += val;
                        valSquaredSum += (val*val);
                        ++counter;
                    }
                }
                double n = static_cast<double>(counter);
                double mean = valSum / n;
                double stdDev = std::sqrt((valSquaredSum - n * mean * mean) / (n - 1.));
                
                *idxCalc.first << psi << "\t" << khi << "\t" << mean << "\t" << stdDev << std::endl;
            }
        }
        ++pg;
    }
}

} // end namespace Helper

#endif // SignaturesGenerator_h