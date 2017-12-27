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

#include <fstream>

#include "otbPolarimetricSynthesisStokesMatrixFilter.h"
#include <otbImage.h>
#include <otbVectorImage.h>
#include "otbImageFileReader.h"
#include "itkTimeProbe.h"

typedef double   RealType;

typedef otb::Image<RealType>         OutputRealImageType;
typedef otb::VectorImage<RealType>   InputVectorImageType;

void CalculateSignature(const InputVectorImageType::Pointer image, int mode, std::ofstream& out);

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " ImageFilename" << std::endl
                << "\twhere \'ImageFilename\' name of the input vector image file \n"
                <<"\twith pixel representing Stokes scattering operator (matrix).";
        return EXIT_SUCCESS;
    }

    typedef otb::ImageFileReader<InputVectorImageType> ImageReaderType;
    
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(argv[1]);
    try{
        reader->Update();
    }catch(itk::ExceptionObject& e){
        std::cerr << "Error occurred! \n" << e.what();
        return EXIT_FAILURE;
    }
    
    const InputVectorImageType::Pointer inputImage = reader->GetOutput();

    itk::TimeProbe clock;
    clock.Start();
    
    std::string results_file_1(argv[1]);
    results_file_1 += "_copolar";
    std::ofstream outCo1(results_file_1.c_str());
    CalculateSignature(inputImage, 1 /*co_polar*/, outCo1);
    clock.Stop();
    std::cout << "Co-polar calculation time:\t" << clock.GetMean() << std::endl;
    clock.Start();
    
    std::string results_file_2(argv[1]);
    results_file_2 += "_crosspolar";
    std::ofstream outCo2(results_file_2.c_str());
    CalculateSignature(inputImage, 2 /*cross_polar*/, outCo2);
    clock.Stop();
    std::cout << "Cross-polar calculation time:\t" << clock.GetMean() << std::endl;
    std::cout << "Total calculation time:\t" << clock.GetTotal() << std::endl;
    
    return EXIT_SUCCESS;
}

void CalculateSignature(const InputVectorImageType::Pointer image, int mode, std::ofstream& out)
{
    typedef otb::PolarimetricSynthesisStokesMatrixFilter<InputVectorImageType, 
            OutputRealImageType> FilterType;

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(image);

    size_t maxXSize = image->GetLargestPossibleRegion().GetSize(0);
    size_t maxYSize = image->GetLargestPossibleRegion().GetSize(1);
    itk::Index<2> pixelIndex = {{maxXSize/2,maxYSize/2}};

    out << "# psi \t khi \t ps \n ";
    for (float psi = 0; psi <= 180.; psi+=1.) {
        filter->SetPsiI( psi );
        for (float khi = -45.; khi <= 45.; khi+=1.) {
            filter->SetKhiI( khi);
            (mode == 1) ? filter->ForceCoPolar() : filter->ForceCrossPolar();
            filter->Update();
            OutputRealImageType::Pointer m_SynthesizedImage = filter->GetOutput();
            RealType pixelVal = m_SynthesizedImage->GetPixel(pixelIndex);
            out << psi << "\t" << khi << "\t" << pixelVal << "\n";
        }
    }
}

