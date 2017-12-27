/*=========================================================================

  Program:   StokesMatrixImageMultilooking
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) IPMS SB RAS


     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. 

=========================================================================*/

#include <iostream>

#include "otbImage.h"
#include "otbVectorImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbPerBandVectorImageFilter.h"
#include "itkBinShrinkImageFilter.h"
#include "itkTimeProbe.h"


int main(int argc, char* argv[]) 
{
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " InputImage RangeLooks AzimuthLooks" << std::endl
                << "\twhere \'InputImage\' - file name of the input vector image," << std::endl
                << "\t \'RangeLooks\' - number of looks in range directions," << std::endl
                << "\t \'AzimuthLooks\' - number of looks in azimuth directions.";
        return EXIT_SUCCESS;
    }

    const unsigned int Dimension = 2;
    typedef double InputPixelType;
    typedef double OutputPixelType;

    typedef otb::Image<InputPixelType, Dimension> PerBandImageType;
    typedef otb::VectorImage<InputPixelType, Dimension> InputVectorImageType;
    typedef otb::VectorImage<OutputPixelType, Dimension> OutputVectorImageType;
    
    typedef otb::ImageFileReader<InputVectorImageType> ReaderType;
    typedef otb::ImageFileWriter<OutputVectorImageType> WriterType;

    std::string inputFileName = argv[1];
    int rangeLooks = atoi(argv[2]);
    int azimuthLooks = atoi(argv[3]);
    
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(inputFileName);
    
    typedef itk::BinShrinkImageFilter<PerBandImageType, PerBandImageType> BinShrinkFilterType;
    typedef otb::PerBandVectorImageFilter<InputVectorImageType, OutputVectorImageType, 
            BinShrinkFilterType> PerBandShrinkFilterType;
    PerBandShrinkFilterType::Pointer perBandFilter = PerBandShrinkFilterType::New();

    perBandFilter->GetFilter()->SetShrinkFactor(0, rangeLooks);
    perBandFilter->GetFilter()->SetShrinkFactor(1, azimuthLooks);
    perBandFilter->SetInput(reader->GetOutput());
    
    std::string outputFilename = inputFileName.substr(0, inputFileName.find_first_of("_"));
    outputFilename += ("_mls.hdr");

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFilename);
    writer->SetInput(perBandFilter->GetOutput());
    
    itk::TimeProbe clock;
    clock.Start();
    try
    {
        writer->Update();
    }catch(itk::ExceptionObject& e)
    {
        std::cerr << "Error occurred! \n" << e.what();
    }
    clock.Stop();
    std::cout << "Total working time: " << clock.GetTotal() << std::endl;
    
    return EXIT_SUCCESS;
}


