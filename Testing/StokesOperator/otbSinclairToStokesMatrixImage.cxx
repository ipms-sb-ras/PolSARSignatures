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


#include <otbImage.h>
#include <otbVectorImage.h>
#include "otbImageFileWriter.h"
#include "otbSinclairToStokesMatrixFunctor.h"
#include "otbSinclairImageFilter.h"

typedef std::complex <double>   InputPixelType;
typedef double                  OutputRealPixelType;

typedef otb::Image<InputPixelType>              InputImageType;
typedef otb::VectorImage<OutputRealPixelType>   OutputRealImageType;


void CreateImage(const InputPixelType fillVal, InputImageType::Pointer image);

int main(int argc, char* argv[]) 
{
    typedef otb::Functor::SinclairToStokesMatrixFunctor<InputImageType::PixelType,
                InputImageType::PixelType,
                InputImageType::PixelType,
                InputImageType::PixelType,
                OutputRealImageType::PixelType> FunctorType;
    
    typedef otb::SinclairImageFilter<InputImageType, InputImageType, 
            InputImageType, InputImageType, 
            OutputRealImageType, FunctorType> FilterType;

    FilterType::Pointer filter = FilterType::New();
    InputImageType::Pointer imageHH = InputImageType::New();
    InputImageType::Pointer imageHV = InputImageType::New();
    InputImageType::Pointer imageVV = InputImageType::New();
    
    CreateImage(InputPixelType(1.0,0.0), imageHH);
    CreateImage(InputPixelType(0.0,0.0), imageHV);
    CreateImage(InputPixelType(1.0,0.0), imageVV);

    filter->SetInputHH(imageHH);
    filter->SetInputHV(imageHV);
    filter->SetInputVH(imageHV);
    filter->SetInputVV(imageVV);

    filter->UpdateOutputInformation();

    typedef otb::ImageFileWriter<OutputRealImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    std::string outputFilename = "~/stokes_trihedral.hdr";
    writer->SetFileName(outputFilename);
    writer->SetInput(filter->GetOutput());
    try{
        writer->Update();
    }catch(itk::ExceptionObject &e){
        std::cout << "Exception caught!" << std::endl
                << e.what() << std::endl;
    }
    return EXIT_SUCCESS;
}

void CreateImage(const InputPixelType fillVal, InputImageType::Pointer image)
{
    InputImageType::IndexType start;
    start.Fill(0);

    InputImageType::SizeType size;
    size.Fill(200);

    InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(fillVal);
}
