/*=========================================================================

  Program:   Fractal Polarization Signature
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) IPMS SB RAS

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  

=========================================================================*/


#include "otbImageFileReader.h"

#include "SignaturesCalculationHelper.h"
#include "SignaturesGenerator.h"

std::string GetInputImagePixelType(const std::string& filename);

int main(int argc, char* argv[]) {
    typedef otb::ImageFileReader<Helper::InputSinclairImageType> SinclairImageReaderType;
    typedef otb::ImageFileReader<Helper::InputStokesImageType> StokesImageReaderType;

    const double step = 3.0; // calculation step

    Helper::CommandLineHelper clParser(argc, argv);

    std::string imageFilename = clParser.GetInputImageFilename();
    Helper::IndexVectorType poisVector;
    clParser.ReadPOIs(poisVector);

    try {
        std::string imagePixelTypeString = GetInputImagePixelType(imageFilename);

        Helper::IndexCalculationDescriptorVectorType indexCalculationVector;
        clParser.PrepareOutput(indexCalculationVector);
        // Case of Sinclair 4 channels (2x2 matrix) complex image
        if (imagePixelTypeString == "complex") {
            SinclairImageReaderType::Pointer reader =
                    SinclairImageReaderType::New();
            reader->SetFileName(imageFilename);
            reader->Update();
            Helper::InputSinclairImageType::Pointer image = reader->GetOutput();
            if (image->GetNumberOfComponentsPerPixel() != 4){
                itkGenericExceptionMacro( << "Number of pixel components must be equal 4!");
            } 
            Helper::SignaturesGenerator<Helper::InputSinclairImageType,
                    Helper::SinclairFilterType>(image, poisVector, indexCalculationVector,
                    clParser.GetPolarizationMode(), step);
            return EXIT_SUCCESS;
        }
        // Case of Stokes 16 channels (4x4 matrix) real image
        if (imagePixelTypeString == "vector") {
            StokesImageReaderType::Pointer reader =
                    StokesImageReaderType::New();
            reader->SetFileName(imageFilename);
            reader->Update();
            Helper::InputStokesImageType::Pointer image = reader->GetOutput();
            if (image->GetNumberOfComponentsPerPixel() != 16){
                itkGenericExceptionMacro( << "Number of pixel components must be equal 16!");
            } 
            Helper::SignaturesGenerator<Helper::InputStokesImageType,
                    Helper::StokesFilterType>(image, poisVector, indexCalculationVector,
                    clParser.GetPolarizationMode(), step);
            return EXIT_SUCCESS;
        } else {
            itkGenericExceptionMacro( << "Unsupported type of the image:\t" << imageFilename);
        }
    } catch (itk::ExceptionObject& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

std::string GetInputImagePixelType(const std::string& filename) {
    
    std::string fname(filename);
    if (itksys::SystemTools::GetFilenameLastExtension(filename) == ".hdr"){ // ENVI file. 
        fname = itksys::SystemTools::GetFilenamePath(filename) + '/'
               +itksys::SystemTools::GetFilenameWithoutLastExtension(filename);
    }
    otb::ImageIOBase::Pointer imageIO =
            otb::ImageIOFactory::CreateImageIO(fname.c_str(),
            otb::ImageIOFactory::ReadMode);

    if (!imageIO) {
        itkGenericExceptionMacro( << "Unknown format of the input image:\t" << fname);
    }
    imageIO->SetFileName(fname);
    imageIO->ReadImageInformation();
    return imageIO->GetPixelTypeAsString(imageIO->GetPixelType());
}

