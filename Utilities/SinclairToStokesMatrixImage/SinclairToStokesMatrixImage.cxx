/*=========================================================================

  Program:   
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) IPMS SD RAS


     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <map>

#include "itksys/SystemTools.hxx"
#include "itkDirectory.h"
#include "itkMacro.h"

#include <otbImage.h>
#include <otbVectorImage.h>
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbSinclairToStokesMatrixFunctor.h"
#include "otbSinclairImageFilter.h"

typedef std::vector<std::string> InputImagesList;

typedef std::map < std::string, std::string > FilePolarizationMap;

size_t ScanInputDirectoryForImages(const std::string& inputDir, FilePolarizationMap *imgMap) {
    static const std::string fileExt(".hdr");

    if (!itksys::Directory::GetNumberOfFilesInDirectory(inputDir))
        return 0; // Empty directory

    itksys::Directory dir;
    if (!dir.Load(inputDir.c_str())) {
        itkGenericExceptionMacro( << "Cannot open directory " << inputDir << std::endl);
    }
    std::string dirElement;
    for (unsigned long i = 2; i < dir.GetNumberOfFiles(); i++) {
        dirElement = inputDir + std::string("/") + std::string(dir.GetFile(i));
        dirElement = itksys::SystemTools::ConvertToOutputPath(dirElement);
        if (!itksys::SystemTools::FileIsDirectory(dirElement)) {
            if (itksys::SystemTools::GetFilenameLastExtension(dirElement) == fileExt) {
                std::string filename = itksys::SystemTools::GetFilenameName(dirElement);
                if (filename.find("HH") != std::string::npos) {
                    (*imgMap)["HH"] = dirElement;
                }
                if (filename.find("HV") != std::string::npos) {
                    (*imgMap)["HV"] = dirElement;
                }
                if (filename.find("VH") != std::string::npos) {
                    (*imgMap)["VH"] = dirElement;
                }
                if (filename.find("VV") != std::string::npos) {
                    (*imgMap)["VV"] = dirElement;
                }
            }
        }
    }
    return imgMap->size();
}

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " inputDir" << std::endl
                << "where \'inputDir\' is a directory were images are located";
        return EXIT_SUCCESS;
    }

    const unsigned int Dimension = 2;
    typedef std::complex <double> InputPixelType;
    typedef double OutputRealPixelType;

    typedef otb::Image<InputPixelType, Dimension> InputImageType;
    typedef otb::VectorImage<OutputRealPixelType, Dimension> OutputRealImageType;
    
    typedef otb::ImageFileReader<InputImageType> ReaderType;
    typedef otb::ImageFileWriter<OutputRealImageType> WriterType;
    
    typedef otb::Functor::SinclairToStokesMatrixFunctor<InputImageType::PixelType,
                InputImageType::PixelType,
                InputImageType::PixelType,
                InputImageType::PixelType,
                OutputRealImageType::PixelType> FunctorType;
    
    typedef otb::SinclairImageFilter<InputImageType, InputImageType, 
            InputImageType, InputImageType, 
            OutputRealImageType, FunctorType> FilterType;

    FilterType::Pointer filter = FilterType::New();
    ReaderType::Pointer readerHH = ReaderType::New();
    ReaderType::Pointer readerHV = ReaderType::New();
    ReaderType::Pointer readerVH = ReaderType::New();
    ReaderType::Pointer readerVV = ReaderType::New();

    const std::string inputDir = argv[1];
    FilePolarizationMap imageFilenames;
    
    size_t numFiles = ScanInputDirectoryForImages(inputDir, &imageFilenames);
    if (numFiles != 4){
        std::cerr << "Not enough images for calculation!" << std::endl;
        return EXIT_FAILURE;
    }
    
    readerHH->SetFileName(imageFilenames["HH"]);
    readerHV->SetFileName(imageFilenames["HV"]);
    readerVH->SetFileName(imageFilenames["VH"]);
    readerVV->SetFileName(imageFilenames["VV"]);
    filter->SetInputHH(readerHH->GetOutput());
    filter->SetInputHV(readerHV->GetOutput());
    filter->SetInputVH(readerVH->GetOutput());
    filter->SetInputVV(readerVV->GetOutput());

    filter->UpdateOutputInformation();

    typename WriterType::Pointer writer = WriterType::New();
    std::string outputFilename = itksys::SystemTools::GetFilenameName(imageFilenames["HH"]);
    outputFilename = inputDir + "/" + outputFilename.substr(0, outputFilename.find_first_of("-"));
    outputFilename += ("-QP_slc.hdr");
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
