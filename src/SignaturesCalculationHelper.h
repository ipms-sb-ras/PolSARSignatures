/*=========================================================================

  Program:   SignaturesCalculationHelper
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) IPMS SB RAS

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  

=========================================================================*/

#ifndef SignaturesCalculationHelper_h
#define	SignaturesCalculationHelper_h

#include <vector>
#include <string>
#include <fstream>
#include <map>

#include "otbImage.h"
#include "otbVectorImage.h"
#include "otbWrapperApplicationFactory.h"
#include "otbMultiChannelsPolarimetricSynthesisFilter.h"
#include "otbPolarimetricSynthesisStokesMatrixFilter.h"
#include "otbWrapperApplication.h"
#include "otbCommandLineArgumentParser.h"

// Spatial indexes
#include "itkLocalStatisticImageFunctionBase.hxx"
#include "itkLocalMeanImageFunction.hxx"
#include "itkLocalSigmaImageFunction.hxx"
#include "itkStochasticFractalImageFunction.hxx"
// Neighborhood types
#include "itkRegularNeighborhood.hxx"


namespace Helper
{

#ifdef  USE_DOUBLE_PRECISION 
typedef double  RealPixelType;
#else
typedef float   RealPixelType;
#endif

typedef std::complex<RealPixelType>     ComplexPixelType;

typedef otb::Image<RealPixelType, 2U>           RealImageType;
typedef otb::VectorImage<RealPixelType, 2U>     InputStokesImageType;
typedef otb::VectorImage<ComplexPixelType, 2U>  InputSinclairImageType;

typedef RealImageType::SizeType SizeType;
typedef RealImageType::PixelType PixelType;
typedef RealImageType::ValueType ValueType;
typedef RealImageType::IndexType IndexType;
typedef RealImageType::IndexType::IndexValueType IndexValueType;

typedef std::vector<RealPixelType>  RealVectorType;
typedef std::vector<std::string>    StringListType;
typedef std::vector<IndexType>      IndexVectorType;

typedef itk::NumericTraits<RealImageType>::AccumulateType AccumulateType;
typedef itk::NumericTraits<PixelType>::RealType RealType;

typedef itk::LocalStatisticImageFunctionBase<RealImageType, RealType> StatisticImageFunctionType;
typedef StatisticImageFunctionType::Pointer StatisticImageFunctionPointerType;
typedef itk::LocalMeanImageFunction < RealImageType, RealType > MeanFunctionType;
typedef itk::LocalSigmaImageFunction < RealImageType, RealType > SigmaFunctionType;
typedef itk::StochasticFractalImageFunction< RealImageType, RealType > FractalFunctionType;

typedef itk::RegularNeighborhood<PixelType, 2> RegularNeighborhoodType;

typedef otb::MultiChannelsPolarimetricSynthesisFilter<InputSinclairImageType,
    RealImageType> SinclairFilterType;

typedef otb::PolarimetricSynthesisStokesMatrixFilter<InputStokesImageType,
    RealImageType> StokesFilterType;

typedef std::shared_ptr<std::ofstream> FileStreamPointerType;
typedef std::pair<FileStreamPointerType, StatisticImageFunctionPointerType> IndexCalculationDescriptorType;
typedef std::vector<IndexCalculationDescriptorType> IndexCalculationDescriptorVectorType;

typedef std::map<char, std::string> IndexNameType;

typedef enum {
    ANY_POLAR,
    CO_POLAR,
    CROSS_POLAR,
} PolarizationMode;


class CommandLineHelper {
    typedef otb::CommandLineArgumentParser ParserType;
    typedef otb::CommandLineArgumentParseResult ParserResultType;
    
public:
    CommandLineHelper(int argc, char* argv[]);
    
    virtual ~CommandLineHelper(){};
    
    size_t GetRadius() const  { return m_Radius;};
    PolarizationMode GetPolarizationMode() const { return m_PolMode;};
    
    std::string GetInputImageFilename() const
    {
        return m_InputImageFile;
    }
    
    size_t ReadPOIs(IndexVectorType& poiVector);
    size_t PrepareOutput(IndexCalculationDescriptorVectorType& IndexCalculationVector);

private:
    CommandLineHelper(const CommandLineHelper& orig){};
        
    std::string GetOptionFromInput(const std::string option, const std::string synonym,
        ParserResultType* parserResult);
    void InitializeVariables(ParserResultType* parserResult);

    std::string m_InputImageFile;
    std::string m_InputPOIFile;
    std::string m_OutputDir;
    std::string m_IndexesString;
    size_t m_Radius;
    PolarizationMode m_PolMode;
    IndexNameType m_IndexNames;
};

} // end namespace Helper
 
#endif	/* SignaturesCalculationHelper_h */

