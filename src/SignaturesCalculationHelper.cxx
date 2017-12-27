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

#include <boost/algorithm/string.hpp>
#include "SignaturesCalculationHelper.h"

namespace Helper
{

std::string 
CommandLineHelper::GetOptionFromInput(const std::string option, 
            const std::string synonym, ParserResultType* parserResult)
{
        std::string optString;
        if (parserResult->IsOptionPresent(option))
        {
                optString = parserResult->GetParameterString(option);
        }else
        {
                if (parserResult->IsOptionPresent(synonym))
                {
                        optString = parserResult->GetParameterString(synonym);
                }
        }
        return optString;
}

void CommandLineHelper::InitializeVariables(ParserResultType* parserResult)
{
    m_InputImageFile = this->GetOptionFromInput("--InputFile", "-i", parserResult);
    if( !itksys::SystemTools::FileExists(m_InputImageFile.c_str(), true))
    {
        std::cerr << "File "<< m_InputImageFile <<" does not exist!" << std::flush;
        exit(-1);
    }

    m_InputPOIFile = this->GetOptionFromInput("--InputPOIFile", "-p", parserResult);
    if( !itksys::SystemTools::FileExists(m_InputPOIFile.c_str(), true))
    {
        std::cerr << "File "<< m_InputPOIFile <<" does not exist!" << std::flush;
        exit(-2);
    }

    if ( !this->GetOptionFromInput("--Radius", "-r", parserResult).empty() )
    {
        std::istringstream s(this->GetOptionFromInput("--Radius", "-r", parserResult));
        if (std::isdigit(s.str()[0])){
            s >> m_Radius;
        }else{
            std::cout << "Warning! Radius not represented as digit. ";
            std::cout << "Setting default value: 2.\n";
            m_Radius = 2U;
        }
        
    }

    if ( !this->GetOptionFromInput("--PolMode", "-m", parserResult).empty() )
    {
        std::string polModeString = this->GetOptionFromInput("--PolMode", "-m", parserResult);
        if (polModeString == "crosspolar"){
            m_PolMode = PolarizationMode::CROSS_POLAR;
        }
	   if (polModeString == "any"){
                m_PolMode = PolarizationMode::ANY_POLAR;
		}
    }

    m_OutputDir = this->GetOptionFromInput("--OutputDir", "-o", parserResult);
    if (!m_OutputDir.empty())
    {
        if (!itksys::SystemTools::FileExists(m_OutputDir.c_str(),false))
        {
            itksys::SystemTools::MakeDirectory(m_OutputDir.c_str());
        }
    }

    if ( !this->GetOptionFromInput("--IndexList", "-l", parserResult).empty() )
    {
        m_IndexesString = this->GetOptionFromInput("--IndexList", "-l", parserResult);
    }
    std::transform(m_IndexesString.begin(), m_IndexesString.end(), m_IndexesString.begin(), toupper);
    if (m_IndexesString == "A")
    {   // Calculate all indexes
        m_IndexesString = "FPS";
    }
}

CommandLineHelper::CommandLineHelper(int argc, char* argv[]):
        m_Radius(2),
        m_PolMode(PolarizationMode::CO_POLAR),
        m_IndexesString("FP")
{
    m_IndexNames['F'] = "Fractal";
    m_IndexNames['P'] = "Signature";
    m_IndexNames['S'] = "Sigma";
    
    ParserType::Pointer parser = ParserType::New();
    parser->AddOption("--InputFile", "Input filename.", "-i", 1, true);
    parser->AddOption("--InputPOIFile", "Input POI's filename.", "-p", 1, true);
    parser->AddOption("--Radius", "Neighborhood radius. Default value is 2", "-r", 1, false);
    parser->AddOption("--PolMode", "Polarization mode. Default value is \"copolar\"", "-m", 1, false);
    parser->AddOption("--IndexList", "Index calculation list. Default value is \"FP\"", "-l", 1, false);
    parser->AddOption("--OutputDir", "Output directory name for results. \
        Default value is a directory where input file is located", "-o", 1, false);
    
    ParserResultType::Pointer parserResult = ParserResultType::New();
    try
    {
    	parser->ParseCommandLine(argc, argv, parserResult);
    }catch(CommandLineArgumentParserArgumentErrorException& ex)
    {
    	std::cout << ex.what() << std::endl;
    	exit(-1);
    }
    InitializeVariables(parserResult);
}

size_t CommandLineHelper::ReadPOIs(IndexVectorType& poiVector)
{
    StringListType poiStringList;
    std::ifstream pois;
    pois.open(m_InputPOIFile);
    
    while(!pois.eof())
    {
        std::string poiStr;
        std::getline (pois,poiStr);
        poiStringList.push_back(poiStr);
    }
    pois.close();
    
    for (auto poi : poiStringList) 
    {
        // Extracting POIs
        if ( (poi[0]=='#') || (poi[0]==';') ) continue; //Stripping comments
        if (poi[0] == ' ') boost::trim_left(poi); // trimming leading whitespaces

        std::vector<std::string> poiXY;
        boost::split(poiXY, poi, boost::is_any_of(" ,;\t"),boost::token_compress_on);
        if (poiXY.size() < 2) continue; // Skipping invalid poi
        IndexValueType poiX = atol(poiXY[0].c_str());
        IndexValueType poiY = atol(poiXY[1].c_str());
        poiVector.push_back({{poiX, poiY}});
    }
    return poiVector.size();
}

size_t CommandLineHelper::PrepareOutput(IndexCalculationDescriptorVectorType& IndexCalculationVector)
{
    std::string outputBaseFilename; 
    if (m_OutputDir.empty())
    {
        outputBaseFilename = itksys::SystemTools::GetFilenamePath(m_InputImageFile);
    }else{
        outputBaseFilename = m_OutputDir;
    }
    outputBaseFilename.append("/");
    outputBaseFilename.append(itksys::SystemTools::GetFilenameWithoutLastExtension(m_InputImageFile));

    std::string modeName;
    switch(m_PolMode)
    {
        case CO_POLAR:          modeName = "copolar"; break;
        case CROSS_POLAR:       modeName = "crosspolar"; break;
        case ANY_POLAR:         modeName = "any"; break;
        default:
            break;
    }

    for (auto s : m_IndexesString) {
        std::ostringstream oss;
        oss << outputBaseFilename << "_" << modeName << "_r"
                << m_Radius << "_" << m_IndexNames[s] << ".avrg";

        FileStreamPointerType  newFilePtr = std::make_shared<std::ofstream>(oss.str());
        newFilePtr->width(12);
        newFilePtr->setf(std::ios_base::fixed, std::ios_base::floatfield);

        std::ostringstream ossInfo;
        ossInfo << "#Input image filename: " << m_InputImageFile
                << "\n#Input POI filename: " << m_InputPOIFile
                << "\n#Neighborhood radius is: " << m_Radius
                << "\n#psi\t" << "khi\t" << m_IndexNames[s] << "\t" << "StdDev";
        *newFilePtr << ossInfo.str() << std::endl;

        if (m_IndexNames[s] == "Signature")
        {
            IndexCalculationVector.push_back(std::make_pair(std::move(newFilePtr), 
                static_cast<StatisticImageFunctionPointerType>(MeanFunctionType::New())));
            continue;
        }
        if (m_IndexNames[s] == "Fractal")
        {
            IndexCalculationVector.push_back(std::make_pair(std::move(newFilePtr),
                static_cast<StatisticImageFunctionPointerType>(FractalFunctionType::New())));
            continue;
        }
        if (m_IndexNames[s] == "Sigma")
        {
            IndexCalculationVector.push_back(std::make_pair(std::move(newFilePtr),
                static_cast<StatisticImageFunctionPointerType>(SigmaFunctionType::New())));
            continue;
        }
	}
    std::shared_ptr<RegularNeighborhoodType> neighborhood = 
        std::make_shared<RegularNeighborhoodType>();
    neighborhood->SetRadius(m_Radius);
    for (auto idx : IndexCalculationVector) 
	{
        idx.second->SetNeighborhood(neighborhood);
    }
    return IndexCalculationVector.size();
}

} // end namespace Helper