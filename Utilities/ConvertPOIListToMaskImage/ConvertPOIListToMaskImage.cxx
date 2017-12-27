#include <boost/algorithm/string.hpp>
#include "otbImage.h"
#include "otbImageFileWriter.h"


typedef otb::Image<unsigned char, 2U> MaskImageType;
typedef itk::Index<2> IndexType;
typedef IndexType::IndexValueType IndexValueType;
typedef std::vector<std::string> StringListType;
typedef std::vector<IndexType> IndexVectorType;

size_t ReadPOIs(const std::string& inputPOIFile, IndexVectorType& poiVector)
{
    std::ifstream pois;
    pois.open(inputPOIFile);
	while(!pois.eof())
    {
        std::string poiStr;
        std::getline (pois,poiStr);
		// Extracting POIs
		if ((poiStr[0] == '#') || (poiStr[0] == ';')) continue; //Stripping comments
		if (poiStr[0] == ' ') boost::trim_left(poiStr); // trimming leading whitespaces

		StringListType poiXY;
		boost::split(poiXY, poiStr, boost::is_any_of(" ,;\t"), boost::token_compress_on);
		if (poiXY.size() < 2) continue; 
		IndexValueType poiX = atol(poiXY[0].c_str());
		IndexValueType poiY = atol(poiXY[1].c_str());
		poiVector.push_back({ { poiX, poiY } });
    }
    pois.close();
    return poiVector.size();
}

IndexType GetImageSize(const std::string& inputPOIFile)
{
	std::ifstream pois;
	pois.open(inputPOIFile);

	std::string poiStr;
	while (!pois.eof())
	{
		std::getline(pois, poiStr);
		std::size_t found = poiStr.find("File Dimension");
		if (found != std::string::npos) {
			std::cout << poiStr << std::endl;
			break;
		}
	}
	pois.close();
	StringListType stringToParse;
	boost::split(stringToParse, poiStr, boost::is_any_of(":"), boost::token_compress_on);
	StringListType sizes;
	boost::split(sizes, stringToParse[1], boost::is_any_of("x"), boost::token_compress_on);
	boost::trim_left(sizes[0]);
	boost::trim_left(sizes[1]);
	return IndexType{ atol(sizes[0].c_str()), atol(sizes[1].c_str()) };
}

int main(int argc, char* argv[])
{
	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " InputPOIFile\n";
		return EXIT_FAILURE;
	}
	IndexType imageSize = GetImageSize(argv[1]);
	IndexVectorType poiList;
	if (!ReadPOIs(argv[1], poiList))
	{
		std::cerr << "No POIs?\n";
		return EXIT_FAILURE;
	}
		
	MaskImageType::IndexType start;
	start.Fill(0);

	MaskImageType::SizeType size;
	size[0] = imageSize[0];
	size[1] = imageSize[1];

	MaskImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);

	MaskImageType::Pointer image = MaskImageType::New();
	image->SetRegions(region);
	image->Allocate();
	image->FillBuffer(0);

	for (auto index : poiList) {
		image->SetPixel(index, 1);
	}

	std::string outputFilename(argv[1]);
	outputFilename = outputFilename.substr(0, outputFilename.rfind('.'));
	outputFilename.append(".hdr");
	std::cout << outputFilename << std::endl;
	otb::ImageFileWriter<MaskImageType>::Pointer writer = otb::ImageFileWriter<MaskImageType>::New();
	writer->SetFileName(outputFilename);
	writer->SetInput(image);
	writer->Update();

	return EXIT_SUCCESS;
}
