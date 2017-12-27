#include <iostream>
#include <memory>

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"

#include "itkRegularNeighborhood.hxx"
#include "itkQueenNeighborhood.hxx"
#include "itkRookNeighborhood.hxx"
#include "itkBishopNeighborhood.hxx"

#include "itkLacunarityDongImageFunction.hxx"

typedef float PixelType;
typedef otb::Image<PixelType, 2> ImageType;
typedef otb::Image<PixelType, 2> OutputImageType;

typedef otb::ImageFileReader<ImageType> ReaderType;
typedef otb::ImageFileWriter<OutputImageType> WriterType;

typedef ImageType::RegionType RegionType;
typedef itk::ImageRegionConstIterator<ImageType> ImageRegionConstIteratorType;
typedef RegionType::SizeType SizeType;
typedef RegionType::IndexType IndexType;
typedef SizeType::SizeValueType SizeValueType;

typedef itk::LocalStatisticNeighborhoodBase<PixelType, 2>	NeighborhoodType;
typedef itk::RegularNeighborhood<ImageType::PixelType, 2>	RegularNeighborhoodType;
typedef itk::QueenNeighborhood<ImageType::PixelType, 2>		QueenNeighborhoodType;
typedef itk::RookNeighborhood<ImageType::PixelType, 2>		RookNeighborhoodType;
typedef itk::BishopNeighborhood<ImageType::PixelType, 2>	BishopNeighborhoodType;

typedef itk::LacunarityDongImageFunction<ImageType>			LacunarityDongImageFunctionType;

typedef std::shared_ptr<NeighborhoodType> NeighborhoodPointerType;

void CreateTestImage(ImageType::Pointer image, const size_t size_);

void PrintRegion(ImageType::Pointer image, const RegionType& region);
void PrintImagePixelValues(ImageType::Pointer image);

void TestNeighborhoodSelection(SizeValueType windowRadius, SizeValueType glBoxRadius, const IndexType& index);


int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " InputImage " << std::endl
			<< "\twhere \'InputImage\' - file name of the input image." << std::endl;
		return EXIT_FAILURE;
	}

	const unsigned int radius[] = { 1, 2, 3, 4, 5 };
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);
	reader->Update();
	ImageType* image = reader->GetOutput();

	typedef LacunarityDongImageFunctionType ImageFunctionType;
	std::shared_ptr<NeighborhoodType> neighborhood = std::make_shared<BishopNeighborhoodType>();

	ImageFunctionType::Pointer function = ImageFunctionType::New();
	function->SetInputImage(image);
	function->SetNeighborhood(neighborhood);
	function->SetWindowRadius(20);
	std::cout << "\n\nNeighborhood name: " << function->GetNeighborhoodPointer()->GetNameOfClass() << std::endl;
	ImageType::RegionType::IndexType index = { { 110, 110 } };
	for (auto i : radius)
	{
		function->SetGlidingBoxRadius(i);
		std::cout << i << ":\t" << function->EvaluateAtIndex(index) << "\n";
		// system("pause");
	}
	//IndexType index = { {5, 5} };
	//TestNeighborhoodSelection(5, 5, index);
	return EXIT_SUCCESS;
}


void CreateTestImage(ImageType::Pointer image, const size_t size_)
{
	ImageType::IndexType start;
	start.Fill(0);

	ImageType::SizeType size;
	size.Fill(size_);

	ImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);

	image->SetRegions(region);
	image->Allocate();

	itk::ImageRegionIterator<ImageType> imageIterator(image, region);

	size_t counter = 1;
	while (!imageIterator.IsAtEnd())
	{
		imageIterator.Set(counter);
		++imageIterator;
		++counter;
	}
	std::cout << "Pixels count: " << counter << std::endl;
}

void CreateEmptyImage(ImageType::Pointer image, const size_t size_)
{
	ImageType::IndexType start;
	start.Fill(0);

	ImageType::SizeType size;
	size.Fill(size_);

	ImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);

	image->SetRegions(region);
	image->Allocate();
}


void PrintRegion(ImageType::Pointer image, const RegionType& region)
{
	ImageRegionConstIteratorType regionIterator(image, region);
	size_t counter = 0;
	ImageType::SizeValueType regionSize = region.GetSize(0);

	while (!regionIterator.IsAtEnd())
	{
		std::cout << regionIterator.Get() << "  ";
		++regionIterator;
		++counter;
		if (counter / regionSize)
		{
			counter = 0;
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

void PrintImagePixelValues(ImageType::Pointer image)
{
	ImageRegionConstIteratorType inputIt(image, image->GetLargestPossibleRegion());
	size_t counter = 0;
	ImageType::SizeValueType regionSize = image->GetLargestPossibleRegion().GetSize()[0];
	inputIt.GoToBegin();
	while (!inputIt.IsAtEnd())
	{
		std::cout << inputIt.Get() << "  ";
		++inputIt;
		++counter;
		if (counter / regionSize)
		{
			counter = 0;
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

void TestNeighborhoodSelection(SizeValueType windowRadius, SizeValueType glBoxRadius, const IndexType& index)
{
	const PixelType SIZE = 15;
	ImageType::Pointer image = ImageType::New();

	CreateTestImage(image, SIZE);
	PrintImagePixelValues(image);
	std::cout << "Pixel value at index " << index << " is " << image->GetPixel(index) << std::endl;

	//Top left corner of the calculation region (window)
	IndexType regionIndex =
	{ { index[0] - static_cast<itk::IndexValueType>(windowRadius),
		index[1] - static_cast<itk::IndexValueType>(windowRadius) } };

	SizeType regionSize;
	regionSize.Fill(2*windowRadius+1);

	//Set a region for gliding box movements which centered at the index
	RegionType region(regionIndex, regionSize);
	std::cout << "region for gliding box movements:\n ";
	region.Print(std::cout);

	// Check if it lies within buffered region
	if (region.IsInside(image->GetLargestPossibleRegion()))
	{	// Normally it should not occurs
		std::cerr << "Window region is outside of buffered region!\n"
			<< "Window center location:\t" << index
			<< "\nWindow size:\t\t\t" << region.GetSize()
			<< "\nWindow corner:\t\t" << region.GetIndex() << "\n";
	}

	SizeType glidingBoxRadius;
	glidingBoxRadius.Fill(glBoxRadius);

	SizeValueType glBoxMovementWindowRadius = windowRadius - glBoxRadius;
	IndexType glBoxMovementWindowIndex =
	{ { index[0] - static_cast<itk::IndexValueType>(glBoxMovementWindowRadius),
		index[1] - static_cast<itk::IndexValueType>(glBoxMovementWindowRadius) } };

	SizeType glBoxMovementWindowSize;
	glBoxMovementWindowSize.Fill(2*glBoxMovementWindowRadius+1);
	RegionType glBoxMovementWindow(glBoxMovementWindowIndex, glBoxMovementWindowSize);
	std::cout << "glBoxMovementWindow region:\n ";
	std::cout << glBoxMovementWindow << std::endl;

	itk::ConstNeighborhoodIterator<ImageType> iterator(glidingBoxRadius, image, glBoxMovementWindow);
	iterator.GoToBegin();
	
	size_t count = 0;
	while (!iterator.IsAtEnd())
	{
		std::cout << "Iterator region \n";
		std::cout << iterator.GetBoundingBoxAsImageRegion() << std::endl;
		PrintRegion(image, iterator.GetBoundingBoxAsImageRegion());
		double min = itk::NumericTraits<double>::max();
		double max = itk::NumericTraits<double>::min();
		//Find min and max value in gliding box
		SizeValueType numPixelsInGlidingBox = (2 * glBoxRadius + 1)*(2 * glBoxRadius + 1);
		for (SizeValueType i = 0; i < numPixelsInGlidingBox; ++i)
		{
			double val = static_cast<double>(iterator.GetPixel(i));
			if (val < min) { min = val; }
			if (val > max) { max = val; }
		}
		std::cout << "min = " << min << "\tmax = " << max << "\n";
		++iterator;
		++count;
		system("pause");
	}
	std::cout << "Number of box movements " << count << std::endl;
}