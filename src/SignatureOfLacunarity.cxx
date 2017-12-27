/*=========================================================================

Program:   SignatureOfLacunarity
Language:  C++
Date:      2017/08/16
Version:   0.1

Copyright (c) IPMS SB RAS

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.
=========================================================================*/

#include <ctype.h>
#include <vector>
#include <fstream>

#include "itkHistogram.h"
#include "itkListSample.h"
#include "itkSampleToHistogramFilter.h"
#include "itkLocalStatisticImageFunctionBase.hxx"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkRegularNeighborhood.hxx"

#include "otbImageFileReader.h"
#include "otbCommandLineArgumentParser.h"
#include "otbWrapperTypes.h"
#include "otbMultiChannelsPolarimetricSynthesisStokesFilter.h"

/*
 * Данная программа предназначена для построения сигнатуры лакунарности фрагмента радиолокационного изображения.
 * В качестве исходных используются изображения, полученные поляриметрическими радиолокаторами с синтезированной 
 * апертурой космического базирования, и прошедшие процедуру некогерентного накопления сигнала. Формирование сигнатуры 
 * основано на анализе лакунарности таких изображений при различных состояниях поляризационного эллипса.
 * Результат работы программы представляется в виде текстового файла с разделителями, в котором в первых двух столбцах 
 * записаны углы ориентации и эллиптичности поляризационного эллипса, в третьем и четвертом - рассчитанные для заданного 
 * фрагмента исходного изображения средние значения лакунарности и соответствующие им стандартные отклонения. 
 * Данный файл можно использовать для графического отображения поляризационной сигнатуры лакунарности во внешних программах, 
 * поддерживающих построение трехмерной графики. Программа предоставляет интерфейс командной строки, обеспечивающий настройку 
 * параметров расчета, задания расположения входных данных и каталога для записи результатов расчёта. Программа может быть 
 * использована для тематической обработки данных дистанционного зондирования земли с целью выявления неоднородностей 
 * пространственного распределения рассеивающих дискретных объектов, в частности, ветвей деревьев
*/

/**
* class LacunarityDongImageFunction
* Calculate the lacunarity over the neighborhood of a pixel
*
* Calculate the lacunarity in accordance with paper
* Dong P. Test of a new lacunarity estimation method for image texture
* analysis. Int. J. Remote Sensing. 2000, vol. 21, no. 17, 3369-3373
*/
namespace itk
{
	template< typename TInputImage, typename TCoordRep = float >
	class LacunarityDongImageFunction :
		public LocalStatisticImageFunctionBase < TInputImage, TCoordRep >
	{
	public:
		/** Standard class typedefs. */
		typedef LacunarityDongImageFunction Self;
		typedef LocalStatisticImageFunctionBase< TInputImage, TCoordRep > Baseclass;

		typedef typename Baseclass::InputImageType InputImageType;
		typedef typename Baseclass::RealType RealType;
		typedef typename Baseclass::IndexType IndexType;

		typedef itk::ImageRegionConstIterator<InputImageType> ImageRegionConstIteratorType;

		typedef typename InputImageType::RegionType RegionType;
		typedef typename RegionType::SizeType SizeType;
		typedef typename SizeType::SizeValueType SizeValueType;
		typedef typename InputImageType::PixelType  PixelType;
		typedef typename InputImageType::PointType  PointType;
		typedef typename std::vector<RegionType> VectorRegionType;

		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Run-time type information (and related methods). */
		itkTypeMacro(LacunarityDongImageFunction, LocalStatisticImageFunctionBase);

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Get/Set the radius of the window over which the
		lacunarity is evaluated */
		itkSetMacro(WindowRadius, unsigned int);
		itkGetConstReferenceMacro(WindowRadius, unsigned int);

		/** Get/Set the radius of the gliding box which moves through the window */
		itkGetConstReferenceMacro(GlidingBoxRadius, unsigned int);
		void SetGlidingBoxRadius(unsigned int radius = 1U)
		{
			m_GlidingBoxRadius = radius;
			if (m_GlidingBoxRadius > m_WindowRadius)
			{
				itkWarningMacro(<< "Gliding box radius is greater than window radius.\n"
					<< "Setting GlidingBoxRadius = WindowRadius.");
				m_GlidingBoxRadius = m_WindowRadius;
			}
		}

		/** Evaluate the function at specified index */
		virtual RealType EvaluateAtIndex(const IndexType & index) const ITK_OVERRIDE
		{
			if (!this->GetInputImage())			return (NumericTraits< RealType >::max());
			if (!this->IsInsideBuffer(index))	return (NumericTraits< RealType >::max());
			if (!this->m_NeighborhoodPointer)	itkExceptionMacro(<< "Neighborhood is not set!");

			// Gliding box will walks through a smaller region in order do not
			// cross window boundaries
			SizeType regionSize;
			regionSize.Fill(2 * (m_WindowRadius - m_GlidingBoxRadius) + 1);

			//Top left corner of the region to walk
			IndexType regionIndex =
			{ { index[0] - m_WindowRadius, index[1] - m_WindowRadius } };

			//Set a region for gliding box movements which centered at the index
			RegionType iteratingRegion(regionIndex, regionSize);

			// Check if it lies within buffered region
			if (!this->GetInputImage()->GetBufferedRegion().IsInside(iteratingRegion))
			{	// Normally it should not occurs
				itkExceptionMacro(<< "Window region is outside of buffered region!\n"
					<< "Window center location:\t" << index
					<< "\nWindow size:\t\t\t" << iteratingRegion.GetSize()
					<< "\nWindow corner:\t\t" << iteratingRegion.GetIndex() << "\n");
			}

			SizeType radius;
			radius.Fill(m_GlidingBoxRadius);

			ConstNeighborhoodIterator< TInputImage > nbIter(radius, this->GetInputImage(), iteratingRegion);

			MeasurementVectorType mv;
			while (!nbIter.IsAtEnd())
			{	//Find min and max value in gliding box
				auto minmax = std::minmax_element(nbIter.Begin(), nbIter.End());
				//Calculate the relative height of the column
				mv[0] = (*(*minmax.second) - *(*minmax.first)) / m_GlidingBoxRadius - 1.0;
				m_Sample->PushBack(mv);
				++nbIter;
			}
			double lacunarity = ComputeLacunarity();
			return static_cast<RealType>(lacunarity);
		}

		void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE
		{
			Baseclass::PrintSelf(os, indent);
			os << indent << "Window radius (size): " << this->m_WindowRadius << " ("
				<< 2 * this->m_WindowRadius + 1 << "x" << 2 * this->m_WindowRadius + 1 << ")" << std::endl;
			os << indent << "Gliding box radius (box size): " << this->m_GlidingBoxRadius << " ("
				<< 2 * this->m_GlidingBoxRadius + 1 << "x" << 2 * this->m_GlidingBoxRadius + 1 << ")" << std::endl;
		}


	protected:
		typedef itk::Vector<double, 1> MeasurementVectorType;
		typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;

		typedef itk::Statistics::Histogram< float,
			itk::Statistics::DenseFrequencyContainer2 > HistogramType;

		typedef itk::Statistics::SampleToHistogramFilter<SampleType, HistogramType>
			SampleToHistogramFilterType;

		LacunarityDongImageFunction() :
			Baseclass(),
			m_WindowRadius(10U),
			m_GlidingBoxRadius(1U)
		{
			m_Sample = SampleType::New();
		};

		~LacunarityDongImageFunction() {};

		virtual double ComputeLacunarity() const
		{
			SampleToHistogramFilterType::Pointer sampleToHistogramFilter =
				SampleToHistogramFilterType::New();
			sampleToHistogramFilter->SetInput(m_Sample);

			auto sampleSize = m_Sample->Size();
			SampleToHistogramFilterType::HistogramSizeType histogramSize(1);
			histogramSize.Fill(sampleSize);
			sampleToHistogramFilter->SetHistogramSize(histogramSize);
			sampleToHistogramFilter->Update();

			const HistogramType* histogram = sampleToHistogramFilter->GetOutput();

			double secondMoment = 0.0;
			double firstMoment = 0.0;
			for (unsigned int i = 0; i < histogram->GetSize()[0]; i++)
			{
				double probability = static_cast<double>(histogram->GetFrequency(i)) / sampleSize;
				double mass = histogram->GetMeasurement(i, 0);
				secondMoment += (mass*mass*probability);
				firstMoment += (mass*probability);
			}

			double lacunarity = secondMoment / (firstMoment*firstMoment);

			if (lacunarity < 1e-7)
				lacunarity = 0.0;
			m_Sample->Clear();
			return lacunarity;
		}

	private:
		LacunarityDongImageFunction(const Self &); //purposely not implemented
		void operator=(const Self &);    //purposely not implemented

		unsigned int m_WindowRadius;
		unsigned int m_GlidingBoxRadius;
		SampleType::Pointer m_Sample;
	};
} // end namespace itk


typedef otb::CommandLineArgumentParser ParserType;
typedef otb::CommandLineArgumentParseResult ParserResultType;

std::string GetOptionFromInput(const std::string option,
	const std::string synonym, ParserResultType* parserResult)
{
	std::string optString;
	if (parserResult->IsOptionPresent(option)) {
		optString = parserResult->GetParameterString(option);
	}
	else {
		if (parserResult->IsOptionPresent(synonym)) {
			optString = parserResult->GetParameterString(synonym);
		}
	}
	return optString;
}

typedef enum {
	CO_POLAR,
	CROSS_POLAR,
} PolarizationMode;

std::string GetInputImagePixelType(const std::string& filename) {

	std::string fname(filename);
	if (itksys::SystemTools::GetFilenameLastExtension(filename) == ".hdr") { // ENVI file. 
		fname = itksys::SystemTools::GetFilenamePath(filename) + '/'
			+ itksys::SystemTools::GetFilenameWithoutLastExtension(filename);
	}
	otb::ImageIOBase::Pointer imageIO =
		otb::ImageIOFactory::CreateImageIO(fname.c_str(),
			otb::ImageIOFactory::ReadMode);

	if (!imageIO) {
		itkGenericExceptionMacro(<< "Unknown format of the input image:\t" << fname);
	}
	imageIO->SetFileName(fname);
	imageIO->ReadImageInformation();
	return imageIO->GetPixelTypeAsString(imageIO->GetPixelType());
}

typedef itk::Index<2> IndexType;
typedef std::vector<IndexType> IndexVectorType;

size_t ConvertMaskImageToPOIList(const std::string& inputMaskImage, IndexVectorType& poiVector)
{
	typedef otb::Wrapper::UInt8ImageType MaskImageType;
	/* Открываем изображение-маску. */
	otb::ImageFileReader<MaskImageType>::Pointer readerMask =
		otb::ImageFileReader<MaskImageType>::New();
	readerMask->SetFileName(inputMaskImage);
	try
	{
		readerMask->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << "Cannot open mask image!\n";
		std::cerr << e.what() << std::endl;
		exit(-1);
	}
	MaskImageType::Pointer maskImage = readerMask->GetOutput();
	itk::ImageRegionConstIteratorWithIndex<MaskImageType> imageIterator(readerMask->GetOutput(),
		readerMask->GetOutput()->GetLargestPossibleRegion());
	while (!imageIterator.IsAtEnd()) {
		if (imageIterator.Get()) {
			poiVector.push_back(imageIterator.GetIndex());
		}
		++imageIterator;
	}
	return poiVector.size();
}

int main(int argc, char* argv[])
{
	typedef otb::Wrapper::ComplexFloatVectorImageType InputImageType;
	typedef otb::Wrapper::FloatImageType FloatImageType;

	/* Инициализация входных параметров */
	ParserType::Pointer parser = ParserType::New();
	parser->AddOption("--InputFile", "Input filename.", "-i", 1, true);
	parser->AddOption("--MaskImage", "Mask image filename.", "-m", 1, true);
	parser->AddOption("--WindowRadius", "Radius of image window. Default value is 10 (window size 21x21)", "-w", 1, false);
	parser->AddOption("--GlidingBoxRadius", "Radius of gliding box. Default value is 1 (box size 3x3)", "-r", 1, false);
	parser->AddOption("--PolMode", "Polarization mode. Default value is \"copolar\"", "-p", 1, false);
	parser->AddOption("--Step", "Polarization ellipse angles incremental step. Default value is 3", "-s", 1, false);
	parser->AddOption("--OutputDir", "Output directory name for results. Default value is a directory where input file is located", "-o", 1, false);

	/* Разбор командной строки */
	ParserResultType::Pointer parserResult = ParserResultType::New();
	try
	{
		parser->ParseCommandLine(argc, argv, parserResult);
	}
	catch (CommandLineArgumentParserArgumentErrorException& ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

	std::string inputImageFile = GetOptionFromInput("--InputFile", "-i", parserResult);
	if (!itksys::SystemTools::FileExists(inputImageFile.c_str(), true))
	{
		std::cerr << "File " << inputImageFile << " does not exist!" << std::flush;
		return EXIT_FAILURE;
	}

	std::string inputMaskImage = GetOptionFromInput("--MaskImage", "-m", parserResult);
	if (!itksys::SystemTools::FileExists(inputMaskImage.c_str(), true))
	{
		std::cerr << "File " << inputMaskImage << " does not exist!" << std::flush;
		return EXIT_FAILURE;
	}

	size_t windowRadius = 10;
	if (!GetOptionFromInput("--WindowRadius", "-w", parserResult).empty())
	{
		std::istringstream s(GetOptionFromInput("--WindowRadius", "-w", parserResult));
		if (isdigit(s.str()[0])) {
			s >> windowRadius;
		}
		else {
			std::cout << "Warning! Radius is not represented as digit. ";
			std::cout << "Setting default value: 10.\n";
		}
	}

	size_t glidingBoxRadius = 1;
	if (!GetOptionFromInput("--GlidingBoxRadius", "-r", parserResult).empty())
	{
		std::istringstream s(GetOptionFromInput("--GlidingBoxRadius", "-r", parserResult));
		if (isdigit(s.str()[0])) {
			s >> glidingBoxRadius;
		}
		else {
			std::cout << "Warning! Gliding box radius is not represented as digit. ";
			std::cout << "Setting default value: 1.\n";
		}
	}

	PolarizationMode polMode = PolarizationMode::CO_POLAR;
	std::string modeName = "copolar";
	if (!GetOptionFromInput("--PolMode", "-p", parserResult).empty())
	{
		std::string polModeString = GetOptionFromInput("--PolMode", "-p", parserResult);
		if (polModeString == "crosspolar") {
			polMode = PolarizationMode::CROSS_POLAR;
			modeName = "crosspolar";
		}
	}

	float step = 3.0;
	if (!GetOptionFromInput("--WindowRadius", "-w", parserResult).empty())
	{
		std::istringstream s(GetOptionFromInput("--Step", "-s", parserResult));
		if (isdigit(s.str()[0])) {
			s >> step;
		}
		else {
			std::cout << "Warning! Angles incrmental step is not represented as digit. ";
			std::cout << "Setting default value: 3.\n";
		}
	}

	std::string outputDir = GetOptionFromInput("--OutputDir", "-o", parserResult);
	if (!outputDir.empty())
	{
		if (!itksys::SystemTools::FileExists(outputDir.c_str(), false))
		{
			itksys::SystemTools::MakeDirectory(outputDir.c_str());
		}
	}
	/* Проверяем тип входного изображения. Оно должно содержать 4 канала, по одному
	 * для каждой из поляризаций HH, HV, VH, VV, и иметь комплексные значения пикселей
	 */
	if (GetInputImagePixelType(inputImageFile).compare("complex") != 0)
	{
		std::cerr << "Input image must consist of complex type pixels!\n";
		return EXIT_FAILURE;
	}
	/* Открываем файл изображения. */
	otb::ImageFileReader<InputImageType>::Pointer reader =
		otb::ImageFileReader<InputImageType>::New();
	reader->SetFileName(inputImageFile);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << "Cannot open image!\n";
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	InputImageType::Pointer inputImage = reader->GetOutput();
	if (inputImage->GetNumberOfComponentsPerPixel() != 4)
	{
		std::cerr << "Number of image channels must be equal to 4 (HH, HV, VH, VV)!\n";
		return EXIT_FAILURE;
	}
	/* Считываем изображение-маску в массив точек,
	 * для которых будет производиться расчёт. */
	IndexVectorType poiList;
	if (!ConvertMaskImageToPOIList(inputMaskImage, poiList)) {
		std::cerr << "Nothing to compute!\n";
		return EXIT_FAILURE;
	}

	/* Создание файла для записи результатов расчета */
	std::string outputBaseFilename;
	if (outputDir.empty())
	{
		outputBaseFilename = itksys::SystemTools::GetFilenamePath(inputImageFile);
	}
	else {
		outputBaseFilename = outputDir;
	}
	outputBaseFilename.append("/");
	outputBaseFilename.append(itksys::SystemTools::GetFilenameWithoutLastExtension(inputImageFile));

	std::ostringstream oss;
	oss << outputBaseFilename << "_" << modeName << "_w"
		<< windowRadius << "_r" << glidingBoxRadius << "_Lacunarity.avrg";

	std::ofstream outFileStream(oss.str());
	outFileStream.width(12);
	outFileStream.setf(std::ios_base::fixed, std::ios_base::floatfield);

	std::ostringstream ossInfo;
	ossInfo << "#Input image filename: " << inputImageFile
		<< "\n#Input mask image filename: " << inputMaskImage
		<< "\n#Image window radius is: " << windowRadius
		<< "\n#Gliding box radius is: " << glidingBoxRadius
		<< "\n#psi\tkhi\tLacunarity\tStdDev";
	outFileStream << ossInfo.str() << std::endl;

	/* Углы эллиптичности и ориентации поляризационного эллипса*/
	float psiStart, psiStop;
	float khiStart, khiStop;
	psiStart = 0.;
	psiStop = 180.;
	khiStart = -45.;
	khiStop = 45.;

	typedef otb::MultiChannelsPolarimetricSynthesisStokesFilter<InputImageType,
		FloatImageType> PolarimetricSynthesisFilterType;
	typedef itk::LacunarityDongImageFunction< FloatImageType, float >  LacunarityDongFunctionType;
	typedef itk::RegularNeighborhood<float, 2> RegularNeighborhoodType;

	std::shared_ptr<RegularNeighborhoodType> neighborhood = std::make_shared<RegularNeighborhoodType>();
	neighborhood->SetRadius(glidingBoxRadius);

	LacunarityDongFunctionType::Pointer function = LacunarityDongFunctionType::New();
	function->SetNeighborhood(neighborhood);
	function->SetWindowRadius(windowRadius);
	function->SetGlidingBoxRadius(glidingBoxRadius);
	/* Основной расчет */
	for (float psi = psiStart; psi <= psiStop; psi += step)
	{
		for (float khi = khiStart; khi <= khiStop; khi += step)
		{
			PolarimetricSynthesisFilterType::Pointer
				polarimetricSynthesisFilter = PolarimetricSynthesisFilterType::New();
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
			}

			polarimetricSynthesisFilter->Update();
			FloatImageType::Pointer synthesizedImage =
				polarimetricSynthesisFilter->GetOutput();

			function->SetInputImage(synthesizedImage);
			double valSum = 0.;
			double valSquaredSum = 0.;
			size_t counter = 0;
			float val;
			for (auto index : poiList) {
				val = function->EvaluateAtIndex(index);
				if (!std::isnan(val)) {
					valSum += val;
					valSquaredSum += (val*val);
					++counter;
				}
			}
			float n = static_cast<float>(counter);
			float mean = valSum / n;
			float stdDev = std::sqrt((valSquaredSum - n * mean * mean) / (n - 1.));

			outFileStream << psi << "\t" << khi << "\t" << mean << "\t" << stdDev << std::endl;
		}
	}
	return EXIT_SUCCESS;
}
