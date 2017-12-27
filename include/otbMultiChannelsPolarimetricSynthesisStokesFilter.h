/*=========================================================================

  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See OTBCopyright.txt for details.


  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

  =========================================================================*/
#ifndef __otbMultiChannelsPolarimetricSynthesisStokesFilter_h
#define __otbMultiChannelsPolarimetricSynthesisStokesFilter_h

#include "itkInPlaceImageFilter.h"
#include "otbPolarimetricData.h"
#include "otbSinclairToStokesMatrixFunctor.h"
#include <complex>

namespace otb {

    /** \class MultiChannelsPolarimetricSynthesisStokesFilter
     * \brief This class computes the polarimetric synthesis from two to four radar images, depening on the polarimetric architecture.
     *
     * It has the same behaviour as the PolarimetricSynthesisImageFilter expect the fact that it
     * considers a VectorImage as input which each channels is HH, HV, VH and VV (in this particular order).
     *
     * \ingroup SARPolarimetry
     * \sa PolarimetricSynthesisFilter
     * \sa PolarimetricSynthesisFunctor
     *
     */

    template <class TInputImage, class TOutputImage>
    class ITK_EXPORT MultiChannelsPolarimetricSynthesisStokesFilter : 
        public itk::InPlaceImageFilter < TInputImage, TOutputImage > 
    {
    public:
        /** Standard class typedefs. */
        typedef MultiChannelsPolarimetricSynthesisStokesFilter Self;
        typedef itk::InPlaceImageFilter<TInputImage, TOutputImage> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(MultiChannelsPolarimetricSynthesisStokesFilter, InPlaceImageFilter);

        /** Some typedefs. */
        typedef TInputImage                             InputImageType;
        //typedef typename InputImageType::InputPixelType InputPixelType;
        typedef typename InputImageType::ConstPointer   InputImagePointer;
        typedef typename InputImageType::RegionType     InputImageRegionType;
        typedef typename InputImageType::InternalPixelType InputImagePixelType;
        typedef TOutputImage                            OutputImageType;
        typedef typename OutputImageType::Pointer       OutputImagePointer;
        typedef typename OutputImageType::RegionType    OutputImageRegionType;
        typedef typename OutputImageType::PixelType     OutputImagePixelType;
        typedef typename std::complex <double>          ComplexType;

        //typedef typename itk::FixedArray<ComplexType, 2> ComplexArrayType;
        //typedef typename itk::FixedArray<int, 4>        IndexArrayType;

        /** Set/Get PsiI */
        itkSetMacro(PsiI, double);
        itkGetMacro(PsiI, double);
        /** Set/Get KhiI */
        itkSetMacro(KhiI, double);
        itkGetMacro(KhiI, double);
        /** Set/Get PsiR */
        itkSetMacro(PsiR, double);
        itkGetMacro(PsiR, double);
        /** Set/Get KhiR */
        itkSetMacro(KhiR, double);
        itkGetMacro(KhiR, double);
        /** Set/Get EmissionH */
        itkSetMacro(EmissionH, bool);
        itkGetMacro(EmissionH, bool);
        /** Set/Get EmissionV */
        itkSetMacro(EmissionV, bool);
        itkGetMacro(EmissionV, bool);
        /** Set/Get Mode */
        itkSetMacro(Mode, int);
        itkGetMacro(Mode, int);
        /** Set the gain */
        itkSetMacro(Gain, double);
        itkGetMacro(Gain, double);
        /** Force the copolar mode */
        void ForceCoPolar();
        /** Force the crosspolar mode */
        void ForceCrossPolar();

    protected:
        /** Constructor */
        MultiChannelsPolarimetricSynthesisStokesFilter();

        /** Destructor */
        virtual ~MultiChannelsPolarimetricSynthesisStokesFilter() {
        }

        /** MultiChannelsPolarimetricSynthesisStokesFilter can produce an image
         * which is a synthesis of channels HH, HV, VH and VV.
         *
         * As such, MultiChannelsPolarimetricSynthesisStokesFilter
         * needs to provide an implementation for
         * GenerateOutputInformation() in order to inform the pipeline
         * execution model.  The original documentation of this method is
         * below.
         *
         * \sa ProcessObject::GenerateOutputInformaton()  */
        virtual void GenerateOutputInformation();

        virtual void BeforeThreadedGenerateData();

        /** MultiChannelsPolarimetricSynthesisStokesFilter can be implemented as a multithreaded filter.
         * Therefore, this implementation provides a ThreadedGenerateData() routine
         * which is called for each processing thread. The output image data is
         * allocated automatically by the superclass prior to calling
         * ThreadedGenerateData().  ThreadedGenerateData can only write to the
         * portion of the output image specified by the parameter
         * "outputRegionForThread"
         *
         * \sa ImageToImageFilter::ThreadedGenerateData(),
         *     ImageToImageFilter::GenerateData()  */
        void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                itk::ThreadIdType threadId);


        /** Verify and force the inputs, if only  2 or 3 channels are present */
        void VerifyAndForceInputs();

        void PrintSelf(std::ostream& os, itk::Indent indent) const;
        
        /** */
        inline double CalculatePolarimetricResponse(const InputImagePixelType Shh,
                const InputImagePixelType Shv, const InputImagePixelType Svh, 
                const InputImagePixelType Svv);
        

    private:
        MultiChannelsPolarimetricSynthesisStokesFilter(const Self &); //purposely not implemented

        /** Psi Incident */
        double m_PsiI;
        /** Khi Incident */
        double m_KhiI;
        /** Psi Reflected */
        double m_PsiR;
        /** Khi Reflected */
        double m_KhiR;

        /** Gain */
        double m_Gain;

        /** None = 0 , copolar = 1 , crosspolar = 2 */
        int m_Mode;

        /** Architecture Type */
        PolarimetricData::Pointer m_ArchitectureType;

        /** Emission mode */
        bool m_EmissionH;
        bool m_EmissionV;

    };

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbMultiChannelsPolarimetricSynthesisStokesFilter.txx"
#endif

#endif
