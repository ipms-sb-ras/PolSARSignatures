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
#ifndef otbPolarimetricSynthesisStokesMatrixFilter_h
#define otbPolarimetricSynthesisStokesMatrixFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkVariableLengthVector.h"

namespace otb {

    /** \class PolarimetricSynthesisStokesMatrixFilter
     * \brief This class computes the polarimetric synthesis from image where pixels 
     *  are the elements of Stokes scattering operator (matrix).
     *
     *
     * \ingroup SARPolarimetry
     * \sa PolarimetricSynthesisFilter
     *
     */

    template <class TInputImage, class TOutputImage>
    class ITK_EXPORT PolarimetricSynthesisStokesMatrixFilter : 
        public itk::InPlaceImageFilter < TInputImage, TOutputImage > 
    {
    public:
        /** Standard class typedefs. */
        typedef PolarimetricSynthesisStokesMatrixFilter Self;
        typedef itk::InPlaceImageFilter<TInputImage, TOutputImage> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(PolarimetricSynthesisStokesMatrixFilter, InPlaceImageFilter);

        /** Some typedefs. */
        typedef TInputImage                             InputImageType;
        typedef typename InputImageType::ConstPointer   InputImagePointer;
        typedef typename InputImageType::RegionType     InputImageRegionType;
        typedef typename InputImageType::PixelType      InputImagePixelType;
        typedef TOutputImage                            OutputImageType;
        typedef typename OutputImageType::Pointer       OutputImagePointer;
        typedef typename OutputImageType::RegionType    OutputImageRegionType;
        typedef typename OutputImageType::PixelType     OutputImagePixelType;
        typedef double                                  RealType;
        //typedef typename itk::VariableLengthVector<RealType> PixelType;


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
        /** Set/Get Mode */
        itkSetMacro(Mode, int);
        itkGetMacro(Mode, int);

        /** Force the copolar mode */
        void ForceCoPolar();
        /** Force the crosspolar mode */
        void ForceCrossPolar();

    protected:
        /** Constructor */
        PolarimetricSynthesisStokesMatrixFilter();

        /** Destructor */
        virtual ~PolarimetricSynthesisStokesMatrixFilter() {
        }

        virtual void GenerateOutputInformation();

        virtual void BeforeThreadedGenerateData();

        void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                itk::ThreadIdType threadId);

        void PrintSelf(std::ostream& os, itk::Indent indent) const;
        
        /** */
        inline double CalculatePolarimetricResponse(const InputImagePixelType& m);
        

    private:
        PolarimetricSynthesisStokesMatrixFilter(const Self &); //purposely not implemented

        /** Psi Incident */
        double m_PsiI;
        /** Khi Incident */
        double m_KhiI;
        /** Psi Reflected */
        double m_PsiR;
        /** Khi Reflected */
        double m_KhiR;

        /** None = 0 , copolar = 1 , crosspolar = 2 */
        int m_Mode;

    };

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbPolarimetricSynthesisStokesMatrixFilter.txx"
#endif

#endif // otbPolarimetricSynthesisStokesMatrixFilter_h
