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
#ifndef __otbMultiChannelsPolarimetricSynthesisStokesFilter_txx
#define __otbMultiChannelsPolarimetricSynthesisStokesFilter_txx

#include <complex>

#include "itkVariableLengthVector.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "otbMath.h"
#include "otbMultiChannelsPolarimetricSynthesisStokesFilter.h"
#include "otbSinclairToStokesMatrixFunctor.h"

namespace otb {

    /**
     * Constructor
     */
    template <class TInputImage, class TOutputImage>
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::MultiChannelsPolarimetricSynthesisStokesFilter() {
        this->SetNumberOfRequiredInputs(1);
        this->InPlaceOff();
        SetEmissionH(false);
        SetEmissionV(false);
        SetGain(1);
        m_ArchitectureType = PolarimetricData::New();
    }

    /**
     * GenerateOutputInformation()
     */
    template <class TInputImage, class TOutputImage>
    void
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::GenerateOutputInformation() {
        // do not call the superclass' implementation of this method since
        // this filter allows the input the output to be of different dimensions

        // get pointers to the input and output
        typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
        typename Superclass::InputImageConstPointer inputPtr = this->GetInput();

        if (!outputPtr || !inputPtr) {
            return;
        }

        // Set the output image largest possible region.  Use a RegionCopier
        // so that the input and output images can be different dimensions.
        OutputImageRegionType outputLargestPossibleRegion;
        this->CallCopyInputRegionToOutputRegion(outputLargestPossibleRegion,
                inputPtr->GetLargestPossibleRegion());
        outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);

        // Set the output spacing and origin
        const itk::ImageBase<Superclass::InputImageDimension> *phyData;

        phyData
                = dynamic_cast<const itk::ImageBase<Superclass::InputImageDimension>*> (this->GetInput());

        if (phyData) {
            // Copy what we can from the image from spacing and origin of the input
            // This logic needs to be augmented with logic that select which
            // dimensions to copy
            unsigned int i, j;
            const typename InputImageType::SpacingType&
                    inputSpacing = inputPtr->GetSpacing();
            const typename InputImageType::PointType&
                    inputOrigin = inputPtr->GetOrigin();
            const typename InputImageType::DirectionType&
                    inputDirection = inputPtr->GetDirection();

            typename OutputImageType::SpacingType outputSpacing;
            typename OutputImageType::PointType outputOrigin;
            typename OutputImageType::DirectionType outputDirection;

            // copy the input to the output and fill the rest of the
            // output with zeros.
            for (i = 0; i < Superclass::InputImageDimension; ++i) {
                outputSpacing[i] = inputSpacing[i];
                outputOrigin[i] = inputOrigin[i];
                for (j = 0; j < Superclass::OutputImageDimension; ++j) {
                    if (j < Superclass::InputImageDimension) {
                        outputDirection[j][i] = inputDirection[j][i];
                    } else {
                        outputDirection[j][i] = 0.0;
                    }
                }
            }
            for (; i < Superclass::OutputImageDimension; ++i) {
                outputSpacing[i] = 1.0;
                outputOrigin[i] = 0.0;
                for (j = 0; j < Superclass::OutputImageDimension; ++j) {
                    if (j == i) {
                        outputDirection[j][i] = 1.0;
                    } else {
                        outputDirection[j][i] = 0.0;
                    }
                }
            }

            // set the spacing and origin
            outputPtr->SetSpacing(outputSpacing);
            outputPtr->SetOrigin(outputOrigin);
            outputPtr->SetDirection(outputDirection);

        } else {
            // pointer could not be cast back down
            itkExceptionMacro( << "otb::MultiChannelsPolarimetricSynthesisStokesFilter::GenerateOutputInformation "
                    << "cannot cast input to "
                    << typeid (itk::ImageBase<Superclass::InputImageDimension>*).name());
        }
    }

    /**
     * ThreadedGenerateData Performs the pixel-wise addition
     */
    template <class TInputImage, class TOutputImage>
    void
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
            itk::ThreadIdType threadId) {

        InputImagePointer inputPtr = this->GetInput();
        OutputImagePointer outputPtr = this->GetOutput(0);

        // Define the portion of the input to walk for this thread, using
        // the CallCopyOutputRegionToInputRegion method allows for the input
        // and output images to be different dimensions
        InputImageRegionType inputRegionForThread;
        this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);

        // Define the iterators
        itk::ImageRegionConstIterator<TInputImage> inputIt(inputPtr, inputRegionForThread);
        itk::ImageRegionIterator<TOutputImage> outputIt(outputPtr, outputRegionForThread);

        itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

        inputIt.GoToBegin();
        outputIt.GoToBegin();

        ArchitectureType val = m_ArchitectureType->GetArchitectureType();

        // Computation with 4 channels
        switch (val) {
            case HH_HV_VH_VV:
                while (!inputIt.IsAtEnd()) {
                    outputIt.Set( static_cast<OutputImagePixelType>(
                        m_Gain * CalculatePolarimetricResponse(inputIt.Get()[0], 
                        inputIt.Get()[1], inputIt.Get()[2], inputIt.Get()[3])));
                    ++inputIt;
                    ++outputIt;
                    progress.CompletedPixel(); // potential exception thrown here
                }
                break;

                // With 3 channels : HH HV VV ou HH VH VV
            case HH_HV_VV:
                while (!inputIt.IsAtEnd()) {
                    outputIt.Set( static_cast<OutputImagePixelType>(
                        m_Gain * CalculatePolarimetricResponse(inputIt.Get()[0], 
                        inputIt.Get()[1],inputIt.Get()[1], inputIt.Get()[2])));
                    ++inputIt;
                    ++outputIt;
                    progress.CompletedPixel(); // potential exception thrown here
                }
                break;

                // Only HH and HV are present
            case HH_HV:
                while (!inputIt.IsAtEnd()) {
                    outputIt.Set( static_cast<OutputImagePixelType>(
                        m_Gain * CalculatePolarimetricResponse(inputIt.Get()[0], 
                        inputIt.Get()[1], 0, 0)));
                    ++inputIt;
                    ++outputIt;
                    progress.CompletedPixel(); // potential exception thrown here
                }
                break;

                // Only VH and VV are present
            case VH_VV:
                while (!inputIt.IsAtEnd()) {
                    /** CHANGED: � �������� ���� ��������������, ��� ������� �����������
                                     * ����� 4 ������ � ���  VH_VV ���������� ��������� 2 ������.
                                     *  �� ����� ���� �� ����� ����� ����������� � 2 ��������, �������
                                     *  ������� ��� ������� �������� ������� �� ������ ������ ���������� � 0.
                                    */
                    outputIt.Set(static_cast<OutputImagePixelType>(
                        m_Gain * CalculatePolarimetricResponse(0, 0, 
                        inputIt.Get()[0/*2*/], inputIt.Get()[1/*3*/])));
                    ++inputIt;
                    ++outputIt;
                    progress.CompletedPixel(); // potential exception thrown here
                }
                break;

            default:
                itkExceptionMacro("Unknown architecture : Polarimetric synthesis is impossible !");
                return;
        }

    }

    /**
     * Verify and force the inputs, if only  2 or 3 channels are present
     */
    template <class TInputImage, class TOutputImage>
    void
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::VerifyAndForceInputs() {

        ArchitectureType val = m_ArchitectureType->GetArchitectureType();

        switch (val) {

            case HH_HV_VH_VV:
                //std::cout << "HH_HV_VH_VV polarimetric architecture." << std::endl;
                break;
            case HH_HV_VV:
                //std::cout << "HH_HV_VV polarimetric architecture." << std::endl;
                break;
            case HH_VH_VV:
                //std::cout << "HH_VH_VV polarimetric architecture." << std::endl;
                break;
                // Only HH and HV are present
            case HH_HV:
                //std::cout << "HH_HV polarimetric architecture." << std::endl;
                // Forcing KhiI=0 PsiI=0
                this->SetKhiI(0);
                this->SetPsiI(0);
                break;

                // Only VH and VV are present
            case VH_VV:
                //std::cout << "VH_VV polarimetric architecture." << std::endl;
                // Forcing KhiI=0 PsiI=90
                this->SetKhiI(0);
                this->SetPsiI(90);
                break;

            default:
                itkExceptionMacro("Unknown architecture : Polarimetric synthesis is impossible !!");
                return;
        }

        if (GetMode() == 1) ForceCoPolar();
        else if (GetMode() == 2) ForceCrossPolar();

    }

    /**
     * BeforeThreadedGenerateData
     */
    template <class TInputImage, class TOutputImage>
    void
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::BeforeThreadedGenerateData() {

        int NumberOfImages = this->GetInput()->GetNumberOfComponentsPerPixel();

        // First Part. Determine the kind of architecture of the input picture
        m_ArchitectureType->DetermineArchitecture(NumberOfImages, GetEmissionH(), GetEmissionV());

        // Second Part. Verify and force the inputs
        VerifyAndForceInputs();
    }

    /**
     * Force Copolar mode
     */
    template <class TInputImage, class TOutputImage>
    void
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::ForceCoPolar() {
        SetPsiR(m_PsiI);
        SetKhiR(m_KhiI);
    }

    /**
     * Force Crosspolar mode
     */
    template <class TInputImage, class TOutputImage>
    void
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::ForceCrossPolar() {
        SetPsiR(m_PsiI + 90);
        SetKhiR(-m_KhiI);
        SetMode(2);
    }

    /**
     * Printself
     */
    template <class TInputImage, class TOutputImage>
    void
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::PrintSelf(std::ostream& os, itk::Indent indent) const {
        this->Superclass::PrintSelf(os, indent);
        os << indent << "PsiI: " << m_PsiI << std::endl;
        os << indent << "KhiI: " << m_KhiI << std::endl;
        os << indent << "PsiR: " << m_PsiR << std::endl;
        os << indent << "KhiR: " << m_KhiR << std::endl;
        os << indent << "Gain: " << m_Gain << std::endl;
    }

    template <class TInputImage, class TOutputImage>
    inline double
    MultiChannelsPolarimetricSynthesisStokesFilter<TInputImage, TOutputImage>
    ::CalculatePolarimetricResponse(const InputImagePixelType Shh,
            const InputImagePixelType Shv, const InputImagePixelType Svh, 
            const InputImagePixelType Svv) {
        
        typedef otb::Functor::SinclairToStokesMatrixFunctor<InputImagePixelType, 
                InputImagePixelType, InputImagePixelType, InputImagePixelType, 
                itk::VariableLengthVector<OutputImagePixelType> >  StokesFunctorType;


        // Define the incident normalized Stokes vector
        const double si0 = 1.0;
        const double si1 = cos(2.*m_PsiI*CONST_PI_180) * cos(2.*m_KhiI*CONST_PI_180);
        const double si2 = sin(2.*m_PsiI*CONST_PI_180) * cos(2.*m_KhiI*CONST_PI_180);
        const double si3 = sin(2.*m_KhiI*CONST_PI_180);
        // Define the recieved normalized Stokes vector
        const double sr0 = 1.0;
        const double sr1 = cos(2.*m_PsiR*CONST_PI_180) * cos(2.*m_KhiR*CONST_PI_180);
        const double sr2 = sin(2.*m_PsiR*CONST_PI_180) * cos(2.*m_KhiR*CONST_PI_180);
        const double sr3 = sin(2.*m_KhiR*CONST_PI_180);

        // Calculating polarimetric response
        StokesFunctorType functor;
        auto m = functor.operator ()(Shh, Shv, Svh, Svv);
        double sigma;
        sigma  = sr0 * ( m[0]*si0 + m[1]*si1 + m[2]*si2 + m[3]*si3 );
        sigma += sr1 * ( m[4]*si0 + m[5]*si1 + m[6]*si2 + m[7]*si3 );
        sigma += sr2 * ( m[8]*si0 + m[9]*si1 + m[10]*si2 + m[11]*si3 );
        sigma += sr3 * ( m[12]*si0 + m[13]*si1 + m[14]*si2 + m[15]*si3 );
        //sigma *= (4. * CONST_PI); 

        return (sigma);    
    }

} // end namespace otb

#endif
