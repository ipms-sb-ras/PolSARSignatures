cmake_minimum_required(VERSION 2.8)

set(NAME TestStokesMatrix)

project(${NAME})

add_executable(otbSinclairToStokesMatrixFunctor  otbSinclairToStokesMatrixFunctor.cxx)
target_link_libraries(otbSinclairToStokesMatrixFunctor ${OTB_LIBRARIES}) 

add_executable(PolarSynthesisFromStokesMatrix  PolarSynthesisFromStokesMatrix.cxx)
target_link_libraries(PolarSynthesisFromStokesMatrix ${OTB_LIBRARIES}) 

add_executable(otbSinclairToStokesMatrixImage  otbSinclairToStokesMatrixImage.cxx)
target_link_libraries(otbSinclairToStokesMatrixImage ${OTB_LIBRARIES}) 

add_executable(otbPolarimetricSynthesisStokesMatrixFilter  otbPolarimetricSynthesisStokesMatrixFilter.cxx)
target_link_libraries(otbPolarimetricSynthesisStokesMatrixFilter ${OTB_LIBRARIES}) 

