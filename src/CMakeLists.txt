cmake_minimum_required(VERSION 2.8)

set(NAME PolSig)

project(${NAME})

set(SOURCES FractalPolarizationSignature.cxx 
	SignaturesCalculationHelper.cxx
	)
include_directories (${PROJECT_SOURCE_DIR})

add_executable(FractalPolarizationSignature  ${SOURCES})
target_link_libraries(FractalPolarizationSignature ${OTB_LIBRARIES}) 

add_executable(SignatureOfLacunarity  SignatureOfLacunarity.cxx)
target_link_libraries(SignatureOfLacunarity ${OTB_LIBRARIES}) 


