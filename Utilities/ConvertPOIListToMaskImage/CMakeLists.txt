# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8)

project(ConvertPOIListToMaskImage)


add_executable(ConvertPOIListToMaskImage  ConvertPOIListToMaskImage.cxx)
target_link_libraries(ConvertPOIListToMaskImage ${OTB_LIBRARIES}) 
