cmake_minimum_required(VERSION 2.8)

set(NAME StokesMatrixImage)
project(${NAME})

add_executable(SinclairToStokesMatrixImage SinclairToStokesMatrixImage.cxx)
target_link_libraries(SinclairToStokesMatrixImage ${OTB_LIBRARIES}) 

add_executable(StokesMatrixImageMultilooking StokesMatrixImageMultilooking.cxx)
target_link_libraries(StokesMatrixImageMultilooking ${OTB_LIBRARIES}) 

