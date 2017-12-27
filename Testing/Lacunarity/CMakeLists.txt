cmake_minimum_required(VERSION 2.8)

set(NAME LacunarityTest)
#set(PROJECT ${NAME})

project(${NAME})

message(STATUS "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}\n")
#include_directories (${PROJECT_SOURCE_DIR})

add_executable(${NAME}  ${NAME}.cxx)
target_link_libraries(${NAME} ${OTB_LIBRARIES}) 

