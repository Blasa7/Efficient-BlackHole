cmake_minimum_required(VERSION 3.16)

project(Networks LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

# Gather all files recursively
file(GLOB_RECURSE PROJECT_FILES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/*.cpp"
    "${CMAKE_SOURCE_DIR}/*.h"
    "${CMAKE_SOURCE_DIR}/*.hpp"
)

# Mirror actual filesystem layout in Visual Studio
source_group(
    TREE ${CMAKE_SOURCE_DIR}
    FILES ${PROJECT_FILES}
)

add_executable(${PROJECT_NAME} ${PROJECT_FILES})

target_include_directories(${PROJECT_NAME}
    PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
    OpenMP::OpenMP_CXX
)