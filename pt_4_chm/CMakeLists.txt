cmake_minimum_required(VERSION 3.10)
project(IntegrationSchemes)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(integration
    "main.cpp"
    "point.cpp"
    "integration_scheme_interval.cpp"
)
target_include_directories(integration PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})