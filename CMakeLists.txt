cmake_minimum_required(VERSION 3.14)
project(SurfaceExt)

include(FetchContent)

#GLM
FetchContent_Declare(glm URL https://github.com/g-truc/glm/archive/refs/tags/0.9.9.8.tar.gz)
FetchContent_MakeAvailable(glm)

#SurfaceExt executable
add_executable(SurfaceExt)
set_target_properties(SurfaceExt PROPERTIES
    CXX_STANDARD 20
)
target_link_libraries(SurfaceExt PRIVATE
    glm
)
target_include_directories(SurfaceExt PRIVATE .)
add_custom_command(TARGET SurfaceExt POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/test_data/ .
)

#Add subdirectories
add_subdirectory(SurfaceExtractor)