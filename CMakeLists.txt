cmake_minimum_required(VERSION 3.20)
project(botany LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
include(FetchContent)

FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG 1.0.1
)

FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.2
)

FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
)

FetchContent_Declare(
    glad
    GIT_REPOSITORY https://github.com/Dav1dde/glad.git
    GIT_TAG v2.0.6
    SOURCE_SUBDIR cmake
)

FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.8
)

FetchContent_MakeAvailable(glm catch2 glfw)

# glad
FetchContent_GetProperties(glad)
if(NOT glad_POPULATED)
    FetchContent_Populate(glad)
endif()
add_subdirectory("${glad_SOURCE_DIR}/cmake" glad_cmake)
glad_add_library(glad_gl REPRODUCIBLE API gl:core=4.1)

# imgui - build as a static library
FetchContent_GetProperties(imgui)
if(NOT imgui_POPULATED)
    FetchContent_Populate(imgui)
endif()
add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
target_link_libraries(imgui PUBLIC glfw glad_gl)

# Engine library (no graphics dependencies)
add_library(botany_engine STATIC
    src/engine/node.cpp
    src/engine/plant.cpp
    src/engine/hormone.cpp
    src/engine/meristem.cpp
    src/engine/sugar.cpp
    src/engine/engine.cpp
    src/engine/world_params.cpp
)
target_include_directories(botany_engine PUBLIC src)
target_link_libraries(botany_engine PUBLIC glm::glm)

# Serialization library
add_library(botany_serialization STATIC
    src/serialization/serializer.cpp
)
target_include_directories(botany_serialization PUBLIC src)
target_link_libraries(botany_serialization PUBLIC botany_engine)

# Renderer library
add_library(botany_renderer STATIC
    src/renderer/shader.cpp
    src/renderer/camera.cpp
    src/renderer/renderer.cpp
)
target_include_directories(botany_renderer PUBLIC src)
target_link_libraries(botany_renderer PUBLIC botany_engine glad_gl glfw glm::glm)

# Executables
add_executable(botany_realtime src/app_realtime.cpp)
target_link_libraries(botany_realtime PRIVATE botany_engine botany_renderer imgui)

add_executable(botany_headless src/app_headless.cpp)
target_link_libraries(botany_headless PRIVATE botany_engine botany_serialization)

add_executable(botany_playback src/app_playback.cpp)
target_link_libraries(botany_playback PRIVATE botany_engine botany_renderer botany_serialization imgui)

add_executable(botany_dump src/app_dump.cpp)
target_link_libraries(botany_dump PRIVATE botany_engine botany_serialization)

# Tests
enable_testing()
add_executable(botany_tests
    tests/test_genome.cpp
    tests/test_node.cpp
    tests/test_plant.cpp
    tests/test_hormone.cpp
    tests/test_meristem.cpp
    tests/test_engine.cpp
    tests/test_sugar.cpp
    tests/test_serializer.cpp
)
target_link_libraries(botany_tests PRIVATE botany_engine botany_serialization Catch2::Catch2WithMain)
include(CTest)
include(Catch)
catch_discover_tests(botany_tests)
