cmake_minimum_required(VERSION 3.20)
project(SFML_Lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)   # вимірювання часу мають сенс лише з оптимізацією
endif()

# 1) Консольні експерименти (без SFML): створюють CSV для графіків
add_executable(experiments experiments.cpp)

# 1b) Побудова рисунків (SVG) з CSV: теж чистий C++, без залежностей
add_executable(make_plots make_plots.cpp)

# 2) Візуалізація на SFML (потребує vcpkg: sfml)
option(BUILD_SFML_DEMO "Build SFML visualisation" ON)
if(BUILD_SFML_DEMO)
  find_package(SFML COMPONENTS Graphics Window System CONFIG REQUIRED)
  add_executable(SFML_Lab main.cpp)
  target_link_libraries(SFML_Lab PRIVATE SFML::Graphics SFML::Window SFML::System)
endif()

# Автоматичне копіювання папки ресурсів у каталог з виконуваним файлом
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
