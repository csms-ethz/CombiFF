project(cnv)

add_subdirectory(lib)

include_directories(lib ${shared_path})

add_executable(cnv main.cpp)
target_link_libraries(cnv cnvLib sharedLib)

