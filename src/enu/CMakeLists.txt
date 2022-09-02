project(enu)

add_subdirectory(lib)

include_directories(lib ${shared_path})

add_executable(enu main.cpp)
target_link_libraries(enu enuLib sharedCombiFFLib)

