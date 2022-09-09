project(tbl)

add_subdirectory(lib)

include_directories(lib ${shared_path})
   
add_executable(tbl main.cpp)
target_link_libraries(tbl tblLib sharedCombiFFLib)

install(TARGETS tbl RUNTIME DESTINATION bin)
