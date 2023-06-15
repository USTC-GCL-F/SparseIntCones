cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)

add_definitions("/DNOMINMAX")

add_subdirectory(SparseCones)

add_executable(exeSparseCones main.cpp)
target_link_libraries(exeSparseCones PUBLIC
  SparseCones
)
