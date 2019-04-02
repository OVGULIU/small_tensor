cmake_minimum_required(VERSION 3.7.2)

project(vonMises VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 17)
include_directories(${PROJECT_SOURCE_DIR}/../../)


# build tsmap cpp library
file(GLOB VM_SRC
    ${PROJECT_SOURCE_DIR}/*.cpp
)

add_executable(vm
	${VM_SRC}
)

