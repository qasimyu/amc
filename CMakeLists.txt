cmake_minimum_required(VERSION 3.0)
project(AMC)

SET(CMAKE_BUILD_TYPE Release)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11")

add_subdirectory(lib)

include_directories (
	${AMC_SOURCE_DIR}/lib/config
	${AMC_SOURCE_DIR}/lib/dkm/include
	${AMC_SOURCE_DIR}/lib/hetcaller
	${AMC_SOURCE_DIR}/lib/hclust
	${AMC_SOURCE_DIR}/lib/inputparser
	${AMC_SOURCE_DIR}/lib/mathfunc
	${AMC_SOURCE_DIR}/lib/matrix
	${AMC_SOURCE_DIR}/lib/mydefine
	${AMC_SOURCE_DIR}/lib/pca
	${AMC_SOURCE_DIR}/lib/split
	${AMC_SOURCE_DIR}/lib/threadpool
)

add_executable(amc src/amc.cpp)

target_link_libraries(amc mydefine)

SET_TARGET_PROPERTIES(amc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${AMC_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
	set_property(TARGET amc PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET amc PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
	target_link_libraries(amc "${CMAKE_THREAD_LIBS_INIT}")
endif()

INSTALL(TARGETS amc
	RUNTIME DESTINATION bin
)
