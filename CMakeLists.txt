# CMake simple example

## [main]

# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.15...3.31)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(
	RG
	VERSION 1.0
	LANGUAGES CXX)


# Find packages go here.

# find ROOT

find_package(ROOT CONFIG REQUIRED)
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)


# SMDR package libraries
set(SMDR_INSTALL_DIR "/usr/local")


# static libraries 

add_library(smdr STATIC IMPORTED)
set_target_properties(smdr PROPERTIES
	IMPORTED_LOCATION "${SMDR_INSTALL_DIR}/lib/libsmdr.a"
	INTERFACE_INCLUDE_DIRECTORIES "${SMDR_INSTALL_DIR}/include/"
	)

add_library(tsil STATIC IMPORTED)
set_target_properties(tsil PROPERTIES
	IMPORTED_LOCATION "${SMDR_INSTALL_DIR}/lib/libtsil.a"
	INTERFACE_INCLUDE_DIRECTORIES "${SMDR_INSTALL_DIR}/include/"
	)
add_library(3vil STATIC IMPORTED)
set_target_properties(3vil PROPERTIES
	IMPORTED_LOCATION "${SMDR_INSTALL_DIR}/lib/lib3vil.a"
	INTERFACE_INCLUDE_DIRECTORIES "${SMDR_INSTALL_DIR}/include/"
	)


# define executables and the corresponding source files

add_executable(fit fit.cpp)
add_executable(init_random init_random.cpp)
add_executable(generate_random generate_random.cpp)
add_executable(generate_net generate_net.cpp)
add_executable(theor_err theor_err.cpp)


# define how to link executables 

target_link_libraries(fit  PUBLIC   ${ROOT_LIBRARIES}
				    ${ROOT_EXE_LINKER_FLAGS}
				    smdr
				    tsil
				    3vil)	
target_link_libraries(init_random  PUBLIC ${ROOT_LIBRARIES}
					  ${ROOT_EXE_LINKER_FLAGS}
				          smdr
	        			  tsil
					  3vil)	
target_link_libraries(generate_random  PUBLIC   ${ROOT_LIBRARIES}
						${ROOT_EXE_LINKER_FLAGS}
				          	smdr
	        			  	tsil
					  	3vil)	
target_link_libraries(generate_net  PUBLIC   ${ROOT_LIBRARIES}
						${ROOT_EXE_LINKER_FLAGS}
				          	smdr
	        			  	tsil
					  	3vil)	
target_link_libraries(theor_err  PUBLIC   smdr
	        			  tsil
					  3vil)	
