# FindCPLEX.cmake -- this module finds the IBM CPLEX optimizer
#
# Users need to set/overwrite CPLEX_ROOT_DIR and possibly overwrite CPLEX_LIBFORMAT

# Guess CPLEX_LIBFORMAT
if (NOT CPLEX_LIBFORMAT)
	if (UNIX)
		if (APPLE)
			set(CPLEX_LIBFORMAT "x86-64_osx/static_pic" CACHE PATH "CPLEX lib format")
		else()
			set(CPLEX_LIBFORMAT "x86-64_linux/static_pic" CACHE PATH "CPLEX lib format")
		endif()
	endif()
endif()

# Find CPLEX headers and library
find_path(CPLEX_INCLUDE_DIR "ilcplex/cplex.h" PATHS $ENV{CPLEX_ROOT_DIR}/include)
find_library(CPLEX_LIBRARY "libcplex.a" PATHS $ENV{CPLEX_ROOT_DIR}/lib/${CPLEX_LIBFORMAT})

mark_as_advanced(CPLEX_INCLUDE_DIR CPLEX_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

if (CPLEX_FOUND)
	# Create imported target Cplex::Cplex
	add_library(Cplex::Cplex STATIC IMPORTED)
	set_target_properties(Cplex::Cplex PROPERTIES
		IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
		IMPORTED_LOCATION "${CPLEX_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_INCLUDE_DIR}"
		INTERFACE_LINK_LIBRARIES "Threads::Threads;${CMAKE_DL_LIBS}"
	)
endif()
