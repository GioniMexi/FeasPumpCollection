# FindSCIP.cmake -- this module finds SCIP

find_path(SCIP_INCLUDE_DIR "scip/scip.h" PATHS $ENV{SCIP_DIR}/include)

find_library(SCIP_LIBRARY scip PATHS $ENV{SCIP_DIR}/lib)

mark_as_advanced(SCIP_INCLUDE_DIR SCIP_LIBRARY)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCIP DEFAULT_MSG SCIP_LIBRARY SCIP_INCLUDE_DIR)


if (SCIP_FOUND)
	# Create imported target Scip::Scip
	add_library(Scip::Scip SHARED IMPORTED)
	set_target_properties(Scip::Scip PROPERTIES
		IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
		IMPORTED_LOCATION "${SCIP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SCIP_INCLUDE_DIR}"
	)
endif()