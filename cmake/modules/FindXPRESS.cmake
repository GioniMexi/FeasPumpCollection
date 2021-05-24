# FindXPRESS.cmake -- this module finds the FICO XPRESS optimizer
#
# Users might need to set/overwrite XPRESSDIR

# Find XPRESS headers and library
find_path(XPRESS_INCLUDE_DIR "xprs.h" PATHS ${XPRESSDIR}/include)
find_library(XPRESS_LIBRARY NAMES "xprs" PATHS ${XPRESSDIR}/lib)
find_library(XPRESS_XPRL_LIBRARY NAMES "xprl" PATHS ${XPRESSDIR}/lib)

mark_as_advanced(XPRESS_INCLUDE_DIR XPRESS_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XPRESS DEFAULT_MSG XPRESS_LIBRARY XPRESS_INCLUDE_DIR)

if (XPRESS_FOUND)
	# Create imported target Xpress::Xpress
	add_library(Xpress::Xpress SHARED IMPORTED)
	set_target_properties(Xpress::Xpress PROPERTIES
		IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
		IMPORTED_LOCATION "${XPRESS_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${XPRESS_INCLUDE_DIR}"
		INTERFACE_LINK_LIBRARIES "${XPRESS_XPRL_LIBRARY}"
	)
endif()
