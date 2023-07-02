find_path(ORTOOLS_INCLUDE_DIRS
    NAMES ortools/
    HINTS ${ORTOOLS_DIR} $ENV{ORTOOLS_DIR}
    PATH_SUFFIXES include)

set(ORTOOLS_INCLUDE_DIRS
    ${ORTOOLS_INCLUDE_DIRS}
    )

# todo: enable recursive search
find_library(ORTOOLS_LIBRARY
    NAMES ortools
    HINTS ${ORTOOLS_DIR} $ENV{ORTOOLS_DIR}
    PATH_SUFFIXES lib)

find_library(GLOG_LIBRARY
    NAMES glog
    HINTS ${ORTOOLS_DIR}/dependencies/install/lib $ENV{ORTOOLS_DIR}/dependencies/install/lib
    HINTS ${ORTOOLS_DIR} $ENV{ORTOOLS_DIR}
    PATH_SUFFIXES lib)

if(GLOG_LIBRARY)
    set(ORTOOLS_LIBRARIES ${ORTOOLS_LIBRARY} ${GLOG_LIBRARY})
else()
    set(ORTOOLS_LIBRARIES ${ORTOOLS_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ORTOOLS DEFAULT_MSG ORTOOLS_INCLUDE_DIRS ORTOOLS_LIBRARIES)

if (ORTOOLS_FOUND)
	# Create imported target Pdlp::Pdlp
	add_library(Pdlp::Pdlp SHARED IMPORTED)
	set_target_properties(Pdlp::Pdlp PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "${ORTOOLS_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${ORTOOLS_INCLUDE_DIRS}"
)
endif()