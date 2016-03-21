# search path WignerSymbols_ROOT
#

find_path(WignerSymbols_INCLUDE_DIRS
NAME wignerSymbols.h
HINTS /usr/local/include /usr/include ${WignerSymbols_ROOT}/include
DOC "Header files of Wigner Symbols library")

find_library(WignerSymbols_LIBRARIES
NAME wignerSymbols
HINTS /usr/local/lib/ usr/lib ${WignerSymbols_ROOT}/lib
DOC "Libraries of Wigner Symbols")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(WignerSymbols WignerSymbols_INCLUDE_DIRS WignerSymbols_LIBRARIES)

mark_as_advanced(WignerSymbols_INCLUDE_DIRS WignerSymbols_LIBRARIES)
