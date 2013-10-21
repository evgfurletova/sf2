cmake_minimum_required(VERSION 2.6)


find_package(PythonLibs 2.6)
find_package(PythonInterp 2.6)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_PATH})

add_definitions(-Dnullptr=0)

set(CMAKE_CXX_FLAGS "-std=c++0x")
set(HEADERS
d_hhm_prob.h
hmm_ovgraf.h
h_m_node.h
maindata.h
mmodel_prob.h
m_ovgraf.h
m_trtree.h
nd_hhm_prob.h
nodeac.h
nodebern.h
nodeov.h
statedata.h
type_malloc.h
sf_main.h
)

set(SOURCES
ac_trie.cpp
bprob_graf.cpp
d_hhm_prob.cpp
hmm_ovgraf.cpp
hmm_preprocessing.cpp
h_m_ovgraf.cpp
maindata.cpp
mmodel_prob.cpp
m_ovgraf.cpp
nd_hhm_prob.cpp
nodeinit.cpp
ov_graf.cpp
sufpref2.cpp
sf_main.cpp
)

set(PY_SOURCES
    sfmodule.c
)

add_executable(sufpref ${SOURCES})
PYTHON_ADD_MODULE(sf ${SOURCES} ${PY_SOURCES})
set_target_properties(sf PROPERTIES PREFIX "")

if((${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64") AND NOT (EXISTS /etc/ubuntu-release) AND NOT (EXISTS /etc/debian_version))
    set(LIB_BASENAME lib64)
else()
    set(LIB_BASENAME lib)
endif()

exec_program(
	${PYTHON_EXECUTABLE}
	ARGS -c \"import sys\; print '%d.%d' % (sys.version_info[0], sys.version_info[1])\"
	OUTPUT_VARIABLE PYTHON_VERSION
)

install(
	TARGETS sufpref sf
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION "${LIB_BASENAME}/python${PYTHON_VERSION}/site-packages"
)


