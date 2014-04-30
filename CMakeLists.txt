cmake_minimum_required(VERSION 2.6)
set(CMAKE_VERBOSE_MAKEFILE on)


find_package(PythonLibs 2.6)
find_package(PythonInterp 2.6)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_PATH})

add_definitions(-Dnullptr=0)

set(CMAKE_CXX_FLAGS "-std=c++0x -O0")
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


add_executable(sufpref ${SOURCES})
target_link_libraries(sufpref -lstdc++ -lm)

if(NOT WIN32)

set(PY_SOURCES
    sfmodule.c
)


PYTHON_ADD_MODULE(sf ${SOURCES} ${PY_SOURCES})
set_target_properties(sf PROPERTIES PREFIX "")
target_link_libraries(sf -lstdc++ -lm)

exec_program(
	${PYTHON_EXECUTABLE}
	ARGS -c \"from distutils.sysconfig import get_python_lib\; print\(get_python_lib\(1\)[4:]\)\"
	OUTPUT_VARIABLE PYTHON_SITEARCH
)

install(
	TARGETS sufpref sf
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}${PYTHON_SITEARCH}"
)

if(EXISTS /etc/SuSE-release)
    set(DOCDIR /share/doc/packages/sf2)
else()
    set(DOCDIR /share/doc/sf2-2.0.0)
endif()

install(
	FILES py_module.html
	DESTINATION "${CMAKE_INSTALL_PREFIX}${DOCDIR}"
)

endif()

