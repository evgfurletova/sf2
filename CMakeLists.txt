cmake_minimum_required(VERSION 2.8)


find_package(PythonLibs 2.7)
include_directories(${PYTHON_INCLUDE_DIRS})

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

add_executable(sufpref2 ${SOURCES})
PYTHON_ADD_MODULE(sf ${SOURCES} ${PY_SOURCES})
#add_library(sf MODULE ${SOURCES} ${PY_SOURCES})
#target_link_libraries(sf ${PYTHON_LIBRARIES})
