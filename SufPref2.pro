TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11

SOURCES += \
    sufpref2.cpp \
    ov_graf.cpp \
    nodeinit.cpp \
    nd_hhm_prob.cpp \
    mmodel_prob.cpp \
    maindata.cpp \
    m_ovgraf.cpp \
    hmm_preprocessing.cpp \
    hmm_ovgraf.cpp \
    h_m_ovgraf.cpp \
    d_hhm_prob.cpp \
    bprob_graf.cpp \
    ac_trie.cpp

HEADERS += \
    type_malloc.h \
    statedata.h \
    nodeov.h \
    nodebern.h \
    nodeac.h \
    nd_hhm_prob.h \
    mmodel_prob.h \
    maindata.h \
    m_trtree.h \
    m_ovgraf.h \
    hmm_ovgraf.h \
    h_m_node.h \
    d_hhm_prob.h

