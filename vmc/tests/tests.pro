TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lunittest++

SOURCES += main.cpp

include(../defaults.pri)

HEADERS += $$SRC_DIR/Slater.h \
           $$SRC_DIR/Orbitals.h \
           $$SRC_DIR/Jastrow.h

SOURCES += $$SRC_DIR/Slater.cpp \
           $$SRC_DIR/Orbitals.cpp \
           $$SRC_DIR/Jastrow.cpp
