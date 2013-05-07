TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lunittest++

SOURCES += main.cpp

include(../defaults.pri)

HEADERS += $$SRC_DIR/Slater.h \
           $$SRC_DIR/Orbitals.h \
           $$SRC_DIR/Jastrow.h \
           $$SRC_DIR/Wavefunction.h \
           $$SRC_DIR/Solver/Solver.h \
           $$SRC_DIR/Solver/SolverMCBF.h \
           $$SRC_DIR/Solver/SolverMCIS.h \
           $$SRC_DIR/lib.h \
           $$SRC_DIR/Datalogger.h \


SOURCES += $$SRC_DIR/Slater.cpp \
           $$SRC_DIR/Orbitals.cpp \
           $$SRC_DIR/Jastrow.cpp \
           $$SRC_DIR/Wavefunction.cpp \
           $$SRC_DIR/Solver/Solver.cpp \
           $$SRC_DIR/Solver/SolverMCBF.cpp \
           $$SRC_DIR/Solver/SolverMCIS.cpp \
           $$SRC_DIR/lib.cpp \
           $$SRC_DIR/Datalogger.cpp \
