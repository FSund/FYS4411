TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../defaults.pri)

TARGET = MyProject

SOURCES += main.cpp \
    lib.cpp \
    Jastrow.cpp \
    Orbitals.cpp \
    Slater.cpp \
    Wavefunction.cpp \
    Solver/SolverMCBF.cpp \
    Solver/Solver.cpp \
    Solver/SolverMCIS.cpp \
    Minimizer.cpp \
    VMCApp.cpp

HEADERS += \
    lib.h \
    Jastrow.h \
    Slater.h \
    Orbitals.h \
    Wavefunction.h \
    Solver/Solver.h \
    Solver/SolverMCIS.h \
    Solver/SolverMCBF.h \
    Minimizer.h \
    VMCApp.h

