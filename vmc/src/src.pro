TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../defaults.pri)

TARGET = vmc

SOURCES += main.cpp \
    lib.cpp \
    Jastrow.cpp \
    Slater.cpp \
    Wavefunction.cpp \
    Solver/SolverMCBF.cpp \
    Solver/Solver.cpp \
    Solver/SolverMCIS.cpp \
    VMCApp.cpp \
    Datalogger.cpp \
    Localenergy/LocalEnergy.cpp \
    Localenergy/SingleAtomLocalEnergy.cpp \
    Localenergy/DiatomicLocalEnergy.cpp \
    Orbitals/Orbitals.cpp \
    Orbitals/Hydrogenic.cpp \
    Orbitals/Diatomic.cpp \
    Minimizer/Minimizer.cpp \
    Minimizer/BruteForce.cpp \
    Minimizer/SteepestDescent.cpp

HEADERS += \
    lib.h \
    Jastrow.h \
    Slater.h \
    Wavefunction.h \
    Solver/Solver.h \
    Solver/SolverMCIS.h \
    Solver/SolverMCBF.h \
    VMCApp.h \
    Datalogger.h \
    Localenergy/LocalEnergy.h \
    Localenergy/SingleAtomLocalEnergy.h \
    Localenergy/DiatomicLocalEnergy.h \
    Orbitals/Hydrogenic.h \
    Orbitals/Orbitals.h \
    Orbitals/Diatomic.h \
    Minimizer/Minimizer.h \
    Minimizer/BruteForce.h \
    Minimizer/SteepestDescent.h

