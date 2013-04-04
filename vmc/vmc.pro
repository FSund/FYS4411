TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/lib.cpp \
    src/Jastrow.cpp \
    src/Orbitals.cpp \
    src/Slater.cpp \
    src/Wavefunction.cpp \
    src/Solver/SolverMCBF.cpp \
    src/Solver/Solver.cpp \
    src/Solver/SolverMCIS.cpp \
    src/Minimizer.cpp \
    src/VMCApp.cpp

HEADERS += \
    src/lib.h \
    src/Jastrow.h \
    src/Slater.h \
    src/Orbitals.h \
    src/Wavefunction.h \
    src/Solver/Solver.h \
    src/Solver/SolverMCIS.h \
    src/Solver/SolverMCBF.h \
    src/Minimizer.h \
    src/VMCApp.h

LIBS += -larmadillo

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}
