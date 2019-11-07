#-------------------------------------------------
#
# Project created by QtCreator 2019-09-05T23:28:14
#
#-------------------------------------------------
QMAKE_CXXFLAGS += -std=gnu++0x
QT       -= gui

TARGET = ISolverImpl1
TEMPLATE = lib

DEFINES += ISOLVERIMPL1_LIBRARY SHARED_EXPORTS

SOURCES += \
    ISolverImpl1.cpp

HEADERS += SHARED_EXPORT.h \
    IVector.h \
    ISolver.h \
    ISet.h \
    IProblem.h \
    ILog.h \
    ICompact.h \
    IBrocker.h \
    error.h

LIBS += debug/Ilog.dll
LIBS += debug/IVectorImpl.dll
LIBS += debug/ICompactImpl.dll


symbian {
    MMP_RULES += EXPORTUNFROZEN
    TARGET.UID3 = 0xE2EB6D4F
    TARGET.CAPABILITY = 
    TARGET.EPOCALLOWDLLDATA = 1
    addFiles.sources = ISolverImpl1.dll
    addFiles.path = !:/sys/bin
    DEPLOYMENT += addFiles
}

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}









