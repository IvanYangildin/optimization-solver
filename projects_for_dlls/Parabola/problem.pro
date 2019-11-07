#-------------------------------------------------
#
# Project created by QtCreator 2019-09-05T23:02:19
#
#-------------------------------------------------

QMAKE_CXXFLAGS += -std=gnu++0x
QT       -= gui

TARGET = Parabola
TEMPLATE = lib

DEFINES += PARABOLA_LIBRARY SHARED_EXPORTS

SOURCES += \
    IProblemImplParabola.cpp

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
    TARGET.UID3 = 0xE2945B5F
    TARGET.CAPABILITY = 
    TARGET.EPOCALLOWDLLDATA = 1
    addFiles.sources = Parabola.dll
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










