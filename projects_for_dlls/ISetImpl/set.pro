#-------------------------------------------------
#
# Project created by QtCreator 2019-09-09T18:04:18
#
#-------------------------------------------------
QMAKE_CXXFLAGS += -std=gnu++0x
QT       -= gui

TARGET = ISetImpl
TEMPLATE = lib

DEFINES += ISETIMPL_LIBRARY SHARED_EXPORTS


LIBS += debug/Ilog.dll
LIBS += debug/IVectorImpl.dll

SOURCES += \
    ISetImpl.cpp

HEADERS += \
    IVector.h \
    ISet.h \
    ILog.h \
    error.h \
    SHARED_EXPORT.h \
    ISolver.h \
    IProblem.h \
    ICompact.h \
    IBrocker.h \

symbian {
    MMP_RULES += EXPORTUNFROZEN
    TARGET.UID3 = 0xE24D5C49
    TARGET.CAPABILITY = 
    TARGET.EPOCALLOWDLLDATA = 1
    addFiles.sources = ISetImpl.dll
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








