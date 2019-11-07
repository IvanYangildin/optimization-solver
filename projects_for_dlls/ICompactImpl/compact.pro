#-------------------------------------------------
#
# Project created by QtCreator 2019-09-09T18:19:23
#
#-------------------------------------------------
QMAKE_CXXFLAGS += -std=gnu++0x
QT       -= gui

TARGET = ICompactImpl
TEMPLATE = lib

DEFINES += ICOMPACTIMPL_LIBRARY SHARED_EXPORTS

LIBS += debug/Ilog.dll
LIBS += debug/IVectorImpl.dll

SOURCES += \
    ICompactImpl.cpp \
    main.cpp

HEADERS += \
    SHARED_EXPORT.h \
    IVector.h \
    ISolver.h \
    ISet.h \
    IProblem.h \
    ILog.h \
    ICompact.h \
    IBrocker.h \
    error.h

symbian {
    MMP_RULES += EXPORTUNFROZEN
    TARGET.UID3 = 0xE2F00700
    TARGET.CAPABILITY = 
    TARGET.EPOCALLOWDLLDATA = 1
    addFiles.sources = ICompactImpl.dll
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







