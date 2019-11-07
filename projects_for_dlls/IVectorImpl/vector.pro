#-------------------------------------------------
#
# Project created by QtCreator 2019-09-09T17:45:24
#
#-------------------------------------------------
QMAKE_CXXFLAGS += -std=gnu++0x
QT       -= gui

TARGET = IVectorImpl
TEMPLATE = lib

DEFINES += IVECTORIMPL_LIBRARY SHARED_EXPORTS

SOURCES += \
    IVectorImpl.cpp

LIBS += debug/Ilog.dll

HEADERS += \
    SHARED_EXPORT.h \
    IVector.h \
    ILog.h \










