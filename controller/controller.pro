QMAKE_CXXFLAGS += -std=gnu++0x
QT       += gui sql printsupport
QT       += core declarative
QT       += qml

HEADERS += \
    ui_sqlconnectiondialog.h \
    ui_qplot.h \
    ui_dialog.h \
    sqlconnectiondialog.h \
    SHARED_EXPORT.h \
    qplot.h \
    paramsdialog.h \
    IVector.h \
    ISolver.h \
    ISet.h \
    IProblem.h \
    ILog.h \
    ICompact.h \
    IBrocker.h \
    error.h \
    controller.h \
    qcustomplot.h \

SOURCES += \
    sqlconnectiondialog.cpp \
    qplot.cpp \
    paramsdialog.cpp \
    main.cpp \
    controller.cpp \
    qcustomplot.cpp

FORMS += \
    sqlconnectiondialog.ui \
    qplot.ui \
    paramsdialog.ui \
    controller.ui

LIBS += debug/ICompactImpl.dll
LIBS += debug/ILog.dll
LIBS += debug/IVectorImpl.dll
LIBS += debug/ISetImpl.dll






