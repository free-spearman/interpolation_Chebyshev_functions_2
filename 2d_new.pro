QT       += core gui

QT += opengl
LIBS += -lGLU
LIBS += -lpthread

#win32:LIBS += -lOpenGL32 -lglu32

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs de$

SOURCES += \
    chebyshev.cpp \
    help.cpp \
    main.cpp \
    mainwindow.cpp \
    myglwidget.cpp \
    auto_parallelization.cpp 

HEADERS += \
    chebyshev.hpp \
    help.hpp \
    mainwindow.h \
    myglwidget.h \
    auto_parallelization.hpp

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
