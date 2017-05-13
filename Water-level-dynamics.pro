
CONFIG += console debug_and_release
CONFIG -= app_bundle
QT -= core gui
TEMPLATE = app

CONFIG(release, debug|release) {
  DEFINES += NDEBUG
}


unix:!macx{
  # Linux only
  message("Console application, built for Linux")
  message(Host name: $$QMAKE_HOST.name)
  QMAKE_CXX = g++-5
  QMAKE_LINK = g++-5
  QMAKE_CC = gcc-5
  QMAKE_CXXFLAGS += -Wall -Wextra -Weffc++ -std=c++14
  #QMAKE_CXXFLAGS += -Wall -Wextra -Weffc++ -Werror -std=c++14

  # gcov
  QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
  LIBS += -lgcov
}

HEADERS += \
    GetParams.h \
    Gillespie.h \
    randomc.h

SOURCES += \
    GetParams.cpp \
    Gillespie.cpp \
    mersenne.cpp \
    theta.cpp


