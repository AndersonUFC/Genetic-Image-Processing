#-------------------------------------------------
#
# Project created by QtCreator 2017-11-20T14:18:12
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = GIP
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    gip.cpp \

HEADERS  += mainwindow.h \
    gip.h \
    plugins/add_fileformat.h \
    plugins/bayer.h \
    plugins/chlpca.h \
    plugins/cvMat.h \
    plugins/draw_gradient.h \
    plugins/inpaint.h \
    plugins/ipl.h \
    plugins/ipl_alt.h \
    plugins/jpeg_buffer.h \
    plugins/loop_macros.h \
    plugins/matlab.h \
    plugins/nlmeans.h \
    plugins/skeleton.h \
    plugins/tiff_stream.h \
    plugins/tinymatwriter.h \
    plugins/vrml.h \
    plugins/vtk.h \
    CImg.h \
    cimg/examples/img/CImg_demo.h \
    cimg/examples/img/odykill.h \
    cimg/examples/img/tetris.h \
    cimg/plugins/add_fileformat.h \
    cimg/plugins/bayer.h \
    cimg/plugins/chlpca.h \
    cimg/plugins/cvMat.h \
    cimg/plugins/draw_gradient.h \
    cimg/plugins/inpaint.h \
    cimg/plugins/ipl.h \
    cimg/plugins/ipl_alt.h \
    cimg/plugins/jpeg_buffer.h \
    cimg/plugins/loop_macros.h \
    cimg/plugins/matlab.h \
    cimg/plugins/nlmeans.h \
    cimg/plugins/skeleton.h \
    cimg/plugins/tiff_stream.h \
    cimg/plugins/tinymatwriter.h \
    cimg/plugins/vrml.h \
    cimg/plugins/vtk.h \
    cimg/CImg.h

FORMS    += mainwindow.ui

DISTFILES += \
    lena.bmp \
    cimg/resources/CImg_reference.pdf \
    cimg/resources/compile_win_icl.bat \
    cimg/resources/compile_win_visualcpp.bat \
    cimg/Licence_CeCILL-C_V1-en.txt \
    cimg/Licence_CeCILL_V2-en.txt \
    cimg/README.txt \

LIBS += -lX11 -ldl -lXext -lz -lGLU -lGL -lglut
