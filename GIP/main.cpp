#include "cimg/CImg.h"
#include "gip.h"
#include <iostream>

#include <GL/glew.h>
#include <GL/glut.h>

using namespace cimg_library;

// DATA
int h,w,size;
unsigned char* data;
Image_Gene* img;

#define PI 3.1415926535

void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glOrtho(-1,1,1,-1,-1,1);

    configure_data(img->data, img->size, &img->data);
    // convert to float
    float* ddata = new float[img->size];
    for(int i = 0 ; i < img->size ; i++)
        ddata[i] = ((float)img->data[i])/255;

    // flip y dimension
    for(int y = 0 ; y < img->height/2 ; y++){
        for(int x = 0 ; x < img->width ; x++){
            int position1 = (y*img->width + x)*3;
            int position2 = ((img->height-y-1)*img->width + x)*3;

            float tr = ddata[position1];
            float tg = ddata[position1+1];
            float tb = ddata[position1+2];

            ddata[position1] = ddata[position2];
            ddata[position1+1] = ddata[position2+1];
            ddata[position1+2] = ddata[position2+2];

            ddata[position2] = tr;
            ddata[position2+1] = tg;
            ddata[position2+2] = tb;
        }// for
    }// for

    glDrawPixels(w,h,GL_RGB,GL_FLOAT,ddata);
    glutSwapBuffers();
}

void init(){
    glClearColor(0,0,0,0);
}

int main(int argc, char *argv[]){

    CImg<unsigned char> src("lena.bmp");

    size = src.size();
    data = new unsigned char[size];
    data = src.data();

    h = src.height();
    w = src.width();

    img = IG_new(data,size,w,h);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);

    //IG_print(img);

    IG_huffman(img);
    IG_run_length(img);

    IG_print(img);

    /*
    //openGL
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(w,h);
    glutCreateWindow("Genetic Image Processing");

    glutDisplayFunc(display);
    init();
    glutMainLoop();

    */

    int k;
    std::cin >> k;

    /*
    unsigned char* d = new unsigned char[12];

    d[0] = 127;
    d[1] = 128;
    d[2] = 129;
    d[3] = 128;

    d[4] = 115;
    d[5] = 118;
    d[6] = 117;
    d[7] = 120;

    d[8] = 0;
    d[9] = 1;
    d[10] = 2;
    d[11] = 3;

    Image_Gene* ig = IG_new(d, 12, 12, 1);

    IG_print(ig);
    IG_predictive(ig);
    IG_print(ig);
    IG_run_length(ig);
    IG_print(ig);

    */

    return 0;

}
