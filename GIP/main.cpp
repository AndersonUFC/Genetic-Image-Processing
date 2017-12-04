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
    CImg<unsigned char> src("/home/nfs/andersonUFC/Pictures/lena.bmp");
    srand(time(NULL));

    size = src.size();
    data = new unsigned char[size];
    data = src.data();

    h = src.height();
    w = src.width();

    img = IG_new(data,size,w,h);

    int f = 1,b = 1;
    if(argc == 3){
        f = atoi(argv[1]);
        b = atoi(argv[2]);
    }
    f = 9;
    b = 9;

    //Image_Gene_Float* pau = IG_int_to_float(img);

    Image_Gene_Float* pau = IG_int_to_float(img);
    mult(50, pau);
    for(int i = 0 ; i < f ; i++)
        IG_haar2D_subdivide_float(pau);

    img = IG_float_to_int(pau);

    IG_predictive(img);
    sum(8500, img);
    IG_run_length(img);
    IG_huffman(img);

    IG_save_file(img, "/home/nfs/andersonUFC/Pictures/lena");
    std::cout << img->size;
    //IG_print(img);
    IG_read_file(img, "/home/nfs/andersonUFC/Pictures/lena");
    std::cout << img->size;

    IG_huffman_inv(img);
    IG_run_length_inv(img);
    sum(-8500, img);
    IG_predictive_inv(img);

    pau = IG_int_to_float(img);
    for(int i = 0 ; i < b ; i++)
            IG_haar2D_inv_subdivide_float(pau);
    mult(0.02, pau);
    img = IG_float_to_int(pau);
    //IG_print(img);


    //IG_haar2D_inv_subdivide(img);

    /*
    IG_haar2D_subdivide(img);
    IG_haar2D_subdivide(img);
    IG_haar2D_subdivide(img);
    IG_haar2D_subdivide(img);
    IG_haar2D_subdivide(img);
    IG_haar2D_subdivide(img);
    IG_haar2D_subdivide(img);

    IG_predictive(img);
    IG_run_length(img);
    IG_huffman(img);

    IG_print(img);

    /*
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);
    IG_haar2D(img);


    IG_huffman(img);
    IG_run_length(img);
    IG_huffman(img);

    IG_print(img);

    std::cout << "compression level: " << img->compression_level << std::endl;
    */

    //openGL
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(w,h);
    glutCreateWindow("Genetic Image Processing");

    glutDisplayFunc(display);
    init();
    glutMainLoop();


    return 0;

}
