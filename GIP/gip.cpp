#include "gip.h"

#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>

// ================================================================================================
// ================================================================================================
// VARIABLES ======================================================================================
// ================================================================================================
// ================================================================================================

std::map<int, IG_Binary_Tree*> bt_coding;

// ================================================================================================
// ================================================================================================
// CONSTRUCTOR ====================================================================================
// ================================================================================================
// ================================================================================================

Image_Gene* IG_new(unsigned char* data, int size, int width, int height){
    Image_Gene* ig = new Image_Gene;

    ig->compression_level = (int)(log2(width));
    ig->max_compression_level = ig->compression_level;
    ig->subdivisions = 0;
    ig->size = size;

    ig->data = new int[size];
    for(int i = 0 ; i < size ; i++){
        ig->data[i] = (int)data[i];
    }// for

    ig->width = width;
    ig->height = height;

    return ig;
}// IG_new

// ================================================================================================
// ================================================================================================
// CONVERSION =====================================================================================
// ================================================================================================
// ================================================================================================

void data_to_RGB(unsigned char* data, int size, unsigned char** R, unsigned char** G, unsigned char** B){
    unsigned char *R_ = new unsigned char[size/3],
         *G_ = new unsigned char[size/3],
         *B_ = new unsigned char[size/3];

    for(int i = 0 ; i < size/3 ; i++){
        R_[i] = data[i*3];
        G_[i] = data[i*3 + 1];
        B_[i] = data[i*3 + 2];
    }//for

    *R = R_;
    *G = G_;
    *B = B_;
}// data_to_RGB

void configure_data(int* data, int size, int** new_data){
    int *d = new int[size];

    for(int i = 0 ; i < size/3 ; i++){
        d[i*3] = data[i];
        d[i*3+1] = data[(size/3)+ i];
        d[i*3+2] = data[((2*size)/3)+ i];
    }// for

    *new_data = d;
}// configure data

Image_Gene_Float* IG_int_to_float(Image_Gene* ig){
    Image_Gene_Float* kek = new Image_Gene_Float;

    kek->compression_level = ig->compression_level;
    kek->height = ig->height;
    kek->max_compression_level = ig->max_compression_level;
    kek->size = ig->size;
    kek->subdivisions = ig->subdivisions;
    kek->width = ig->width;

    kek->data = new float[kek->size];

    for(int i = 0 ; i < kek->size; i++)
        kek->data[i] = (float)ig->data[i];

    return kek;
}

Image_Gene* IG_float_to_int(Image_Gene_Float* ig){
    Image_Gene* kek = new Image_Gene;

    kek->compression_level = ig->compression_level;
    kek->height = ig->height;
    kek->max_compression_level = ig->max_compression_level;
    kek->size = ig->size;
    kek->subdivisions = ig->subdivisions;
    kek->width = ig->width;

    kek->data = new int[kek->size];

    for(int i = 0 ; i < kek->size; i++)
        kek->data[i] = (int)ig->data[i];

    return kek;
}

// ================================================================================================
// ================================================================================================
// TRANSFORMATION =================================================================================
// ================================================================================================
// ================================================================================================

int* IG_haar1D(int* image, int size, int compression_level){
    int* compression = new int[size];
    int comp_size = pow(2, compression_level-1);

    for(int i = 0 ; i < comp_size ; i++){
        compression[i] = (image[2*i] + image[2*i+1])/2;
        compression[comp_size + i] = image[2*i] - compression[i];
    }// for

    for(int i = comp_size*2 ; i < size ; i++)
        compression[i] = image[i];

    return compression;
}// haar1D

int* IG_haar1D_inv(int* image, int size, int compression_level){
    int* compression = new int[size];
    int comp_size = pow(2, compression_level);

    for(int i = 0 ; i < comp_size ; i++){
        compression[2*i] = image[i] + image[comp_size + i];
        compression[2*i + 1] = image[i] - image[comp_size + i];
    }// for

    for(int i = comp_size*2 ; i < size ; i++)
        compression[i] = image[i];

    return compression;
}// haar1D_inv

void IG_haar2D(Image_Gene* ig){

    int bound = pow(2, ig->compression_level);
    int* lineR = new int[bound];
    int* lineG = new int[bound];
    int* lineB = new int[bound];

    int height = ig->height;
    int width = ig->width;
    int* data = ig->data;

    int g_index = ig->size/3;
    int b_index = ((2*ig->size)/3);


    // LINES
    for(int y = 0 ; y < bound ; y++){
        // R sequence, G sequence, B sequence
        // get line
        for(int x = 0 ; x < bound ; x++){
            int position = y*width + x;

            lineR[x] = data[position];
            lineG[x] = data[position+g_index];
            lineB[x] = data[position+b_index];

        }// for-x


        int* cR = IG_haar1D(lineR,bound,ig->compression_level);
        int* cG = IG_haar1D(lineG,bound,ig->compression_level);
        int* cB = IG_haar1D(lineB,bound,ig->compression_level);

        for(int x = 0 ; x < bound ; x++){
            int position = y*width + x;
            data[position] = cR[x];
            data[position+g_index] = cG[x];
            data[position+b_index] = cB[x];
        }// for-x

        free(cR);
        free(cG);
        free(cB);
        cR = NULL;
        cG = NULL;
        cB = NULL;
    }// for-y

    int height_bound = bound;
    // COLUMNS
    for(int x = 0 ; x < bound ; x++){
        // get column
        for(int y = 0 ; y < height_bound ; y++){
            int position = y*width + x;

            lineR[y] = data[position];
            lineG[y] = data[position+g_index];
            lineB[y] = data[position+b_index];

        }// for-x

        int* cR = IG_haar1D(lineR,height_bound,ig->compression_level);
        int* cG = IG_haar1D(lineG,height_bound,ig->compression_level);
        int* cB = IG_haar1D(lineB,height_bound,ig->compression_level);

        for(int y = 0 ; y < height_bound; y++){
            int position = y*width + x;
            data[position] = cR[y];
            data[position+g_index] = cG[y];
            data[position+b_index] = cB[y];
        }// for-x

        free(cR);
        free(cG);
        free(cB);
        cR = NULL;
        cG = NULL;
        cB = NULL;
    }// for-y



    ig->data = data;
    ig->compression_level--;
    ig->subdivisions++;

}// haar2D

float* IG_haar1D_float(float* image, int size, int compression_level){
    float* compression = new float[size];
    int comp_size = pow(2, compression_level-1);

    for(int i = 0 ; i < comp_size ; i++){
        compression[i] = (image[2*i] + image[2*i+1])/2;
        compression[comp_size + i] = image[2*i] - compression[i];
    }// for

    for(int i = comp_size*2 ; i < size ; i++)
        compression[i] = image[i];

    return compression;
}// haar1D

float* IG_haar1D_float_inv(float* image, int size, int compression_level){
    float* compression = new float[size];
    int comp_size = pow(2, compression_level);

    for(int i = 0 ; i < comp_size ; i++){
        compression[2*i] = image[i] + image[comp_size + i];
        compression[2*i + 1] = image[i] - image[comp_size + i];
    }// for

    for(int i = comp_size*2 ; i < size ; i++)
        compression[i] = image[i];

    return compression;
}// haar1D_inv

void IG_haar2D_subdivide_float(Image_Gene_Float* ig){
    int bound = pow(2, ig->compression_level);
    float* lineR = new float[bound];
    float* lineG = new float[bound];
    float* lineB = new float[bound];

    int width = ig->width;
    float* data = ig->data;

    int g_index = ig->size/3;
    int b_index = ((2*ig->size)/3);

    int region_step_index = width/pow(2, ig->subdivisions);

    for(int sub_vert = 0 ; sub_vert < pow(2, ig->subdivisions) ; sub_vert++){
        for(int sub_hor = 0 ; sub_hor < pow(2, ig->subdivisions) ; sub_hor++){
            // LINES
            for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){

                // R sequence, G sequence, B sequence
                // get line
                int i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                    int position = y*width + x;

                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x

                float* cR = IG_haar1D_float(lineR,region_step_index,ig->compression_level);
                float* cG = IG_haar1D_float(lineG,region_step_index,ig->compression_level);
                float* cB = IG_haar1D_float(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                    int position = y*width + x;

                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x
            }// for-y

            // COLUMNS
            for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                int i = 0;
                // get column
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){
                    int position = y*width + x;

                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x

                float* cR = IG_haar1D_float(lineR,region_step_index,ig->compression_level);
                float* cG = IG_haar1D_float(lineG,region_step_index,ig->compression_level);
                float* cB = IG_haar1D_float(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1); y++){
                    int position = y*width + x;
                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x
            }// for-y
        } // horizontal subdivisions
    }// vertical subdivisions

    ig->data = data;
    ig->compression_level--;
    ig->subdivisions++;
}// haar2D

void IG_haar2D_inv_subdivide_float(Image_Gene_Float* ig){
    int bound = pow(2, ig->compression_level+1);
    float* lineR = new float[bound];
    float* lineG = new float[bound];
    float* lineB = new float[bound];
    int width = ig->width;
    float* data = ig->data;

    int g_index = ig->size/3;
    int b_index = ((2*ig->size)/3);
    int region_step_index = width/pow(2, ig->subdivisions-1);

    for(int sub_vert = 0 ; sub_vert < pow(2, ig->subdivisions-1) ; sub_vert++){
        for(int sub_hor = 0 ; sub_hor < pow(2, ig->subdivisions-1) ; sub_hor++){
            // LINES
            for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){

                // R sequence, G sequence, B sequence
                // get line
                int i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){

                    int position = y*width + x;
                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x


                float* cR = IG_haar1D_float_inv(lineR,region_step_index,ig->compression_level);
                float* cG = IG_haar1D_float_inv(lineG,region_step_index,ig->compression_level);
                float* cB = IG_haar1D_float_inv(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                    int position = y*width + x;

                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x

            }// for-y


            // COLUMNS
            for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                // get column
                int i = 0;
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){
                    int position = y*width + x;

                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x

                float* cR = IG_haar1D_float_inv(lineR,region_step_index,ig->compression_level);
                float* cG = IG_haar1D_float_inv(lineG,region_step_index,ig->compression_level);
                float* cB = IG_haar1D_float_inv(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1); y++){
                    int position = y*width + x;
                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x
            }// for-y

        } // horizontal subdivisions
    }// vertical subdivisions

    ig->data = data;
    ig->compression_level++;
    ig->subdivisions--;
}// inv haar2D

// -------------------------------------------------------------------------------------------------


void IG_haar2D_subdivide(Image_Gene* ig){
    int bound = pow(2, ig->compression_level);
    int* lineR = new int[bound];
    int* lineG = new int[bound];
    int* lineB = new int[bound];

    int width = ig->width;
    int* data = ig->data;

    int g_index = ig->size/3;
    int b_index = ((2*ig->size)/3);

    int region_step_index = width/pow(2, ig->subdivisions);

    for(int sub_vert = 0 ; sub_vert < pow(2, ig->subdivisions) ; sub_vert++){
        for(int sub_hor = 0 ; sub_hor < pow(2, ig->subdivisions) ; sub_hor++){
            // LINES
            for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){

                // R sequence, G sequence, B sequence
                // get line
                int i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                    int position = y*width + x;

                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x

                int* cR = IG_haar1D(lineR,region_step_index,ig->compression_level);
                int* cG = IG_haar1D(lineG,region_step_index,ig->compression_level);
                int* cB = IG_haar1D(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                    int position = y*width + x;

                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x
            }// for-y

            // COLUMNS
            for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                int i = 0;
                // get column
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){
                    int position = y*width + x;

                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x

                int* cR = IG_haar1D(lineR,region_step_index,ig->compression_level);
                int* cG = IG_haar1D(lineG,region_step_index,ig->compression_level);
                int* cB = IG_haar1D(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1); y++){
                    int position = y*width + x;
                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x
            }// for-y
        } // horizontal subdivisions
    }// vertical subdivisions

    ig->data = data;
    ig->compression_level--;
    ig->subdivisions++;
}// haar2D

void IG_haar2D_inv_subdivide(Image_Gene* ig){
    int bound = pow(2, ig->compression_level+1);
    int* lineR = new int[bound];
    int* lineG = new int[bound];
    int* lineB = new int[bound];
    int width = ig->width;
    int* data = ig->data;

    int g_index = ig->size/3;
    int b_index = ((2*ig->size)/3);
    int region_step_index = width/pow(2, ig->subdivisions-1);

    for(int sub_vert = 0 ; sub_vert < pow(2, ig->subdivisions-1) ; sub_vert++){
        for(int sub_hor = 0 ; sub_hor < pow(2, ig->subdivisions-1) ; sub_hor++){
            // LINES
            for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){

                // R sequence, G sequence, B sequence
                // get line
                int i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){

                    int position = y*width + x;
                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x


                int* cR = IG_haar1D_inv(lineR,region_step_index,ig->compression_level);
                int* cG = IG_haar1D_inv(lineG,region_step_index,ig->compression_level);
                int* cB = IG_haar1D_inv(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                    int position = y*width + x;

                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x

            }// for-y


            // COLUMNS
            for(int x = region_step_index*sub_hor ; x < region_step_index*(sub_hor+1) ; x++){
                // get column
                int i = 0;
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1) ; y++){
                    int position = y*width + x;

                    lineR[i] = data[position];
                    lineG[i] = data[position+g_index];
                    lineB[i] = data[position+b_index];
                    i++;
                }// for-x

                int* cR = IG_haar1D_inv(lineR,region_step_index,ig->compression_level);
                int* cG = IG_haar1D_inv(lineG,region_step_index,ig->compression_level);
                int* cB = IG_haar1D_inv(lineB,region_step_index,ig->compression_level);

                i = 0;
                for(int y = region_step_index*sub_vert ; y < region_step_index*(sub_vert+1); y++){
                    int position = y*width + x;
                    data[position] = cR[i];
                    data[position+g_index] = cG[i];
                    data[position+b_index] = cB[i];
                    i++;
                }// for-x
            }// for-y

        } // horizontal subdivisions
    }// vertical subdivisions

    ig->data = data;
    ig->compression_level++;
    ig->subdivisions--;
}// inv haar2D

// ------------------------------------------------------------------------------------------------


// ================================================================================================
// ================================================================================================
// COMPRESSION ====================================================================================
// ================================================================================================
// ================================================================================================

void IG_predictive(Image_Gene* ig){
    int* data = ig->data;
    int size = ig->size;
    int* new_data = new int[size];

    int g_index = size/3;
    int b_index = (2*size)/3;
    new_data[0] = data[0];
    new_data[g_index] = data[g_index];
    new_data[b_index] = data[b_index];

    for(int i = 1 ; i < g_index ; i++){
        new_data[i] = data[i] - data[i-1];
        new_data[i+g_index] = data[i+g_index] - data[i-1+g_index];
        new_data[i+b_index] = data[i+b_index] - data[i-1+b_index];
    }// for

    ig->data = new_data;
}// IG_predictive

void IG_predictive_inv(Image_Gene* ig){
    int* data = ig->data;
    int size = ig->size;
    int* new_data = new int[size];

    int g_index = size/3;
    int b_index = (2*size)/3;
    new_data[0] = data[0];
    new_data[g_index] = data[g_index];
    new_data[b_index] = data[b_index];

    for(int i = 1 ; i < g_index ; i++){
        new_data[i] = data[i] + data[i-1];
        new_data[i+g_index] = data[i+g_index] + data[i-1+g_index];
        new_data[i+b_index] = data[i+b_index] + data[i-1+b_index];
    }// for

    ig->data = new_data;
}// IG_predictive_Inverse

void IG_run_length(Image_Gene* ig){
    int* data = ig->data;
    int size = ig->size;

    std::vector<int> rlr, rlg, rlb;
    int current_r, current_g, current_b,
        acc_r = 1, acc_g = 1, acc_b = 1;

    int g_index = size/3;
    int b_index = (2*size)/3;

    current_r = data[0];
    current_g = data[g_index];
    current_b = data[b_index];

    rlr.push_back(current_r);
    rlg.push_back(current_g);
    rlb.push_back(current_b);
    for(int i = 1 ; i < g_index ; i++){

        if(data[i] == current_r){
            acc_r++;
        } else {
            rlr.push_back(acc_r);
            acc_r = 1;
            current_r = data[i];
            rlr.push_back(current_r);
        }// if-else

        if(data[i+g_index] == current_g){
            acc_g++;
        } else {
            rlg.push_back(acc_g);
            acc_g = 1;
            current_g = data[i+g_index];
            rlg.push_back(current_g);
        }// if-else

        if(data[i+b_index] == current_b){
            acc_b++;
        } else {

            rlb.push_back(acc_b);
            acc_b = 1;
            current_b = data[i+b_index];
            rlb.push_back(current_b);
        }// if-else
    } //for

    rlr.push_back(acc_r);
    rlg.push_back(acc_g);
    rlb.push_back(acc_b);
    int rs = rlr.size(), gs = rlg.size(), bs = rlb.size();
    int* new_data = new int[rs + gs + bs];

    int* rdata = rlr.data(), *gdata = rlg.data(), *bdata = rlb.data();


    for(int i = 0 ; i < rs ; i++){
        new_data[i] = rdata[i];
    }

    for(int i = 0 ; i < gs ; i++){
        new_data[i+rs] = gdata[i];
    }

    for(int i = 0 ; i < bs ; i++){
        new_data[i+rs+gs] = bdata[i];
    }

    ig->data = new_data;
    ig->size = rs+gs+bs;

}//IG_run_length

void IG_run_length_inv(Image_Gene* ig){
    int cor = 0;
    std::vector<int> new_data_v;
    for(int i=0, size=ig->size; i<size; i+=2){
        cor = ig->data[i];
        for(int j=0, size_j=ig->data[i+1]; j<size_j; j++){
            new_data_v.push_back(cor);
        }
    }

    int je = new_data_v.size();
    free(ig->data);
    ig->data = new int[je];

    for(int i = 0 ; i < je ; i++){
        ig->data[i] = new_data_v.at(i);
    }
    ig->size = je;

}//IG_run_length_Inverse

void IG_run_length_byte(Image_Gene* ig){
    int* data = ig->data;
    int size = ig->size;

    std::vector<int> rlr, rlg, rlb;
    int current_r, current_g, current_b,
        acc_r = 1, acc_g = 1, acc_b = 1;

    int g_index = size/3;
    int b_index = (2*size)/3;

    current_r = data[0];
    current_g = data[g_index];
    current_b = data[b_index];

    rlr.push_back(current_r);
    rlg.push_back(current_g);
    rlb.push_back(current_b);
    for(int i = 1 ; i < g_index ; i++){

        if(data[i] == current_r){
            acc_r++;
            if(acc_r>255){
                rlr.push_back(255);
                acc_r = 1;
                rlr.push_back(current_r);
            }
        } else {
            rlr.push_back(acc_r);
            acc_r = 1;
            current_r = data[i];
            rlr.push_back(current_r);
        }// if-else

        if(data[i+g_index] == current_g){
            acc_g++;
            if(acc_g>255){
                rlg.push_back(255);
                acc_g = 1;
                rlg.push_back(current_g);
            }
        } else {
            rlg.push_back(acc_g);
            acc_g = 1;
            current_g = data[i+g_index];
            rlg.push_back(current_g);
        }// if-else

        if(data[i+b_index] == current_b){
            acc_b++;
            if(acc_b>255){
                rlb.push_back(255);
                acc_b = 1;
                rlb.push_back(current_b);
            }
        } else {

            rlb.push_back(acc_b);
            acc_b = 1;
            current_b = data[i+b_index];
            rlb.push_back(current_b);
        }// if-else
    } //for

    rlr.push_back(acc_r);
    rlg.push_back(acc_g);
    rlb.push_back(acc_b);
    int rs = rlr.size(), gs = rlg.size(), bs = rlb.size();
    int* new_data = new int[rs + gs + bs];

    int* rdata = rlr.data(), *gdata = rlg.data(), *bdata = rlb.data();


    for(int i = 0 ; i < rs ; i++){
        new_data[i] = rdata[i];
    }

    for(int i = 0 ; i < gs ; i++){
        new_data[i+rs] = gdata[i];
    }

    for(int i = 0 ; i < bs ; i++){
        new_data[i+rs+gs] = bdata[i];
    }

    ig->data = new_data;
    ig->size = rs+gs+bs;

}//IG_run_length

void sum(int value, Image_Gene* ig){

    for(int i = 0 ; i < ig->size; i++)
        ig->data[i] = ig->data[i] + value;
}


void int_to_binary_array(int value, std::vector<int>* new_data_v){
    //std::cout << value << " ----------------------------------------------------------------\n";
    //for(int i=0;i<32;++i)
    for(int i=17;i>=0;--i)
    {
        //std::cout << ((value >> i) & 1);
        new_data_v->push_back((value >> i) & 1);
    }
    //std::cout << "\n";
}

void IG_huffman(Image_Gene* ig){
    std::map<int,int> elements;
    int* data = ig->data;
    int size = ig->size;

    for(int i = 0 ; i < size ; i++){
        if(elements.find(data[i]) == elements.end()){
            elements[data[i]] = 1;
        } else {
            elements[data[i]]++;
        }// if-else
    }// for

    std::map<int, IG_Binary_Tree*> huff_tree;

    for(std::pair<int,int> p : elements){
        float prob = ((float)p.second)/((float)size);
        huff_tree[p.first] = IG_BT_new(p.first, prob);
    }// for

    while(huff_tree.size() > 1){
        int first_index = 0, second_index = 0;
        float first_value = 1, second_value = 1;

        // get minimum values
        for(std::pair<int, IG_Binary_Tree*> p: huff_tree){
            IG_Binary_Tree* node = p.second;
            if(node->probability < first_value){
                second_value = first_value;
                second_index = first_index;

                first_value = node->probability;
                first_index = node->key;
            } else if(node->probability < second_value){
                second_value = node->probability;
                second_index = node->key;
            }// if-else
        }// for

        IG_Binary_Tree* parent = IG_BT_new(first_index, first_value+second_value);
        parent->left = huff_tree[first_index];
        parent->right = huff_tree[second_index];

        huff_tree.erase(first_index);
        huff_tree.erase(second_index);
        huff_tree[parent->key] = parent;
    }// while

    IG_Binary_Tree* ht;
    for(std::pair<int, IG_Binary_Tree*> p : huff_tree)
        ht = p.second;

    // assign binary digits
    IG_BT_assign(ht, NULL, 0, 0);

    // create map <value, binary code>
    IG_BT_create_map(ht);


    std::cout << "SIZE--------------------------------------------------\n";
    std::cout << elements.size();
    std::cout << "SIZE--------------------------------------------------\n";

    std::vector<int> new_data_v;
    IG_Binary_Tree* bt;
    int_to_binary_array(elements.size(), &new_data_v);

    for(std::map<int,int>::iterator iter = elements.begin(); iter != elements.end(); ++iter)
    {
        int_to_binary_array(iter->first, &new_data_v);
        bt = bt_coding[iter->first];

        for(int j = 0 ; j < 12-bt->bincode_size ; j++)
            new_data_v.push_back(0);
        for(int j = 0 ; j < bt->bincode_size ; j++)
            new_data_v.push_back(bt->bincode[j]);
    }

    for(int i = 0 ; i < size ; i++){
        bt = bt_coding[data[i]];
        for(int j = 0 ; j < bt->bincode_size ; j++)
            new_data_v.push_back(bt->bincode[j]);
    }// for

    int je = new_data_v.size();
    free(ig->data);
    ig->data = new int[je];

    for(int i = 0 ; i < je ; i++){
        ig->data[i] = new_data_v.at(i);
    }
    ig->size = je;
}// IG_huffman

void IG_huffman_inv(Image_Gene* ig){
    std::map<std::string,int> mapa;
    std::vector<int> new_data_v;
    std::string kek = "", cod = "";

    int size_val = 18, size_key = 12;

    for(int i=0; i<size_val; i++){
        kek = kek+std::to_string(ig->data[i]);
    }
    int posi, dic_size = std::stoi(kek, nullptr, 2);


    std::cout << "size dessa merda\n";
    std::cout << std::stoull("00000000000000111000001111001011", NULL, 2) << "\n";
    std::cout << dic_size;
    std::cout << "cabô\n";

    for(int j=0; j<dic_size; j++){
        kek = "";
        cod = "";
        for(int i=0; i<size_val; i++){
            posi = size_val + j*(size_val+size_key) + i;
            kek = kek+std::to_string(ig->data[posi]);
        }
        for(int i=0; i<size_key; i++){
            posi = size_val*2 + j*(size_val+size_key) + i;
            cod = cod+std::to_string(ig->data[posi]);
        }
        mapa[std::to_string(std::stoull(cod))] = std::stoull(kek, NULL, 2);
    }
    cod = "";
    for(int i = size_val+dic_size*(size_val+size_key); i<ig->size; i++){
        cod = cod+std::to_string(ig->data[i]);
        if (mapa.count(cod)){
            new_data_v.push_back(mapa[cod]);
            cod = "";
        }
    }

    intArray_to_img_data(new_data_v, ig);

    std::cout << "cabô de vdd agr\n";
}

// ================================================================================================
// ================================================================================================
// MISC ===========================================================================================
// ================================================================================================
// ================================================================================================
void intArray_to_img_data(std::vector<int> new_data_v, Image_Gene* ig){
    int je = new_data_v.size();
    free(ig->data);
    ig->data = new int[je];

    for(int i = 0 ; i < je ; i++){
        ig->data[i] = new_data_v.at(i);
    }
    ig->size = je;
}

void IG_print(Image_Gene* ig){

    std::cout << "----------------------------------------------------------------\n";
    std::cout << "Image:\n";
    std::cout << "Dimensions: " << ig->width << " x " << ig->height << "\n";
    std::cout << "Vector Length: " << ig->size << "\n";
    std::cout << "Current Compression Level: " << ig->compression_level << "\n";
    std::cout << "Max Compression Level: " << ig->max_compression_level << "\n";
    std::cout << "Subdividions: " << ig->subdivisions << "\n";
    std::cout << "Data:\n";

    for(int i = 0 ; i < ig->size ; i++){
        std::cout << (int)ig->data[i] << " ";
    }// for

    std::cout << "\n";
    std::cout << "----------------------------------------------------------------\n";
}// IG_print

// ================================================================================================
// ================================================================================================
// BINARY TREE ====================================================================================
// ================================================================================================
// ================================================================================================

IG_Binary_Tree* IG_BT_new(int key, float p){
    IG_Binary_Tree* tree = new IG_Binary_Tree;
    tree->parent = NULL;
    tree->left = NULL;
    tree->right = NULL;
    tree->probability = p;
    tree->key = key;
    tree->bincode = NULL;

    return tree;
}// IG_BT_new

void IG_BT_print(IG_Binary_Tree* tree,int level){
    for(int i = 0 ; i < level ; i++)
        std::cout << "   ";
    std::cout << "node: " << tree->key << " | " << tree->probability << std::endl;
    for(int i = 0 ; i < level ; i++)
        std::cout << "   ";
    if(level > 0){
        std::cout << "bincode: ";
        for(int i = 0 ; i < tree->bincode_size ; i++){
            std::cout << (int)tree->bincode[i];
        }//for
        std::cout << "\n";
    }//if

    if(tree->left != NULL)
        IG_BT_print(tree->left,level+1);
    if(tree->right != NULL)
        IG_BT_print(tree->right,level+1);
}// IG_BT_print

void IG_BT_print_leafs(IG_Binary_Tree* tree){
    // leaf
    if(tree->left == NULL && tree->right == NULL){
        IG_BT_print(tree, 1);
    } else {
        IG_BT_print_leafs(tree->left);
        IG_BT_print_leafs(tree->right);
    }// if-else
}// IG_BT_print_leafs

void IG_BT_assign(IG_Binary_Tree* tree, char* code, int size, int nextbin){
    // parent
    if(code == NULL){
        tree->bincode = new char[1];
        tree->bincode_size = 0;

        if(tree->left != NULL)
            IG_BT_assign(tree->left, tree->bincode, 0, 0);
        if(tree->right != NULL)
            IG_BT_assign(tree->right, tree->bincode, 0, 1);
    } else {

        tree->bincode = new char[size+1];
        tree->bincode_size = size+1;
        for(int i = 0 ; i < size ; i++)
            tree->bincode[i] = code[i];
        tree->bincode[size] = nextbin;

        if(tree->left != NULL)
            IG_BT_assign(tree->left, tree->bincode, size+1, 0);
        if(tree->right != NULL)
            IG_BT_assign(tree->right, tree->bincode, size+1, 1);
    }// if-else


}// IG_BT_assign

void IG_BT_create_map(IG_Binary_Tree* tree){
    if(tree->left == NULL && tree->right == NULL){
        bt_coding[tree->key] = tree;
    }// if

    if(tree->left != NULL)
        IG_BT_create_map(tree->left);
    if(tree->right != NULL)
        IG_BT_create_map(tree->right);
}// IG_BT_create_map
