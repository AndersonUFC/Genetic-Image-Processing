#ifndef GIP_H
#define GIP_H

#include <stdlib.h>
#include <vector>
#include <map>

// DATA TYPE ======================================================================================
typedef struct Image_Gene{
    // DATA INFO
    /**
     * @brief size: size of vector data
     * @brief height: height of image
     * @biref width: width of image
     */
    int size, height, width;

    /**
     * @brief data: vector of type [R1,G1,B1,R2,G2,B2,...,Rn,Gn,Bn]
     */
    int* data;

    // COMPRESSION INFO
    /**
     * @brief compression_level: level of compression used in Haar transformation
     */
    int compression_level, max_compression_level, subdivisions;
}Image_Gene;

typedef struct Image_Gene_Float{
    // DATA INFO
    /**
     * @brief size: size of vector data
     * @brief height: height of image
     * @biref width: width of image
     */
    int size, height, width;

    /**
     * @brief data: vector of type [R1,G1,B1,R2,G2,B2,...,Rn,Gn,Bn]
     */
    float* data;

    // COMPRESSION INFO
    /**
     * @brief compression_level: level of compression used in Haar transformation
     */
    int compression_level, max_compression_level, subdivisions;
}Image_Gene_Float;

// CONSTRUCTOR ====================================================================================
Image_Gene* IG_new(unsigned char*, int, int, int);

// CONVERSION =====================================================================================
void data_to_RGB(unsigned char*, int, unsigned char**, unsigned char**, unsigned char**);
void configure_data(int*,int, int**);

Image_Gene_Float* IG_int_to_float(Image_Gene*);
Image_Gene* IG_float_to_int(Image_Gene_Float*);

// TRANSFORMATION =================================================================================
int* IG_haar1D(int*, int, int);
int* IG_haar1D_inv(int*, int, int);

void IG_haar2D(Image_Gene*);

void IG_haar2D_subdivide(Image_Gene*);
void IG_haar2D_inv_subdivide(Image_Gene*);

// -----------------------------------------------------------

float* IG_haar1D_float(float*, int, int);
float* IG_haar1D_float_inv(float*, int, int);

void IG_haar2D_subdivide_float(Image_Gene_Float*);
void IG_haar2D_inv_subdivide_float(Image_Gene_Float*);

// COMPRESSION ====================================================================================
void IG_predictive(Image_Gene*);
void IG_run_length(Image_Gene*);

void IG_huffman(Image_Gene*);
void IG_huffman_inv(Image_Gene*);

// MISC ===========================================================================================
void IG_print(Image_Gene*);

// BINARY TREE ====================================================================================

typedef struct IG_Binary_Tree{
    struct IG_Binary_Tree *parent, *left, *right;
    float probability;
    int key, bincode_size;
    char* bincode;
}IG_Binary_Tree;

IG_Binary_Tree* IG_BT_new(int,float);

void IG_BT_print(IG_Binary_Tree*,int);
void IG_BT_print_leafs(IG_Binary_Tree*);

void IG_BT_assign(IG_Binary_Tree*, char*, int, int);
void IG_BT_create_map(IG_Binary_Tree*);

#endif // GIP_H
