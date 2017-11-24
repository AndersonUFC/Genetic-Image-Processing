#include "gip.h"

#include <iostream>
#include <cmath>

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


void IG_haar2D(Image_Gene* ig){

    int width_bound = pow(2, ig->compression_level);
    int* lineR = new int[width_bound];
    int* lineG = new int[width_bound];
    int* lineB = new int[width_bound];

    int height = ig->height;
    int width = ig->width;
    int* data = ig->data;

    int g_index = ig->size/3;
    int b_index = ((2*ig->size)/3);

    // LINES
    for(int y = 0 ; y < height ; y++){

        /* RGB sequential mode
        // get line
        for(int x = 0 ; x < width_bound ; x++){
            int position = (y*width + x)*3;

            //std::cout << position << " ";

            lineR[x] = data[position];
            lineG[x] = data[position+1];
            lineB[x] = data[position+2];
        }// for-x

        int* cR = IG_haar1D(lineR,width_bound,ig->compression_level);
        int* cG = IG_haar1D(lineG,width_bound,ig->compression_level);
        int* cB = IG_haar1D(lineB,width_bound,ig->compression_level);


        for(int x = 0 ; x < width_bound ; x++){
            int position = (y*width + x)*3;
            data[position] = cR[x];
            data[position+1] = cG[x];
            data[position+2] = cB[x];
        }// for-x
        */

        // R sequence, G sequence, B sequence
        // get line
        for(int x = 0 ; x < width_bound ; x++){
            int position = y*width + x;

            lineR[x] = data[position];
            lineG[x] = data[position+g_index];
            lineB[x] = data[position+b_index];
        }// for-x


        int* cR = IG_haar1D(lineR,width_bound,ig->compression_level);
        int* cG = IG_haar1D(lineG,width_bound,ig->compression_level);
        int* cB = IG_haar1D(lineB,width_bound,ig->compression_level);

        for(int x = 0 ; x < width_bound ; x++){
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

    int height_bound = width_bound;
    // COLUMNS
    for(int x = 0 ; x < width ; x++){

        /*
        // get column
        for(int y = 0 ; y < height_bound ; y++){
            int position = (y*width + x)*3;

            lineR[y] = data[position];
            lineG[y] = data[position+1];
            lineB[y] = data[position+2];
        }// for-x

        int* cR = IG_haar1D(lineR,width_bound,ig->compression_level);
        int* cG = IG_haar1D(lineG,width_bound,ig->compression_level);
        int* cB = IG_haar1D(lineB,width_bound,ig->compression_level);

        for(int y = 0 ; y < height_bound; y++){
            int position = (y*width + x)*3;
            data[position] = cR[y];
            data[position+1] = cG[y];
            data[position+2] = cB[y];
        }// for-x
        */

        // get column
        for(int y = 0 ; y < height_bound ; y++){
            int position = y*width + x;

            lineR[y] = data[position];
            lineG[y] = data[position+g_index];
            lineB[y] = data[position+b_index];
        }// for-x

        int* cR = IG_haar1D(lineR,width_bound,ig->compression_level);
        int* cG = IG_haar1D(lineG,width_bound,ig->compression_level);
        int* cB = IG_haar1D(lineB,width_bound,ig->compression_level);

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

// ------------------------------------------------------------------------------------------------

void IG_haar1D_subdivide(Image_Gene* ig){
    IG_print(ig);

}// IG_haar1D_subdivide

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
    /*
    for(int i = 0 ; i < 3 ; i++)
        new_data[i] = data[i];
    for(int i = 3 ; i < size ; i += 3){
        new_data[i] = data[i] - data[i-3];
        new_data[i+1] = data[i+1] - data[i-2];
        new_data[i+2] = data[i+2] - data[i-1];
    }// for
    */
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

    std::vector<int> new_data_v;
    IG_Binary_Tree* bt;
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

// ================================================================================================
// ================================================================================================
// MISC ===========================================================================================
// ================================================================================================
// ================================================================================================

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
