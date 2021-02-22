#pragma once
    // Might need to change this
#include <tuple>

using namespace std;


#define MIN_NUM -pow(1.17, 38)
#define num double
#define pos unsigned


bool Indexable(int i, pos len){
    if( i > len || i < -len) {
        // The negative condition it to make sure out of bounds index is not an issue 
            // while still allowing for negative indexing
        throw invalid_argument("Indexed outside bounds of array.");
    }
    return true;
}


struct Vector {
    num* array;
    pos len;
    Vector(pos len){
        length = len;
        array = new num[len];
    }
    num &operator[](int i) {
        // Add nice indexing
        if (Indexable(i, len)){
            if ( i < 0 ){
                i = len - i;
                // You do not need a -1 here because the negative indexing starts at -1
            }
            return array[i];
        }
      }
};


struct Matrix {
    num** array;
    pos shape[2];
        // We may need to do something here to return shape nicely
            // By reference vs by value
    Matrix(pos rows, pos cols){
        shape[0] = rows;
        shape[1] = cols;
    }
    Vector &operator[](int i){
        // 2D indexing just happens because the 2nd index triggers on the vector that gets returned
        if(Indexable(i, shape[0])) {
            if ( i < 0 ){
                i = len - i;
                // You do not need a -1 here because the negative indexing starts at -1
            }
            return array[i]; 
         }
    }
};