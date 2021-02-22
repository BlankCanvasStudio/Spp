#pragma once
    // Might need to change this
#include <tuple>

using namespace std;


#define MIN_NUM -pow(1.17, 38)
#define num double
#define pos unsigned
#define let auto
    // Too much js



bool Indexable(int i, pos len){
    if( i > len || i < -len) {
        // The negative condition it to make sure out of bounds index is not an issue 
            // while still allowing for negative indexing
        throw invalid_argument("Indexed outside bounds of array.");
        return false;
    }
    return true;
}


struct Vector {
    num* array;
    pos len;
    pos shape;
        // Unnecessary but a nice touch
    Vector(pos length){
        len = length;
        shape = length;
        array = new num[len];
    }
    Vector(num* array_in, pos Length){
        len = Length;
        shape = Length;
        array = new num[len];
        // We could do a pointer swap here but I'm worried that user deletion might interefere 
            // With the logical working of vectors so I went with a copy.
        for( int i=0; i<len; i++){
            array[i] = array_in[i];
        } 
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
        shape = [rows, cols];
        array = new num*[rows];
        for(int i=0; i<rows; i++){
            array[i] = new num[cols];
        }
    }
    Matrix(num** array_in, pos rows, pos cols){
        shape = [rows, cols];
        array = new num*[rows];
        for(int i=0; i<rows; i++){
            array[i] = new num[cols];
            for(int j=0;j<cols;j++){
                array[i][kj] = array_in[i][j];
            }
        }
    }
    Vector operator[](int i){
        // 2D indexing just happens because the 2nd index triggers on the vector that gets returned
        if(Indexable(i, shape[0])) {
            if ( i < 0 ){
                i = shape[0] - i;
                // You do not need a -1 here because the negative indexing starts at -1
            }
            return Vector(array[i], shape[0]);
                // Quick little typecast might want to spped up typecast so this doesn't slow the prog
         }
    }
};


// Some Helpful Constructors 

template<typename T, typename V>
Vector VectorFull(T len, V value) {
	Vector toReturn;
	if (len < 0.0) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Size of the vector must be positive."); }
	for (unsigned i = 0; i<unsigned(len); i++) { toReturn.push_back(num(value)); }
	return toReturn;
}

template<typename T, typename V, typename C>
Matrix MatrixFull(T row, V col, C val) {
    num value = num(val);
    Matrix toReturn(row, col);
    for( pos i=0; i<row; i++) {
        for( pos j=0; j<col; j++) {
            toReturn[i][j] = value;
        }
    }
    return toReturn;
}