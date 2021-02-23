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
        return false;
    }
    return true;
}


struct Vector {
    num* array;
    pos len;
    pos shape;
        // Unnecessary but a nice touch
    template<typename T>
    Vector(T length){
        len = pos(length);
        shape = pos(length);
        array = new num[len];
    }
    template<typename T>
    Vector(num* array_in, T Length){
        len = pos(Length);
        shape = pos(Length);
        array = new num[len];
        // We could do a pointer swap here but I'm worried that user deletion might interefere 
            // With the logical working of vectors so I went with a copy.
        for( int i=0; i<len; i++){
            array[i] = array_in[i];
        } 
    }
    Vector(IndexVector &a){
        len = a.len;
        shape = a.shape;
        array = new num[len];
        for(pos i=0; i<len; i++) { array[i] = a[i]; }
            // This works since we overloaded the indexing operator
    }
    num &operator [] (int i) {
        // Add nice indexing
        if (Indexable(i, len)){
            if ( i < 0 ){
                i = len - i;
                // You do not need a -1 here because the negative indexing starts at -1
            }
            return array[i];
        }
      }
    }
    Vector &operator [] (Vector &b){
        // This will allow us to return the indexes where the condition is true
        pos lenB = b.shape;
        if (lenB != len){ throw invalid_argument("Vectors must have the same dimensions."); }
        pos num_true = 0;
        for(pos i=0; i<lenB; i++){
            if(b[i]){ num_true++; }
        }
        // I'm going to claim that running through to see how long the array should be will be shorter 
            // Than dynamically allocating thr array every time so we do it this way
        Vector toReturn(num_true);
        pos index = 0;
        for(pos i=0; i<lenB; i++){
            if(b[i]){ toReturn[index] = i; index++; }
        }
        return toReturn;
    }
    // This returns index vector so you can do Array(Array>2) = 3
    IndexVector operator () (Vector &b){
        pos lenB = b.shape;
        if (lenB != len){ throw invalid_argument("Vectors must have the same dimensions."); }
        pos num_true = 0;
        for(pos i=0; i<lenB; i++){
            if(b[i]){ num_true++; }
        }
        // I'm going to claim that running through to see how long the array should be will be shorter 
            // Than dynamically allocating thr array every time so we do it this way
        IndexVector toReturn(num_true);
        pos index = 0;
        for(pos i=0; i<lenB; i++){
            if(b[i]){ toReturn[index] = &array[i]; index++; }
        }
        return toReturn;
    }
    ~Vector() {
        delete[] array;
    }
};

// I'm trying to implement foo[bar>4] = 3 and foo[bar>4]
    // One requires pointers, the other requires index number so I'm creating a struct 
    // To hold pointers to the indexes
// This is a fairly beautiful solution. Since we overloaded the index operator it will only ever refer to the value pointed to 
    // at the index. So when you set it, it will set the value at the index. When you read from it, it will read from the 
    // value at the index. Thus, it will behave just like a vector. We may get some reference/value errors when setting 
    // this array later but we have an overload in vector to convert it so hopefully that'll work
struct IndexVector: public Vector {
    num** indicies;
        // This is going to be an array of pointer made dynamically thus a pointer to a pointer 
    pos len;
    pos shape; 
    IndexVector(pos length) {
        len = length;
        shape = length;
    }
    IndexVector(num** index_in, pos length) {
        indicies = new *num[length];
        indicies = index_in;
        len = length;
        shape = length;
    }
    ~IndexVector() {
        delete[] indicies;
            // This should only ever be 1D array but lets allow for expansion if necessary
    }
    num operator [] (int i){
        if (Indexable(i, len)){
            if ( i < 0 ){
                i = len - i;
                // You do not need a -1 here because the negative indexing starts at -1
            }
            return *(indicies[i]);
        }
    }
};


struct Matrix {
    num** array;
    pos shape[2];
        // We may need to do something here to return shape nicely
            // By reference vs by value
    template<typename T, typename C>
    Matrix(T rows, C cols){
        shape = [pos(rows), pos(cols)];
        array = new num*[shape[0]];
        for(int i=0; i<rows; i++){
            array[i] = new num[shape[1]];
        }
    }
    template<typename T>
    Matrix(T in_shape[2]){
        shape = [pos(in_shape[0]), pos(in_shape[1])];
        array = new num*[shape[0]];
        for(int i=0; i<rows; i++){
            array[i] = new num[shape[1]];
        }
    }
    template<typename T, typename C>
    Matrix(num** array_in, T rows, C cols){
        shape = [pos(rows), pos(cols)];
        array = new num*[shape[0]];
        for(int i=0; i<shape[0]; i++){
            array[i] = new num[shape[1]];
            for(int j=0;j<shape[1];j++){
                array[i][kj] = array_in[i][j];
            }
        }
    }
    Matrix(IndexMatrix in_mat){
        shape = in_mat.shape;
        array = new num*[shape[0]];
        for(int i=0; i<shape[0]; i++){
            array[i] = new num[shape[1]];
            for(int j=0;j<shape[1];j++){
                array[i][kj] = array_in[i][j];
            }
        }
    }
    ~Matrix() {
        for(int i=0; i<shape[0]; i++){
            delete[] array[i];
        }
    }
    Vector &operator[](int i){
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
    Matrix &operator [] (Matrix &b){
        
        if (b.shape != shape){ throw invalid_argument("Vectors must have the same dimensions."); }
        Matrix toReturn = new num*[shape[0]];
        // Piggies backing off the Vector function for ease of use
        for(pos i=0; i<shape[0]; i++){
            toReturn[i] = array[i][b[i]];
                //This is indexing a vector using another vector
        }
        return toReturn;
    }
    // This returns index vector so you can do Array(Array>2) = 3
    IndexMatrix operator () (Vector &b){
        // This will allow us to return the indexes where the condition is true
        
        if (b.shape != shape){ throw invalid_argument("Vectors must have the same dimensions."); }
        IndexMatrix toReturn(shape[0]);
        // Piggies backing off the Vector function for ease of use
        for(pos i=0; i<shape[0]; i++){
            toReturn[i] = array[i](b[i]);
            if(toReturn.shape[1] < toReturn[i].shape){
                toReturn.shape[1] = toReturn[i].shape;
                // Just so that shape represents the largest space in the matrix
            }
                //This is indexing a vector using another vector using IndexVector operator
        }
        return toReturn;
    }
};

// Hopefully piggybacking off the IndexVector stuff doesn't slow us down but its nice coding wise
struct IndexMatrix: public Matrix {
    IndexVector *indicies;
    pos shape[2];

    template<typename T>
    IndexMatrix(T rows){
        shape = [pos(rows), 0];
        indicies = new IndexVector[shape[0]];
    }
    template<typename T, typename C>
    IndexMatrix(T rows, C cols){
        shape = [pos(rows), pos(cols)];
        indicies = new IndexVector[shape[0]];
        for (pos i=0; i<shape[0]; i++){
            indicies[i] = new IndexVector(shape[1]);
        }
    }
    IndexMatrix(pos in_shape[2]){
        shape = in_shape;
        indicies = new IndexVector[shape[0]];
        for (pos i=0; i<shape[0]; i++){
            indicies[i] = new num*[shape[1]];
        }
    }
    IndexMatrix(num** index_in, pos in_shape[2]){
        shape = in_shape;
        indicies = new IndexVector[shape[0]];
        for (pos i=0; i<shape[0]; i++){
            indicies[i] = new IndexVector(shape[1]);
            for (pos j=0; j<shape[1]; j++){
                indicies[i][j] = index_in[i][j];
            }
        }
    }
    template<typename T, typename C>
    IndexMatrix(num** index_in, T rows, C cols){
        shape = [pos(rows), pos(cols)];
        indicies = new IndexVector[shape[0]];
        for (pos i=0; i<shape[0]; i++){
            indicies[i] = new IndexVector(shape[1]);
            for (pos j=0; j<shape[1]; j++){
                indicies[i][j] = index_in[i][j];
            }
        }
    }
    ~IndexMatrix() {
        for(pos i=0; i<shape[0]; i++){
            delete indicies[i];
        }
    }
}
