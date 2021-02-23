#include "gloabls"
#pragma once

using namespace std;

template<typename T, typename V>
Vector VectorFull(T len, V value) {
	Vector toReturn;
	if (len < 0.0) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Size of the vector must be positive."); }
	for (unsigned i = 0; i<unsigned(len); i++) { toReturn.push_back(num(value)); }
	return toReturn;
}
template<typename T, typename V>
Vector rep(T len, V value) {
    return VectorFull(len, value);
}
template<typename T, typename V>
Vector full(T len, V val) {
    return VectorFull(len, val);
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
template<typename T, typename V>
Matrix rep(T shape[2], V val) {
    return MatrixFull(shape[0], shape[1], val);
}
Matrix full(T shape[2], V val) {
    return MatrixFull(shape[0], shape[1], val);
}
template<typename T>
Matrix diag(Matrix &toFill, T val){
    pos rows = toFill.shape[0];
    if(rows != toFill.shape[1]){throw invalid_argument("diag must be called on a square matrix.");}
    for(pos i=0; i<rows; i++){
        toFill[i][i] = val
    }
}
Vector seq(num lower, num upper, pos steps) {
	Vector toReturn = new num[steps];
	for (int i = 0; i < steps; i++) {
		toReturn[i] = lower + (i * ((upper - lower) / steps));
	}
	return toReturn;
}