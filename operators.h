#pragma once
#include "globals.h"
#include <iostream>

using namespace std;

//Vector Operators
Vector operator - (Vector &a, Vector &b) {
	pos lenA = a.len;
	if (lenA != b.len) {
		cout << "ERROR ENCOUNTERED" << endl; 
        throw invalid_argument("Vectors must have the same dimensions.");
	}
	Vector result(lenA);
	for (pos i = 0; i < lenA; i++) {
		result[i] = a[i] - b[i];
	}
	return result;
}
Vector operator - (Vector &a) {
	pos lenA = a.len;
	Vector result(lenA);
	for (pos i = 0; i < lenA; i++) {
		result[i] = -a[i];
	}
	return result;
}
template <typename T>
Vector operator - (Vector &a, T b) {
	pos lenA = a.len;
	Vector result(lenA);
	num num_b = num(b);
		// Verify it is a num so c++ doesn't throw run time errors
	for (pos i = 0; i < lenA; i++) {
		result[i] = a[i] - num_b;
	}
	return result;
}
Vector operator + (Vector &a, Vector &b) {
	pos lenA = a.len;
	if (lenA != b.len) {
		cout << "ERROR ENCOUNTERED" << endl; 
        throw invalid_argument("Vectors must have the same dimensions.");
	}
	Vector result(lenA);
	for (pos i = 0; i < lenA; i++) {
		result[i] = a[i] + b[i];
	}
	return result;
}
template <typename T>
Vector operator + (Vector &a, T b) {
	pos lenA = a.len;
	Vector result(lenA);
	num num_b = num(b);
		// Verify it is a num so c++ doesn't throw run time errors
	for (pos i = 0; i < lenA; i++) {
		result[i] = a[i] + num_b;
	}
	return result;
}
template<typename T>
Vector operator > (Vector &a, T b) {
	pos lenA = a.size();
	num num_b = num(b);
	Vector greaterThan(lenA);
	for (pos i = 0; i < lenA; i++) {
		if (a[i] > num_b) { a[i] = 0; }
		else{ a[i] = 1; }
	}
	return greaterThan;
}
template<typename T>
Vector operator < (T b, Vector &a) {
	// Hopefully this isn't too slow
	return a > b;
}
Vector operator >= (Vector &a, T b) {
	pos lenA = a.size();
	num num_b = num(b);
	Vector greaterThan(lenA);
	for (pos i = 0; i < lenA; i++) {
		if (a[i] => num_b) { a[i] = 0; }
		else{ a[i] = 1; }
	}
	return greaterThan;
}
template<typename T>
Vector operator <= (T b, Vector &a) {
	// Hopefully this isn't too slow
	return a => b;
}
template<typename T>
Vector operator < (Vector &a, T b) {
	pos lenA = a.len;
	num num_b = num(b);
	Vector greaterThan(lenA);
	for (pos i = 0; i < lenA; i++) {
		if (a[i] < num_b) { a[i] = 0; }
		else{ a[i] = 1; }
	}
	return greaterThan;
}
template<typename T>
Vector operator > (T b, Vector &a) {
	// Hopefully this isn't too slow
	return a < b;
}
Vector operator <= (Vector &a, T b) {
	pos lenA = a.len;
	num num_b = num(b);
	Vector greaterThan(lenA);
	for (pos i = 0; i < lenA; i++) {
		if (a[i] <= num_b) { a[i] = 0; }
		else{ a[i] = 1; }
	}
	return greaterThan;
}
template<typename T>
Vector operator >= (T b, Vector &a) {
	// Hopefully this isn't too slow
	return a <= b;
}
template <typename T>
Vector operator * (Vector &a, T b) {
	pos lenA = a.len;
	num val = num(b);
	Vector result(lenA);
		// Lets hope that this works
	for (pos i = 0; i < lenA; i++) {
		result[i] = a[i] * val;
	}
	return result;
}
template<typename T>
Vector operator * (T b, Vector &a) {
	return a * b;
}
Vector operator * (Vector &a, Vector &b) {
	pos lenA = a.len;
	pos lenB = b.len;
	if (lenA != lenB) { 
		cout << "ERROR OCCURRED. DIMS DONT MATCH" << endl; 
		throw invalid_argument("Dimensions don't match in vector multiplication."); 
	}
	Vector toReturn(lenA);
	for (int i = 0; i < lenA; i++) {
		toReturn[i] = a[i] * b[i];
	}
	return toReturn;
}
template<typename T>
Vector operator / (Vector &a, T b) {
	pos lenA = a.len;
	num val = num(b);
	Vector toReturn(lenA);
	for (pos i = 0; i < lenA; i++) {
		toReturn[i] = a[i] / val;
	}
	return toReturn;
}
template<typename T>
Vector operator / (T b, Vector &a) {
	pos lenA = a.len;
	num val = num(b);
	Vector toReturn(lenA);
	for (pos i = 0; i < lenA; i++) {
		toReturn[i] = val / a[i];
	}
	return toReturn;
}
Vector operator / (Vector &a, Vector &b){
	pos lenA = a.shape;
	if ( lenA != b.shape ) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Vectors must have the same dimensions."); }
	Vector toReturn(lenA);
	for (pos i=0; i<lenA; i++){
		toReturn[i] = a[i] / b[i];
	}
	return toReturn;
}
template<typename T>
Vector operator == (Vector &a, T value) {
	pos lenA = a.len;
	Vector toReturn(lenA);
	for (pos i = 0; i < lenA; i++) {
		if (a[i] == value) { toReturn[i] = 1; }
		else { toReturn[i] = 0; }
	}
	return toReturn;
}
template<typename T>
Vector operator == (T value, Vector &a) {
	return a == value;
}

bool operator == (Vector &a, Vector &b) {
	pos lenA = a.len;
	if(lenA != b.len){ throw invalid_argument("Vectors must have the same dimensions."); }
	for( int i=0; i<lenA; i++){
		if (a[i]!=b[i]){return false;}
	}
	return true;
}
template<typename T>
Vector operator % (Vector &a, T b) {
	pos lenA = a.shape;
	Vector toReturn(lenA);
	for(pos i=0; i<lenA; i++){ toReturn[i] = a[i] % b; }
	return toReturn;
}
template<typename T>
Vector operator % (T b, Vector &a){
	pos lenA = a.len;
	Vector toReturn(lenA);
	for(pos i=0; i<lenA; i++){
		toReturn[i] = b % a[i];
	}
}
Vector operator % (Vector &a, Vector &b) {
	pos lenA = a.shape;
	if(lenA != b.len){ throw invalid_argument("Vectors must have the same dimensions."); }
	Vector toReturn(lenA);
	for(pos i=0; i<lenA; i++){
		toReturn[i] = a[i] % b[i];
	}
	return toReturn;
}


// Matrix Operators

Matrix operator - (Matrix &a, Matrix &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos rowB = b.shape[0];
	pos colB = b.shape[1];
	if (rowA != rowB || colA != colB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Matricies must have the same dimensions to be subtracted."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = a[i][j] - b[i][j];
		}
	}
	return result;
}
Matrix operator - (Matrix &a, Vector &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (rowA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns for the matrix and vector must be the same."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = a[i][j] - b[i];
		}
	}
	return result;
}
Matrix operator - (Vector &b, Matrix &a) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (rowA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns for the matrix and vector must be the same."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = b[i] - a[i][j];
		}
	}
	return result;
}
Matrix operator - (Matrix &a) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = -a[i][j];
		}
	}
	return result;
}
template <typename T>
Matrix operator - (Matrix &a, T b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num temp = num(b);
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = a[i][j] - temp;
		}
	}
	return result;
}
template <typename T>
Matrix operator - (T b, Matrix &a) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num temp = num(b);
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = temp - a[i][j];
		}
	}
	return result;
}
Matrix operator + (Matrix &a, Matrix &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos rowB = b.shape[0];
	pos colB = b.shape[1];
	if (rowA != rowB || colA != colB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Matricies must have the same dimensions to be subtracted."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = a[i][j] + b[i][j];
		}
	}
	return result;
}
Matrix operator + (Matrix &a, Vector &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (rowA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns for the matrix and vector must be the same."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = a[i][j] + b[i];
		}
	}
	return result;
}
Matrix operator + (Vector &b, Matrix &a) {
	return a + b;
}
template <typename T>
Matrix operator + (Matrix &a, T b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num temp = num(b);
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i][j] = a[i][j] + temp;
		}
	}
	return result;
}
template <typename T>
Matrix operator + (T b, Matrix &a) {
	return a + b;
}
template<typename T>
Matrix operator < (Matrix &a, T b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num val = num(b);
	Matrix toReturn(rowA, colA)
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			if (a[i][j] < val) { toReturn[i][j] = 1; }
			else { toReturn[i][j] = 0; }
		}
	}
	return toReturn;
}
template<typename T>
Matrix operator > (T b, Matrix &a){
	return a < b;
}
template<typename T>
vector<tuple<pos, pos>> operator > (Matrix &a, T b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num val = num(b);
	Matrix toReturn(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			if (a[i][j] > val) { toReturn[i][j] = 1; }
			else { toReturn[i][j] = 0; }
		}
	}
	return toReturn;
}
template<typename T>
Matrix operator < (T b, Matrix &a){
	return a > b;
}
template<typename T>
Matrix operator <= (Matrix &a, T b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num val = num(b);
	Matrix toReturn(rowA, colA)
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			if (a[i][j] <= val) { toReturn[i][j] = 1; }
			else { toReturn[i][j] = 0; }
		}
	}
	return toReturn;
}
template<typename T>
Matrix operator >= (T b, Matrix &a){
	return a <= b;
}
template<typename T>
vector<tuple<pos, pos>> operator >= (Matrix &a, T b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num val = num(b);
	Matrix toReturn(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			if (a[i][j] >= val) { toReturn[i][j] = 1; }
			else { toReturn[i][j] = 0; }
		}
	}
	return toReturn;
}
template<typename T>
Matrix operator <= (T b, Matrix &a){
	return a >= b;
}
Matrix operator * (Matrix &a, Matrix &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos rowB = b.shape[0];
	pos colB = b.shape[1];
	Matrix result(rowA, colB);
	if (colA != colB || rowA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("Matricies were multiplied but not of the same shape."); }
	for (pos i = 0; i < colB; i++) {
		for (pos j = 0; j < colB; j++) {
			result[i][j] = a[i][j] * b[i][j];
		}
	}
	return result;
}
Matrix operator * (Matrix &a, Vector &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (colA != lenB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The column size of Matrix and the length of Vector must be the same."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			result[i] = a[i][j] * b[j];
		}
	}
	return result;
}
Matrix operator * (Vector &b, Matrix &a) {
	return a * b;
}
template<typename T> 
Matrix operator * (Matrix &a, T b){
	num val = num(b);
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	Matrix toReturn(rowA, colA);
	for(pos i=0; i<rowA; i++){
		for(pos j=0; j<colA; j++){
			toReturn[i][j] = a[i][j] * val;
		}
	}
	return toReturn;
}
template<typename T> 
Matrix operator * (T b, Matrix &a){
	return a * b;
}
Matrix operator / (Matrix &a, Vector &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (colA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns must line up for matrix, vector division."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < colA; i++) {
		for (pos j = 0; j < rowA; j++) {
			result[i][j] = a[i][j] / b[j];
		}
	}
	return result;
}
Matrix operator / (Vector &b, Matrix &a) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (colA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns must line up for matrix, vector division."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < colA; i++) {
		for (pos j = 0; j < rowA; j++) {
			result[i][j] = b[j] / a[i][j];
		}
	}
	return result;
}
Matrix operator / (Matrix &a, Matrix &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos rowB = b.shape[0];
	pos colB = b.shape[1];
	Matrix result(rowA, colB);
	if (colA != colB || rowA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("Matricies were multiplied but not of the same shape."); }
	for (pos i = 0; i < colB; i++) {
		for (pos j = 0; j < colB; j++) {
			result[i][j] = a[i][j] / b[i][j];
		}
	}
	return result;
	
}
template<typename T> 
Matrix operator / (Matrix &a, T b){
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num val = num(b);
	Matrix toReturn(rowA, colA);
	for(pos i=0; i<rowA; i++){
		for(pos j=0; j<colB; j++){
			toReturn[i][j] = a[i][j] / val;
		}
	}
	return toReturn;
}
template<typename T> 
Matrix operator / (T b, Matrix &a){
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	num val = num(b);
	Matrix toReturn(rowA, colA);
	for(pos i=0; i<rowA; i++){
		for(pos j=0; j<colB; j++){
			toReturn[i][j] = val / a[i][j];
		}
	}
	return toReturn;
}
template<typename T>
bool operator == (Matrix &a, T value) {
	num rowA = a.shape[0];
	num colA = a.shape[1];
	num val = num(value);
	Matrix toReturn(rowA, colA);
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			if (a[i][j] == val) { toReturn[i][j] = 1; }
			else { toReturn[i][j] = 0; }
		}
	}
	return toReturn;
}
template<typename T>
bool operator == (Matrix a, T value) {
	return a == value;
}
bool operator == (Matrix &a, Matrix &b) {
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos rowB = b.shape[0];
	pos colB = b.shape[1];
	Matrix result(rowA, colB);
	if (colA != colB || rowA != rowB) { throw invalid_argument("Matricies were multiplied but not of the same shape."); }
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			if (a[i][j] != b[i][j]) { return false; }
		}
	}
	return true;
}
template<typename T>
Matrix operator % (Matrix &a, T b){
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	Matrix toReturn(rowA, colA);
	num val = num(b);
	for(pos i=0; i<rowA; i++){
		for(pos j=0; j<colA; j++){
			toReturn[i][j] = a[i][j] % val;
		}
	}
}
template<typename T>
Matrix operator % (T b, Matrix &a){
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	Matrix toReturn(rowA, colA);
	num val = num(b);
	for(pos i=0; i<rowA; i++){
		for(pos j=0; j<colA; j++){
			toReturn[i][j] = val % a[i][j];
		}
	}
}
Matrix operator % (Matrix &a, Vector &b){
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (colA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns must line up for matrix, vector division."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < colA; i++) {
		for (pos j = 0; j < rowA; j++) {
			result[i][j] = a[i][j] % b[j];
		}
	}
	return result;
}
Matrix operator % (Vector &b, Matrix &a){
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos lenB = b.shape;
	if (colA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns must line up for matrix, vector division."); }
	Matrix result(rowA, colA);
	for (pos i = 0; i < colA; i++) {
		for (pos j = 0; j < rowA; j++) {
			result[i][j] = b[j] % a[i][j];
		}
	}
	return result;
}
Matrix operator % (Matrix &a, Matrix &b){
	pos rowA = a.shape[0];
	pos colA = a.shape[1];
	pos rowB = b.shape[0];
	pos colB = b.shape[1];
	Matrix result(rowA, colB);
	if (colA != colB || rowA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("Matricies were multiplied but not of the same shape."); }
	for (pos i = 0; i < colB; i++) {
		for (pos j = 0; j < colB; j++) {
			result[i][j] = a[i][j] % b[i][j];
		}
	}
	return result;
}
