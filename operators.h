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
	for (unsigned i = 0; i < lenA; i++) {
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
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = a[i] - num_b;
	}
	return result;
}
template<typename T>
// This should be true false but with a 0/1 array of ints it'll be fine
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
template<typename T>
Vector operator == (Vector &a, T value) {
	pos lenA = a.len;
	Vector toReturn(lenA);
	for (unsigned i = 0; i < lenA; i++) {
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


// Matrix Operators







