#include "globals.h"

using namespace std;

//Vector Operators
Vector operator - (Vector a, Vector b) {
	pos lenA = a.size();
	if (lenA != b.size()) {
		cout << "ERROR ENCOUNTERED" << endl; 
        throw invalid_argument("Vectors must have the same dimensions.");
	}
	Vector result = new Vector[lenA];
	for (pos i = 0; i < lenA; i++) {
		result[i] = a[i] - b[i];
	}
	return result;
}
Vector operator - (Vector a) {
	pos lenA = a.size();
	Vector result = new Vector[lenA];
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = -a[i];
	}
	return result;
}
template <typename T>
Vector operator - (Vector a, T b) {
	pos lenA = a.size();
	Vector result = new Vector[lenA];
	num num_b = num(b);
		// Verify it is a num so c++ doesn't throw run time errors
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = a[i] - num_b;
	}
	return result;
}
template<typename T>
// This should be true false but with a 0/1 array of ints it'll be fine
Vector operator > (Vector a, T b) {
	pos lenA = a.size();
	num num_b = num(b);
	vector greaterThan = new vector[lenA];
	for (pos i = 0; i < lenA; i++) {
		if (a[i] > num_b) { a[i] = 0; }
		else{ a[i] = 1; }
	}
	return greaterThan;
}
template<typename T>
vector<int> operator < (T b, Vector a) {
	// Hopefully this isn't too slow
	return a > b;
}











