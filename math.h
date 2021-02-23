#include "gloabls.h"
#include <math.h>

using namespace std;

// Vector Functions 

num sum(Vector &a) {
	num toReturn = 0;
    pos size = a.shape;
	for (pos i = 0; i < size; i++) {
		toReturn += a[i];
	}
	return toReturn;
}
num mean(Vector &a) {
	return sum(a) / num(a.shape);
}
num Variance(Vector &a) {
    pos size = a.shape;
	num meanA = mean(a, size);
	num total = 0.0;
	for (int i = 0; i < size; i++) {
		total += pow((a[i] - meanA), 2);
	}
	if (size - 1 == 0) { return 0; }
	return total / num(size - 1);
}
num std_dev(Vector &a) {
	return sqrt(Variance(a, a.shape));
}
num mad(Vector &toMad) {
    pos len = toMad.len;
	num absSum = 0;
	num meanVal = mean(toMad, len);
	for (int i = 0; i < len; i++) {
		absSum += abs(toMad[i] - meanVal);
	}
	return absSum / len;
}
Vector abs(Vector &a) {
    pos size = a.shape;
    Vector toReturn(size);
    for(int i=0; i<size; i++){
        toReturn[i] = abs(a[i]);
    }
    return toReturn;
}
template <typename T>
Vector pow(Vector &a, T degree) {
	pos lenA = a.shape;
	Vector result(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = pow(a[i], degree);
	}
	return result;
}
Vector log(Vector &a) {
	pos lenA = a.shape;
	Vector toReturn(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		toReturn[i] = log(a[i]);
	}
	return toReturn;
}




// Matrix Functions

Vector mean(Matrix toMean) {
    pos n = toMean.shape[0];
    pos p = toMean.shape[1];
	Vector toReturn(p);
	for (int i = 0; i < p; i++) {
		Vector temp = new num[n];
		for (int j = 0; j < n; j++) {
			temp[j] = toMean[j][i];
		}
		toReturn[i] = mean(temp, n);
		delete temp;
	}
	return toReturn;
}
Vector mean(Matrix toMean, pos rowToMean, pos colToMean, int rowCol) {
	Vector toReturn;
	num total;
	if (rowCol == 1) {
		toReturn = new num[rowToMean];
		for (unsigned i = 0; i < rowToMean; i++) {
			total = 0.0;
			for (unsigned j = 0; j < colToMean; j++) { total += toMean[i][j]; }
			toReturn[i] = (total / num(colToMean));
		}
	}
	else if (rowCol == 2) {
		toReturn = new num[colToMean];
		for (unsigned i = 0; i < colToMean; i++) {
			total = 0.0;
			for (unsigned j = 0; j < rowToMean; j++) {
				total += toMean[j][i];
			}
			toReturn[i] = (total / num(rowToMean));
		}
	}
	else { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Select either rows(1) or columns(2)."); }
	return toReturn;
}
Vector std_dev(Matrix a) {
    pos n = a.shape[0];
    pos p = a.shape[1];
	Vector toReturn(p);
	Vector temp(n);
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < n; j++) { temp[j] = a[j][i]; }
		toReturn[i] = std_dev(temp, n);
	}
	delete temp;
	return toReturn;
}
Matrix scale(Matrix X, bool center, Vector s) {
    pos n = X.shape[0];
    pos p = X.shape[1];
    pos lenS = s.len;
	Matrix input(n, p);
	for (pos i = 0; i < n; i++) {
		//cout << i<<" "<<p << endl;
		input[i] = new num[p];
		for (int j = 0; j < p; j++) {
			input[i][j] = X[i][j];
		}
	}
	if (center) {
		Vector meanX = mean(X, n, p);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < p; j++) {
				input[i][j] -= meanX[j];

			}
		}
		delete[] meanX;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p; j++) {
			input[i][j] /= s[j];
			if (s[j] == 0) {
				input[i][j] = 0.0;
			}
		}
	}
	return input;
}
Vector multiply(Matrix a, Vector b) {
    pos rowA = a.shape[0];
    pos colA = a.shape[1];
    pos lenB = B.len;
	Vector toReturn(rowA);
	for (pos i = 0; i < rowA; i++) {
		toReturn[i] = 0.0;
		for (pos j = 0; j < colA; j++) {
			toReturn[i] += a[i][j] * b[j];
		}
	}
	return toReturn;
}
Matrix subtract(Matrix a, Matrix b) {
    pos rowA = a.shape[0];
    pos colA = a.shape[1];
	if (rowA != b.shape[0] && colA != b.shape[1]) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("SIze of A and B must be the same."); }
	Matrix toReturn(rowA, colA);
	for (int i = 0; i < rowA; i++) {
		toReturn[i] = new num[colA];
		for (int j = 0; j < colA; j++) {
			toReturn[i][j] = a[i][j] - b[i][j];
		}
	}
	return toReturn;
}
Matrix square(Matrix toSquare) {
    pos n = toSquare.shape[0];
    pos p = toSquare.shape[1];
	Matrix toReturn(n, p);
	for (int i = 0; i < n; i++) {
		toReturn[i] = new num[p];
		for (int j = 0; j < p; j++) {
			toReturn[i][j] = toSquare[i][j] * toSquare[i][j];
		}
	}
	return toReturn;
}
num arithmeticMean(Matrix toMean) {
	num mean = 0;
    pos row = toMean.shape[0];
    pos col = toMean.shape[1];
	for (pos i = 0; i < row; i++) {
		for (pos j = 0; j < col; j++) {
			mean += toMean[i][j];
		}
	}
	return mean / num(row * col);
}
Vector col_mad(Matrix toMad) {
    pos row = toMad.shape[0];
    pos col = toMad.shape[1];
	Vector toReturn(col);
	Vector temp(row);
	num mean;
	for (int i = 0; i < col; i++) {
		mean = 0;
		for (int j = 0; j < row; j++) {
			temp[j] = toMad[j][i];
		}
		toReturn[i] = mad(temp, row);
	}
	delete temp;
	return toReturn;
}
Matrix multiply(Matrix a, Matrix b) {
    pos rowA = a.shape[0];
    pos colA = a.shape[1];
    pos rowB = b.shape[0];
    pos colB = b.shape[1];
	if (colA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The column size of A and row size of B must be the same."); }
	Matrix result(rowA, colB);
	for (int i = 0; i < rowA; i++) {
		result[i] = new num[colB];
		for (unsigned j = 0; j < colB; j++) {
			result[i][j] = 0.0;
			for (unsigned k = 0; k < colA; k++) {
				result[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return result;
}

