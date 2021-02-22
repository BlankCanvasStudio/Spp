#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <tuple>
#include <math.h>
#include <chrono>
#include <stack>

using namespace std;
#define MIN_NUM -pow(1.17, 38)
#define TRUE true
#define FALSE false
#define num double
#define pos unsigned
#define Vector num*
#define Matrix num**


//Globals
auto TIME = chrono::steady_clock::now();
stack<string> currentTic;



Vector operator - (Vector a, Vector b) {
	pos lenA = a.size();
	if (lenA != b.size()) {
		cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Vectors must have the same dimensions.");
	}
	Vector result(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = a[i] - b[i];
	}
	return result;
}
Matrix operator - (Matrix a, Matrix b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos rowB = b.size();
	Matrix result(rowA, Vector(colA));
	if (rowA != rowB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Matricies must have the same dimensions to be subtracted."); }
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			result[i][j] = a[i][j] - b[i][j];
		}
	}
	return result;
}
Matrix operator - (Matrix a, Vector b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();
	if (rowA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns for the matrix and vector must be the same."); }
	Matrix result(rowA, Vector(colA));
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			result[i][j] = a[i][j] - b[i];
		}
	}
	return result;
}
Matrix subtract(Matrix a, Vector b, int rowCol) {
	Matrix toReturn;
	if (rowCol == 1) { toReturn = a - b; }
	else if (rowCol == 2) {
		pos rowA = a.size();
		pos colA = a[0].size();
		pos lenB = b.size();
		if (colA != lenB) {
			cout << "ERROR. DIMENSIONS DON'T MATCH." << endl; throw invalid_argument("Dimensions don't match.");
		}
		toReturn.resize(rowA, Vector(colA));
		for (int i = 0; i < colA; i++) {
			for (int j = 0; j < rowA; j++) {
				toReturn[j][i] = a[j][i] - b[i];
			}
		}
	}
	else {
		cout << "Enter either rows(0) or columns(1) to subtract." << endl;
		throw invalid_argument("Subtract takes either 0 or 1");
	}
	return toReturn;
}
Vector operator - (Vector a) {
	pos lenA = a.size();
	Vector result(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = -a[i];
	}
	return result;
}
Matrix operator - (Matrix a) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix result(rowA, Vector(colA));
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			result[i][j] = -a[i][j];
		}
	}
	return result;
}
template <typename T>
Vector operator - (Vector a, T b) {
	pos lenA = a.size();
	Vector result(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = a[i] - num(b);
	}
	return result;
}
template <typename T>
Matrix operator - (Matrix a, T b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	num temp = num(b);
	Matrix result(rowA, Vector(colA));
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			result[i][j] = a[i][j] - temp;
		}
	}
	return result;
}
template<typename T>
vector<int> operator > (Vector a, T b) {
	pos lenA = a.size();
	vector<int> trueLoc;
	for (unsigned i = 0; i < lenA; i++) {
		if (a[i] > num(b)) { trueLoc.push_back(i); }
	}
	return trueLoc;
}
template<typename T>
vector<int> operator < (T b, Vector a) {
	return a > b;
}
template <typename T>
vector<int> operator < (Vector a, T b) {
	pos lenA = a.size();
	vector<int> trueLoc;
	for (unsigned i = 0; i < lenA; i++) {
		if (a[i] < num(b)) { trueLoc.push_back(i); }
	}
	return trueLoc;
}
template<typename T>
vector<int> operator > (T b, Vector a) {
	return b < a;
}
template<typename T>
vector<tuple<pos, pos>> operator < (Matrix a, T b) {
	vector<tuple<pos, pos>> trueLoc;
	pos rowA = a.size();
	pos colA = a[0].size();
	num val = num(b);
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			tuple<pos, pos> temp = make_tuple(i, j);
			if (a[i][j] < val) { trueLoc.push_back(temp); }
		}
	}
	return trueLoc;
}
template<typename T>
vector<tuple<pos, pos>> operator > (Matrix a, T b) {
	pos rowA = a.size();
	pos colA = a.size();
	num val = num(b);
	vector<tuple<pos, pos>> trueLoc;
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			tuple<pos, pos> temp = make_tuple(i, j);
			if (a[i][j] > val) { trueLoc.push_back(temp); }
		}
	}
	return trueLoc;
}
Matrix operator * (Matrix a, Matrix b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos rowB = b.size();
	pos colB = b[0].size();
	Matrix result(rowA, Vector(colB));
	if (colA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The column size of A and row size of B must be the same."); }
	for (unsigned i = 0; i < colB; i++) {
		for (unsigned j = 0; j < colB; j++) {
			result[i][j] = 0.0;
			for (unsigned k = 0; k < colA; k++) {
				result[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return result;
}
Vector operator * (Matrix a, Vector b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();
	Vector result(lenB);
	for (unsigned i = 0; i < rowA; i++) {
		result[i] = 0.0;
		for (unsigned j = 0; j < colA; j++) {
			result[i] += a[i][j] * b[j];
		}
	}
	return result;
}
Matrix multiply(Matrix a, Vector b, int rowCol) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();
	Matrix toReturn(rowA, Vector(colA));
	if (rowCol == 1) {
		if (rowA != lenB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The number of rows and elements in vector must agree."); }
		for (pos i = 0; i < rowA; i++) {
			for (pos j = 0; j < colA; j++) {
				toReturn[i][j] = a[i][j] * b[j];
			}
		}
	}
	if (rowCol == 2) {
		if (colA != lenB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The number of columns and elements in vector must agree."); }
		for (pos i = 0; i < colA; i++) {
			for (pos j = 0; j < rowA; j++) {
				toReturn[j][i] = a[j][i] * b[i];
			}
		}
	}
	return toReturn;
}
template <typename T>
Vector operator * (Vector a, T b) {
	pos lenA = a.size();
	num val = num(b);
	Vector result(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = a[i] * val;
	}
	return result;
}
template<typename T>
Vector operator * (T b, Vector a) {
	return a * b;
}
Vector operator * (Vector a, Vector b) {
	pos lenA = a.size();
	pos lenB = b.size();
	if (lenA != lenB) { cout << "ERROR OCCURRED. DIMS DONT MATCH" << endl; throw invalid_argument("Dimensions don't match in vector multiplication."); }
	Vector toReturn(lenA);
	for (int i = 0; i < lenA; i++) {
		toReturn[i] = a[i] * b[i];
	}
	return toReturn;
}
Matrix operator / (Matrix a, Vector b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();

	if (colA != lenB) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("The number of columns must line up for matrix, vector division."); }
	Matrix result(rowA, Vector(colA));
	for (unsigned i = 0; i < colA; i++) {
		for (unsigned j = 0; j < rowA; j++) {
			result[j][i] = a[j][i] / b[i];
		}
	}
	return result;
}
template<typename T>
Vector operator / (Vector a, T b) {
	pos lenA = a.size();
	num val = num(b);
	Vector toReturn(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		toReturn[i] = a[i] / val;
	}
	return toReturn;
}
template<typename T>
Vector operator / (T b, Vector a) {
	pos lenA = a.size();
	num val = num(b);
	Vector toReturn(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		toReturn[i] = val / a[i];
	}
	return toReturn;
}
template<typename T>
bool operator == (Matrix a, T value) {
	num rowA = a.size();
	num colA = a[0].size();
	num val = num(value);
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			if (a[i][j] != val) {
				return false;
			}
		}
	}
	return true;
}
template<typename T>
bool operator == (Vector a, T value) {
	pos lenA = a.size();
	for (unsigned i = 0; i < lenA; i++) {
		if (a[i] != value) {
			return false;
		}
	}
	return true;
}
template<typename T, typename V>
Vector VectorFull(T len, V value) {
	Vector toReturn;
	if (len < 0.0) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Size of the vector must be positive."); }
	for (unsigned i = 0; i<unsigned(len); i++) { toReturn.push_back(num(value)); }
	return toReturn;
}
tuple<pos, pos> dim(Matrix a) {
	tuple<pos, pos> toReturn(a.size(), a[0].size());
	return toReturn;
}
unsigned dim(Vector a) {
	return a.size();
}
template<typename T, typename V, typename C>
Matrix MatrixFull(T n, V p, C value) {
	Matrix toReturn(n, Vector(p));
	for (unsigned i = 0; i<unsigned(n); i++) {
		toReturn[i] = VectorFull(p, value);
	}
	return toReturn;
}
vector<int> R_rank(Vector table) {
	pos len = table.size();
	vector<int> result(len);
	for (unsigned i = 0; i < len; i++) {
		int rank = 1;
		for (unsigned z = 0; z < len; z++) {
			if (table[z] < table[i]) { rank++; }
		}
		result[i] = rank;
	}
	return result;
}
Vector abs(Vector vec) {
	pos len = vec.size();
	Vector abs_vec(len);
	for (unsigned i = 0; i < len; i++) {
		if (vec[i] < 0.0) {
			abs_vec[i] = -1.0 * vec[i];
		}
		else {
			abs_vec[i] = vec[i];
		}
	}
	return abs_vec;
}
num max(Vector input) {
	num max = MIN_NUM;
	pos len = input.size();
	for (unsigned i = 1; i < len; i++) {
		if (max < input[i]) { max = input[i]; }
	}
	return max;

}
Vector mean(Matrix toMean, int rowCol) {
	num rowToMean = toMean.size();
	num colToMean = toMean[0].size();
	Vector toReturn;
	num total;
	if (rowCol == 1) {
		for (unsigned i = 0; i < rowToMean; i++) {
			total = 0.0;
			for (unsigned j = 0; j < colToMean; j++) { total += toMean[i][j]; }
			toReturn.push_back(total / num(colToMean));
		}
	}
	else if (rowCol == 2) {
		toReturn = VectorFull(colToMean, 0.0);
		for (unsigned i = 0; i < colToMean; i++) {
			for (unsigned j = 0; j < rowToMean; j++) { toReturn[i] += toMean[j][i]; }
		}
		for (unsigned i = 0; i < colToMean; i++) {
			toReturn[i] = toReturn[i] / num(rowToMean);
		}
	}
	else { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Select either rows(1) or columns(2)."); }
	return toReturn;
}
Matrix multiply(Matrix a, Matrix b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos rowB = b.size();
	pos colB = b[0].size();
	Matrix result(rowA, Vector(colA));
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			result[i][j] = a[i][j] * b[i][j];
		}
	}
	return result;
}
num sum(Vector a) {
	num total = 0.0;
	pos lenA = a.size();
	for (unsigned i = 0; i < lenA; i++) {
		total += a[i];
	}
	return total;
}
Vector sum(Matrix a, int rowCol) {
	Vector toReturn;
	pos rowA = a.size();
	pos colA = a.size();
	if (rowCol == 1) {
		for (unsigned i = 0; i < colA; i++) {
			toReturn.push_back(sum(a[i]));
		}
	}
	if (rowCol == 2) {
		toReturn = VectorFull(rowA, num(0.0));
		for (unsigned i = 0; i < rowA; i++) {
			for (unsigned j = 0; j < colA; j++) {
				toReturn[j] += a[i][j];
			}
		}
	}
	return toReturn;
}
num mean(Vector a) {
	return sum(a) / num(a.size());
}
num Variance(Vector a) {
	int size = a.size();

	num variance = 0;
	num t = a[0];
	num diff;
	for (int i = 1; i < size; i++) {
		t += a[i];
		diff = ((i + 1) * a[i]) - t;
		variance += (diff * diff) / ((i + 1.0) * i);
	}

	return variance / (size - 1);
}
num std_dev(Vector a) {
	return sqrt(Variance(a));
}
Vector std_dev(Matrix a, int rowCol) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Vector results;
	if (rowCol == 1) {
		for (unsigned i = 0; i < rowA; i++) {
			results.push_back(std_dev(a[i]));
		}
	}
	else if (rowCol == 2) {
		for (unsigned i = 0; i < colA; i++) {
			Vector temp(rowA);
			for (unsigned j = 0; j < rowA; j++) {
				temp[j] = a[j][i];
			}
			results.push_back(std_dev(temp));
		}
	}
	else {
		cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Please enter either rows(1) or columns(2)");
	}
	return results;
}
Matrix multiply(Matrix a, Vector b) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix result(rowA, Vector(colA));
	for (unsigned i = 0; i < a.size(); i++) {
		for (unsigned j = 0; j < colA; j++) {
			result[i][j] = a[i][j] * b[i];
		}
	}
	return result;
}
Matrix scale(Matrix X, bool center, Vector s) {
	Matrix toReturn;
	if (center) { toReturn = subtract(X, mean(X, 2), 2); }
	else { toReturn = X; }
	return toReturn / s;
}
num max(Matrix a) {
	num max = MIN_NUM;
	pos rowA = a.size();
	pos colA = a[0].size();
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			if (max < a[i][j]) { max = a[i][j]; }
		}
	}
	return max;
}
Matrix abs(Matrix a) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix result(rowA, Vector(colA));
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) { result[i][j] = abs(a[i][j]); }
	}
	return result;
}
template <typename T>
Vector pow(Vector a, T degree) {
	pos lenA = a.size();
	Vector result(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		result[i] = pow(a[i], degree);
	}
	return result;
}
template <typename T>
Matrix pow(Matrix a, T degree) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix toReturn(rowA, Vector(colA));
	for (unsigned i = 0; i < rowA; i++) {
		toReturn[i] = pow(a[i], degree);
	}
	return toReturn;
}
Matrix log(Matrix a) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix toReturn(rowA, Vector(colA));
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < colA; j++) {
			toReturn[i][j] = log(a[i][j]);
		}
	}
	return toReturn;
}
Vector log(Vector a) {
	pos lenA = a.size();
	Vector toReturn(lenA);
	for (unsigned i = 0; i < lenA; i++) {
		toReturn[i] = log(a[i]);
	}
	return toReturn;
}
Matrix readMatrix_csv(string filename) {
	Matrix toReturn;
	Vector toReturnV;
	fstream file;
	file.open(filename.c_str(), ios::in);
	string newLine;
	unsigned lastcomma;
	string subs;
	pos len;
	if (!file.is_open()) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Invalid file name."); }
	while (!file.eof()) {
		getline(file, newLine);
		lastcomma = 0;
		toReturnV.clear();
		len = newLine.length();
		for (unsigned i = 0; i < len; i++) {
			if (newLine[i] == ',' || i == len - 1) {
				subs = newLine.substr(lastcomma, i);
				//if(subs !="\n" && subs!="\t"){
				try {
					toReturnV.push_back(stod(subs.c_str()));
				}
				catch (int e) {
					cout << "Invalid number read in. Skipping." << endl;
				}
				lastcomma = i + 1;
				//}
			}
		}
		toReturn.push_back(toReturnV);
	}
	return toReturn;

}
Vector readVector_csv(string filename) {
	Vector toReturn;
	fstream file;
	file.open(filename.c_str(), ios::in);
	string newLine;
	if (!file.is_open()) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Invalid file name."); }
	getline(file, newLine);
	unsigned lastcomma = 0;
	pos len = newLine.length();
	for (unsigned i = 0; i < len; i++) {
		if (newLine[i] == ',' || i == len - 1) {
			string subs = newLine.substr(lastcomma, i);
			if (atof("Nan") != atof(subs.c_str())) {
				toReturn.push_back(atof(subs.c_str()));
				lastcomma = i + 1;
			}

		}
	}
	return toReturn;
}
template <typename T>
Matrix diag(Matrix a, T fill) {
	pos rowA = a.size();
	pos colA = a[0].size();
	num val = num(fill);
	Matrix toReturn(rowA, Vector(colA));
	if (rowA != colA) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Diag can only act on a square matrix."); }
	for (unsigned i = 0; i < rowA; i++) {
		for (unsigned j = 0; j < rowA; j++) {
			toReturn[i][j] = a[i][j];
			if (i == j) { toReturn[i][j] = val; }
		}
	}
	return toReturn;
}
void print(Vector a) {
	pos lenA = a.size();
	for (pos i = 0; i < lenA; i++) {
		cout << a[i] << ", ";
	}
	cout << endl;
}
void print(Matrix a) {
	pos rowA = a.size();
	pos colA = a[0].size();
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			cout << a[i][j] << ", ";
		}
		cout << endl;
	}
	cout << endl;
}
void print(tuple<pos, pos> a) {
	cout << '(' << get<0>(a) << ',' << get<1>(a) << ')' << endl;
}
template <typename T>
void print(T a) {
	cout << a << endl;
}
void tic(string ticName) {
	TIME = chrono::steady_clock::now();
	currentTic.push(ticName);
}

void toc() {
	num diff = chrono::duration<num>(chrono::steady_clock::now() - TIME).count();
	cout << "Time to execute " << currentTic.top() << ": " << diff << endl;
	currentTic.pop();
}
pos len(Vector a) {
	return a.size();
}
template <typename T>
Vector col(Matrix a, T col) {
	pos column = pos(col);
	pos rowA = a.size();
	pos colA = a[0].size();
	Vector toReturn(rowA);
	for (int i = 0; i < rowA; i++) {
		toReturn[i] = a[i][column];
	}
	return toReturn;
}
template <typename T>
Matrix adjust(Matrix a, T val) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix toReturn(rowA, Vector(colA));
	for (int i = 0; i < rowA; i++) {
		for (int j = 0; j < colA; j++) {
			toReturn[i][j] = a[i][j] - copysign(num(val), a[i][j]);
		}
	}
	return toReturn;
}