//R.h
#include <fstream>
#include <string>
#include <iostream>
#include <tuple>
#include <math.h>
#include <chrono>
#include <stack>
#include <unordered_map>
#define MIN_NUM -pow(1.17, 38)
#define TRUE true
#define FALSE false
#define FALSE false
#define num double
#define pos unsigned
#define Matrix vector<vector<num> >
#define Vector vector<num>
#define Matrix2 num**
#define Vector2 num*
using namespace std;


//DataTypes
struct SVector {
	num* base;
	pos len;
};
struct SMatrix {
	num** base;
	pos row;
	pos col;
};
struct cv_return {
	Vector lambda;
	Vector vr;
	num lambda_min;
	num r_min;
	Matrix optloc;
	Matrix cvm;

	cv_return operator()(Vector lambda, Vector vr, num lambda_min, num r_min, Matrix optloc, Matrix cvm) {
		cv_return a;
		a.lambda = lambda;
		a.vr = vr;
		a.lambda_min = lambda_min;
		a.r_min = r_min;
		a.optloc = optloc;
		a.cvm = cvm;
		return a;
	}
	void operator =(cv_return a) {
		lambda = a.lambda;
		vr = a.vr;
		lambda_min = a.lambda_min;
		r_min = a.r_min;
		optloc = a.optloc;
		cvm = a.cvm;
	}
};

struct returnList {
	Vector B0;
	Matrix B;
	Matrix Bini;
	Vector glam;
	Vector weight;
	Vector Wini;
	num gk = 0.0;
	returnList operator ()(Matrix inB, Vector inB0, Vector inweight, Matrix inBini, Vector inWini) {
		returnList a;
		a.B0 = inB0;
		a.B = inB;
		a.Bini = inBini;
		a.weight = inweight;
		a.Wini = inWini;
		return a;
	}
	returnList operator ()(Vector inB0, Matrix inB, Vector inglam, num ingk) {
		returnList a;
		a.B0 = inB0;
		a.B = inB;
		a.glam = inglam;
		a.gk = ingk;
		return a;
	}
	void operator = (returnList toEqual) {
		B0 = toEqual.B0;
		B = toEqual.B;
		Bini = toEqual.Bini;
		glam = toEqual.glam;
		weight = toEqual.weight;
		Wini = toEqual.Wini;
		gk = toEqual.gk;
	}


};

struct ptr_returnList {
	Vector2 B0;
	Matrix2 B;
	Matrix2 Bini;
	Vector2 glam;
	Vector2 weight;
	Vector2 Wini;
	num gk = 0.0;
	ptr_returnList operator ()(Matrix2 inB, Vector2 inB0, Vector2 inweight, Matrix2 inBini, Vector2 inWini) {
		ptr_returnList a;
		a.B0 = inB0;
		a.B = inB;
		a.Bini = inBini;
		a.weight = inweight;
		a.Wini = inWini;
		return a;
	}
	ptr_returnList operator ()(Vector2 inB0, Matrix2 inB, Vector2 inglam, num ingk) {
		ptr_returnList a;
		a.B0 = inB0;
		a.B = inB;
		a.glam = inglam;
		a.gk = ingk;
		return a;
	}
	void operator = (ptr_returnList toEqual) {
		B0 = toEqual.B0;
		B = toEqual.B;
		Bini = toEqual.Bini;
		glam = toEqual.glam;
		weight = toEqual.weight;
		Wini = toEqual.Wini;
		gk = toEqual.gk;
	}
};

//Globals
auto TIME = chrono::steady_clock::now();
stack<string> currentTic;

void tic(string ticName) {
	TIME = chrono::steady_clock::now();
	currentTic.push(ticName);
}

void toc() {
	num diff = chrono::duration<num>(chrono::steady_clock::now() - TIME).count();
	cout << "Time to execute " << currentTic.top() << ": " << diff << endl;
	currentTic.pop();
}

/*
//PROTOTYPES
template <typename T>
Vector pow(Vector &a, T degree);
Matrix scale(Matrix &X, bool center, Vector &s);
num max(Vector &input);
num max(Matrix &a);
Matrix abs(Matrix &a);
Vector sum(Matrix &a);
num sum(Vector &a);
Matrix multiply(Matrix &a, Matrix &b);
Matrix multiply(Matrix &a, Vector &b);
vector<int> R_rank(Vector &table);
Vector abs(Vector &vec);
num std_dev(Vector &a);
Vector std_dev(Matrix &a, int rowCol);
num mean(Vector &a);
Vector mean(Matrix &a, int rowCol);
Vector sum(Matrix &a);
num sum(Vector &a);
Vector mean(Matrix &toMean, int rowCol);
void print(Matrix &a);
void print(Vector &a);
template <typename T, typename V>
Vector VectorFull(T len, V value);
Matrix readMatrix_csv(string filename);
Vector readVector_csv(string filename);
tuple<pos, pos> dim(Matrix &a);
unsigned dim(Vector &a);
template <typename T, typename C>
Vector rep(T val, C len);
vector<int> R_rank(Vector &table, string ties);
vector<int> R_rank(Vector &table);

//ACTUAL CODE
Vector operator - (Vector &a, Vector &b){
	pos lenA = a.size();
	if(lenA != b.size()){
		cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Vectors must have the same dimensions.");
	}
	Vector result (lenA);
	for(unsigned i=0;i<lenA;i++){
		result[i] = a[i] - b[i];
	}
	return result;
}
template <typename T>
bool operator != (Vector &input, T val) {
	pos size = input.size();
	for (int i = 0; i < size; i++) {
		if (input[i] != num(val)) { return true; }
	}
	return false;
}
template <typename T>
bool operator == (Vector &input, T val) {
	pos size = input.size();
	for (int i = 0; i < size; i++) {
		if (input[i] != num(val)) { return false; }
	}
	return true;
}
Matrix operator - (Matrix &a, Matrix &b){
	pos rowA = a.size();
	pos colA = a[0].size();
	pos rowB = b.size();
	Matrix result (rowA, Vector(colA));
	if(rowA!=rowB){cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Matricies must have the same dimensions to be subtracted.");}
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			result[i][j] = a[i][j] - b[i][j];
		}
	}
	return result;
}
Matrix operator - (Matrix &a, Vector &b){
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();
	if(rowA != lenB){cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("The number of columns for the matrix and vector must be the same.");}
	Matrix result (rowA, Vector(colA));
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			result[i][j] = a[i][j] - b[i];
		}
	}
	return result;
}
Matrix subtract(Matrix &a, Vector &b, int rowCol){
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
		throw invalid_argument("Subtract takes either 0 or 1");}
	return toReturn;
}
Vector operator - (Vector &a){
	pos lenA = a.size();
	Vector result(lenA);
	for(unsigned i=0;i<lenA;i++){
		result[i] = -a[i];
	}
	return result;
}
Matrix operator - (Matrix &a){
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix result(rowA, Vector(colA));
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			result[i][j] = -a[i][j];
		}
	}
	return result;
}
template <typename T>
Vector operator - (Vector &a, T b){
	pos lenA = a.size();
	Vector result (lenA);
	for(unsigned i=0;i<lenA;i++){
		result[i] = a[i] - num(b);
	}
	return result;
}
//template <typename T>
Matrix operator - (Matrix &a, num b){
	pos rowA = a.size();
	pos colA = a[0].size();
	num temp = num(b);
	Matrix result (rowA, Vector(colA));
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			result[i][j] = a[i][j] - temp;
		}
	}
	return result;
}
template<typename T>
vector<int> operator > (Vector a, T b){
	pos lenA = a.size();
	vector<int> trueLoc;
	for(unsigned i=0;i<lenA;i++){
		if(a[i]>num(b)){trueLoc.push_back(i);}
	}
	return trueLoc;
}
template<typename T>
vector<int> operator < (T b, Vector &a){
	return a > b;
}
template <typename T>
vector<int> operator < (Vector &a, T b){
	pos lenA = a.size();
	vector<int> trueLoc;
	for(unsigned i=0;i<lenA;i++){
		if(a[i]<num(b)){trueLoc.push_back(i);}
	}
	return trueLoc;
}
template<typename T>
vector<int> operator > (T b, Vector &a){
	return b < a;
}
template<typename T>
vector<tuple<pos, pos>> operator < (Matrix &a, T b){
	vector<tuple<pos, pos>> trueLoc;
	pos rowA = a.size();
	pos colA = a[0].size();
	num val = num(b);
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			tuple<pos, pos> temp = make_tuple(i,j);
			if(a[i][j]<val){trueLoc.push_back(temp);}
		}
	}
	return trueLoc;
}
template<typename T>
vector<tuple<pos, pos>> operator > (Matrix &a, T b){
	pos rowA = a.size();
	pos colA = a.size();
	num val = num(b);
	vector<tuple<pos, pos>> trueLoc;
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			tuple<pos, pos> temp = make_tuple(i, j);
			if(a[i][j]>val){trueLoc.push_back(temp);}
		}
	}
	return trueLoc;
}
Matrix operator * (Matrix &a, Matrix &b){
	pos rowA = a.size();
	pos colA = a[0].size();
	pos rowB = b.size();
	pos colB = b[0].size();
	Matrix result (rowA, Vector(colB));
	if (colA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The column size of A and row size of B must be the same."); }
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0; j<colB;j++){
			result[i][j] = 0.0;
			for(unsigned k=0;k<colA;k++){
				result[i][j]+= a[i][k]*b[k][j];
			}
		}
	}
	return result;
}
Vector operator * (Matrix &a, Vector &b){
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();
	Vector result (lenB);
	for(unsigned i=0; i<rowA;i++){
		result[i] = 0.0;
		for(unsigned j=0;j<colA;j++){
			result[i] += a[i][j] * b[j];
		}
	}
	return result;
}
Matrix multiply(Matrix a, Vector b, int rowCol) {
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();
	Matrix toReturn (rowA, Vector(colA));
	if (rowCol == 1) {
		if (rowA != lenB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The number of rows and elements in vector must agree.");}
		for (pos i = 0; i < rowA; i++) {
			for (pos j = 0; j<colA ; j++) {
				toReturn[i][j] = a[i][j] * b[j];
			}
		}
	}
	if (rowCol == 2) {
		if (colA != lenB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The number of columns and elements in vector must agree.");}
		for (pos i = 0; i < colA; i++) {
			for(pos j=0; j < rowA; j++){
				toReturn[j][i] = a[j][i]*b[i];
			}
		}
	}
	return toReturn;
}
template <typename T>
Vector operator * (Vector &a, num &b){
	pos lenA = a.size();
	num val = num(b);
	Vector result (lenA);
	for(unsigned i=0;i<lenA;i++){
		result[i] = a[i] * val;
	}
	return result;
}
template<typename T>
Vector operator * (num &b, Vector &a){
	return a * b;
}
Vector operator * (Vector &a, Vector &b) {
	pos lenA = a.size();
	pos lenB = b.size();
	if (lenA != lenB) { cout << "ERROR OCCURRED. DIMS DONT MATCH" << endl; throw invalid_argument("Dimensions don't match in vector multiplication."); }
	Vector toReturn(lenA);
	for (int i = 0; i < lenA;i++) {
		toReturn[i] = a[i] * b[i];
	}
	return toReturn;
}
Matrix operator / (Matrix &a, Vector &b){
	pos rowA = a.size();
	pos colA = a[0].size();
	pos lenB = b.size();
	if(colA != lenB){cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("The number of columns must line up for matrix, vector division.");}
	Matrix result(rowA, Vector(colA));
	for(unsigned i=0;i<colA;i++){
		for(unsigned j=0;j<rowA;j++){
			result[j][i] = a[j][i] / b[i];
		}
	}
	return result;
}
template<typename T>
Vector operator / (Vector &a, T b){
	pos lenA = a.size();
	num val = num(b);
	Vector toReturn (lenA);
	for(unsigned i=0;i<lenA;i++){
		toReturn[i] = a[i]/val;
	}
	return toReturn;
}
template<typename T>
Vector operator / (T b, Vector &a){
	pos lenA = a.size();
	num val = num(b);
	Vector toReturn (lenA);
	for(unsigned i=0;i<lenA;i++){
		toReturn[i] = val / a[i];
	}
	return toReturn;
}
template<typename T>
bool operator == (Matrix &a, T value){
	num rowA = a.size();
	num colA = a[0].size();
	num val = num(value);
	for(unsigned i=0; i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			if(a[i][j] != val){
				return false;
			}
		}
	}
	return true;
}

template<typename T, typename V>
Vector VectorFull(T len, V value){
	Vector toReturn;
	if(len < 0.0){cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Size of the vector must be positive.");}
	for(unsigned i=0;i<unsigned(len);i++){toReturn.push_back(num(value));}
	return toReturn;
}
tuple<pos, pos> dim(Matrix &a){
	tuple<pos, pos> toReturn(a.size(), a[0].size());
	return toReturn;
}
unsigned dim(Vector &a){
	return a.size();
}
template<typename T, typename V, typename C>
Matrix MatrixFull(T n, V p, C value){
	Matrix toReturn (n, Vector(p));
	for(unsigned i=0;i<unsigned(n);i++){
		toReturn[i] = VectorFull(p, value);
	}
	return toReturn;
}
vector<int> R_rank(Vector &table) {
	return R_rank(table, "None");
}
vector<int> R_rank(Vector &table, string ties){
	pos len = table.size();
	vector<int> result (len);
	for(unsigned i=0;i<len;i++){
		int rank=1;
		for(unsigned z=0;z<len;z++){
			if(table[z]<table[i]){rank++;}
		}
	result[i]=rank;
	}
	if (ties == "average") {
		Vector unchecked = rep(1, len);
		Vector same = rep(0, len);
		//pos index = 0;
		pos numSame;
		pos totalIndex;
		num avgIndex;
		for (pos index = 0; index < len; index++) {
			while (index < len && unchecked[index] == 0.0) { index++; }
			if (index >= len) { break; }
			numSame = 0;
			totalIndex = index;
			unchecked[index] = 0;
			for (int i = index + 1; i < len; i++) {
				if (result[index] == result[i]) {
					same[index] = 1;
					same[i] = 1;
					numSame++;
					totalIndex += i;
					unchecked[i] = 0;
				}
			}
			if (same[index] == 1) {
				avgIndex = num(totalIndex) / num(numSame);
				for (int i = 0; i < len; i++) {
					if (same[i] == 1) {
						result[i] = avgIndex;
						same[i] = 0;
					}
				}
			}
		}
	}
	return result;
}
Vector abs(Vector &vec){
	pos len = vec.size();
	Vector abs_vec (len);
	for(unsigned i=0; i<len;i++){
		if(vec[i] < 0.0){
			abs_vec[i] = -1.0*vec[i];
		}
		else{
			abs_vec[i] = vec[i];
		}
	}
	return abs_vec;
}
num max(Vector &input){
	num max = MIN_NUM;
	pos len = input.size();
	for(unsigned i=1;i<len;i++){
		if(max<input[i]){max=input[i];}
	}
	return max;

}

Vector mean(Matrix &toMean, int rowCol){
	num rowToMean = toMean.size();
	num colToMean = toMean[0].size();
	Vector toReturn;
	num total;
	if(rowCol==1){
		for(unsigned i=0; i<rowToMean;i++){
			total = 0.0;
			for(unsigned j=0; j<colToMean;j++){total += toMean[i][j];}
			toReturn.push_back(total/num(colToMean));
		}
	}
	else if(rowCol==2){
		for(unsigned i=0; i< colToMean;i++){
			total = 0.0;
			for(unsigned j=0; j< rowToMean;j++){
				total += toMean[j][i];
			}
			toReturn.push_back(total / num(rowToMean));
		}
	}
	else{cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Select either rows(1) or columns(2).");}
	return toReturn;
}
num arithmeticMean(Matrix &toMean) {
	num mean = 0;
	num rowToMean = toMean.size();
	num colToMean = toMean[0].size();
	for (pos i = 0; i < rowToMean; i++) {
		for (pos j = 0; j < colToMean; j++) {
			mean += toMean[i][j];
		}
	}
	return mean / num(rowToMean * colToMean);
}
Matrix multiply(Matrix &a, Matrix &b){
	pos rowA = a.size();
	pos colA = a[0].size();
	pos rowB = b.size();
	pos colB = b[0].size();
	Matrix result(rowA, Vector(colA));
	for(unsigned i=0; i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			result[i][j] = a[i][j] * b[i][j];
		}
	}
	return result;
}
num sum(Vector &a){
	num total=0.0;
	pos lenA = a.size();
	for(unsigned i=0;i<lenA;i++){
		total += a[i];
	}
	return total;
}
Vector sum(Matrix &a, int rowCol){
	Vector toReturn;
	pos rowA = a.size();
	pos colA = a.size();
	if(rowCol==1){
		for(unsigned i=0; i<colA;i++){
			toReturn.push_back(sum(a[i]));
		}
	}
	if(rowCol==2){
		toReturn = VectorFull(rowA, num(0.0));
		for(unsigned i=0;i<rowA;i++){
			for(unsigned j=0;j<colA;j++){
				toReturn[j] += a[i][j];
			}
		}
	}
	return toReturn;
}
num mean(Vector &a){
	return sum(a)/num(a.size());
}
num Variance(Vector &a) {
	uint64_t size = a.size();

	num variance = 0;
	num diff;
	num meanA = mean(a);
	num total = 0.0;
	for (int i =0; i < size; i++) {
		total += pow((a[i] - meanA), 2);
	}
	return total / num(size-1);
}
num std_dev(Vector &a) {
	return sqrt(Variance(a));
}
//num std_dev(Vector a){
//	return sqrt(sum(pow((a - mean(a)), 2.0))/ num(a.size()));
//}
Vector std_dev(Matrix &a, int rowCol){
	pos rowA = a.size();
	pos colA = a[0].size();
		//Flipped these to see if it worked
	Vector results;
	if(rowCol==1){
		for(unsigned i=0;i<rowA;i++){
			results.push_back(std_dev(a[i]));
		}
	}
	else if(rowCol==2){
		for(unsigned i=0;i<colA;i++){
			Vector temp (rowA);
			for(unsigned j=0;j<rowA;j++){
				temp[j] = a[j][i];
			}
			results.push_back(std_dev(temp));
		}
	}
	else{
		cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Please enter either rows(1) or columns(2)");
	}
	return results;
}
Matrix multiply(Matrix &a, Vector &b){
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix result(rowA, Vector(colA));
	for(unsigned i=0;i<a.size();i++){
		for(unsigned j=0; j<colA;j++){
			result[i][j] = a[i][j] * b[i];
		}
	}
	return result;
}
Matrix scale(Matrix &X, bool center, Vector &s){
	Matrix toReturn;
	if(center){toReturn = subtract(X, mean(X, 2), 2);}
	else { toReturn = X; }
	return toReturn / s;
}
num max(Matrix a){
	num max = MIN_NUM;
	pos rowA = a.size();
	pos colA = a[0].size();
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			if(max<a[i][j]){max=a[i][j];}
		}
	}
	return max;
}
Matrix abs(Matrix &a){
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix result(rowA, Vector(colA));
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){result[i][j] = abs(a[i][j]);}
	}
	return result;
}
template <typename T>
Vector pow(Vector &a, T degree){
	pos lenA = a.size();
	Vector result (lenA);
	for(unsigned i=0;i<lenA;i++){
		result[i] = pow(a[i], degree);
	}
	return result;
}
template <typename T>
Matrix pow(Matrix &a, T degree){
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix toReturn (rowA, Vector(colA));
	for(unsigned i=0; i<rowA;i++){
		toReturn[i] = pow(a[i], degree);
	}
	return toReturn;
}
Matrix log(Matrix &a){
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix toReturn(rowA, Vector(colA));
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<colA;j++){
			toReturn[i][j] = log(a[i][j]);
		}
	}
	return toReturn;
}
Vector log(Vector &a){
	pos lenA = a.size();
	Vector toReturn (lenA);
	for(unsigned i=0; i<lenA;i++){
		toReturn[i] = log(a[i]);
	}
	return toReturn;
}
Matrix readMatrix_csv(string filename, bool header){
	Matrix toReturn;
	Vector toReturnV;
	fstream file;
	file.open(filename.c_str(), ios::in);
	string newLine;
	unsigned lastcomma;
	string subs;
	pos len;
	if(!file.is_open()){cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Invalid file name.");}
	if (header) {
		getline(file, newLine);
	}
	while(!file.eof()){
		getline(file, newLine);
		lastcomma=0;
		toReturnV.clear();
		len = newLine.length();
		for(unsigned i=0;i<len;i++){
			if(newLine[i] == ',' || i==len-1){
				subs = newLine.substr(lastcomma, i);
				//if(subs !="\n" && subs!="\t"){
					try {
						toReturnV.push_back(stod(subs.c_str()));
					}
					catch(int e){
						cout << "Invalid number read in. Skipping." << endl;
					}
					lastcomma = i + 1;
				//}
			}
		}
		if (toReturnV.size() != 0) {
			toReturn.push_back(toReturnV);
		}
	}
	return toReturn;
}
Vector readVector_csv(string filename){
	Vector toReturn;
	fstream file;
	file.open(filename.c_str(), ios::in);
	string newLine;
	if(!file.is_open()){cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Invalid file name.");}
	getline(file, newLine);
	unsigned lastcomma=0;
	pos len = newLine.length();
	for(unsigned i=0;i<len;i++){
		if(newLine[i] == ',' || i==len-1){
			string subs = newLine.substr(lastcomma, i);
			if(atof("Nan") != atof(subs.c_str())){
				toReturn.push_back(atof(subs.c_str()));
				lastcomma = i+1;
			}

		}
	}
	return toReturn;
}
template <typename T>
Matrix diag(Matrix &a, T fill){
	pos rowA = a.size();
	pos colA = a[0].size();
	num val = num(fill);
	Matrix toReturn (rowA, Vector(colA));
	if(rowA != colA){cout<<"ERROR ENCOUNTERED"<<endl; throw invalid_argument("Diag can only act on a square matrix.");}
	for(unsigned i=0;i<rowA;i++){
		for(unsigned j=0;j<rowA;j++){
			toReturn[i][j]= a[i][j];
			if(i==j){toReturn[i][j] = val;}
		}
	}
	return toReturn;
}
void print(Vector &a){
	pos lenA = a.size();
	for(pos i=0;i<lenA;i++){
		cout<<a[i]<<", ";
	}
	cout<<endl;
}
void print(Matrix &a){
	pos rowA = a.size();
	pos colA = a[0].size();
	for(pos i=0;i<rowA;i++){
		for(pos j=0;j<colA;j++){
			cout<<a[i][j]<<", ";
		}
		cout<<endl;
	}
	cout<<endl;
}
void print(tuple<pos, pos> &a){
	cout<<'('<<get<0>(a)<<','<<get<1>(a)<<')'<<endl;
}
template <typename T>
void print(T a){
	cout<<a<<endl;
}

pos len(Vector &a) {
	return a.size();
}
template <typename T>
Vector col(Matrix &a, T col) {
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
Matrix adjust(Matrix &a, T val) {
	pos rowA = a.size();
	pos colA = a[0].size();
	Matrix toReturn(rowA, Vector(colA));
	for (pos i = 0; i < rowA; i++) {
		for (pos j = 0; j < colA; j++) {
			toReturn[i][j] = a[i][j] - copysign(num(val), a[i][j]);
		}
	}
	return toReturn;
}
void mean(Matrix2 &toMean, pos n, pos p, int rowCol, Vector2 &result) {
	num total;
	if (rowCol == 1) {
		for (unsigned i = 0; i < n; i++) {
			total = 0.0;
			for (unsigned j = 0; j < p; j++) { total += toMean[i][j]; }
			result[i] = (total / num(p));
		}
	}
	else if (rowCol == 2) {
		for (pos i = 0; i < p; i++) { result[i] = 0.0; }
		for (unsigned i = 0; i < p; i++) {
			for (unsigned j = 0; j < n; j++) { result[i] += toMean[j][i]; }
		}
		for (unsigned i = 0; i < p; i++) {
			result[i] = result[i] / num(n);
		}
	}
	else { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Select either rows(1) or columns(2)."); }
}

void matMul(Matrix2 &a, pos rowA, pos colA, Matrix2 &b, pos rowB, pos colB, Matrix2 &result) {
	//for (int i = 0; i < rowA; i++) { result[i] = new num[colB]; }
	if (colA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The column size of A and row size of B must be the same."); }
	int i = 0;
	while (i < rowA) {
		for (int j = 0; j < colB; j++) {
			result[i][j] = 0.0;
			for (int k = 0; k < colA; k++) {
				result[i][j] += a[i][k] * b[k][j];
			}
		}
		i++;
	}
}
num sum(Vector2 &toSum, pos len) {
	num total = 0.0;
	for (int i = 0; i < len; i++) {
		total += toSum[i];
	}
	return total;
}
void abs(Vector2 &toAbs, pos len) {
	for (int i = 0; i < len; i++) {if (toAbs[i] < 0.0) { toAbs[i] = -toAbs[i]; }}
}
void abs(Matrix2 &toAbs, pos row, pos col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (toAbs[i][j] < 0.0) { toAbs[i][j] = -toAbs[i][j]; }
		}
	}
}
num R_rank_top(Vector2 &table, pos len) {
	vector<int> result(len);
	for (unsigned i = 0; i < len; i++) {
		int rank = 1;
		for (unsigned z = 0; z < len; z++) {
			if (table[z] < table[i]) { rank++; }
		}
		result[i] = rank;
	}
	return result[0];
}
template <typename T, typename C>
Vector rep(T val, C len) {
	Vector toReturn(len);
	for (int i = 0; i < len; i++) { toReturn[i] = val; }
	return toReturn;
}

void saveAsCSV(Matrix &Input, string filename) {
	if (filename.substr(filename.length()- 5) != ".csv") {
		filename.append(".csv");
	}
	string line;
	pos rowInput = Input.size();
	pos columnInput = Input[0].size();
	ofstream outputFile;
	outputFile.open(filename);
	for (int i = 0; i < rowInput; i++) {
		line = "";
		for (int j = 0; j < columnInput; j++) {
			line.append(to_string(Input[i][j]));
			line.append(",");
		}
		line.append("\n");
		outputFile << line;
	}
	outputFile.close();
}
num mad(Vector &toMad) {
	num absSum = 0;
	pos len = toMad.size();
	num meanVal = mean(toMad);
	for (int i = 0; i < len; i++) {
		absSum += abs(toMad[i] - meanVal);
	}
	return absSum / len;
}
Vector col_mad(Matrix &toMad) {
	Vector toReturn;
	Vector temp;
	pos row = toMad.size();
	pos col = toMad[0].size();
	num mean;
	for (int i = 0; i < col; i++) {
		temp.clear();
		temp.clear();
		mean = 0;
		for (int j = 0; j < row; j++) {
			temp.push_back(toMad[j][i]);
		}
		toReturn.push_back(mad(temp));
	}
	return toReturn;
}
Vector seq(num lower, num upper, num steps) {
	Vector toReturn(steps);
	for (int i = 0; i < steps; i++) {
		toReturn[i] = lower + (i * ((upper - lower) / steps));
	}
	return toReturn;
}

Matrix t(Matrix &toTranspose) {
	pos n = toTranspose.size();
	pos q = toTranspose[0].size();
	Matrix toReturn(q, Vector(n));
	for (pos i = 0; i < n; i++) {
		for (pos j = 0; j < q; j++) {
			toReturn[j][i] = toTranspose[i][j];
		}
	}
	return toReturn;
}

num min(Matrix &a) {
	pos row = a.size();
	pos col = a.size();
	num min_num = DBL_MAX;
	for (pos i = 0; i < row; i++) {
		for (pos j = 0; j < col; j++) {
			if (a[i][j] < min_num) {
				min_num = a[i][j];
			}
		}
	}
	return min_num;
}

Vector min_index(Matrix &a) {
	pos row = a.size();
	pos col = a.size();
	num min_num = DBL_MAX;
	Vector min_location(2);
	for (pos i = 0; i < row; i++) {
		for (pos j = 0; j < col; j++) {
			if (a[i][j] < min_num) {
				min_location[0] = i;
				min_location[1] = j;
				min_num = a[i][j];
			}
		}
	}
	return min_location;
}

Matrix mean_over_3D(vector<Matrix> &toMean) {
	Matrix toReturn = MatrixFull(toMean[0].size(), toMean[0][0].size(), 0);
	pos TD = toMean.size();
	pos row = toMean[0].size();
	pos col = toMean[0][0].size();

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			for (int k = 0; k < TD; k++) {
				toReturn[i][j] += toMean[k][i][j];
			}
			toReturn[i][j] /= TD;
		}
	}
	return toReturn;
}
void save_CV_return(cv_return &out, string filename) {
	if (filename.substr(filename.length() - 5) != ".csv") {
		filename.append(".csv");
	}
	string line;
	ofstream outputFile;
	outputFile.open(filename);
	string toAdd = "";
	for (int i = 0; i < out.lambda.size(); i++) {
		toAdd += " " + to_string(out.lambda[i]);
	}
	outputFile << "lambda:" + toAdd+"\n";

	toAdd = "";
	for (int i = 0; i < out.vr.size(); i++) {
		toAdd += " " + to_string(out.vr[i]);
	}
	outputFile << "vr:" + toAdd + "\n";


	outputFile << "r_min: " + to_string(out.r_min)+"\n";
	outputFile << "lambda_min: " + to_string(out.lambda_min) + "\n";

	outputFile << "optloc:\n";
	for (int i = 0; i < out.optloc.size(); i++) {
		toAdd = "";
		for (int j = 0; j < out.optloc[0].size(); j++) {
			toAdd.append(to_string(out.optloc[i][j]));
			toAdd.append(",");
		}
		toAdd.append("\n");
		outputFile << toAdd;
	}
	outputFile << "cvm:\n";
	for (int i = 0; i < out.cvm.size(); i++) {
		toAdd = "";
		for (int j = 0; j < out.cvm[0].size(); j++) {
			toAdd.append(to_string(out.cvm[i][j]));
			toAdd.append(",");
		}
		toAdd.append("\n");
		outputFile << toAdd;
	}
	outputFile.close();
}
*/

//Official Vector2 functions
void print(Matrix2 toPrint, pos n, pos p) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p; j++) {
			cout << toPrint[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
void print(Vector2 toPrint, pos n) {
	for (int i = 0; i < n; i++) {
		cout << toPrint[i] << endl;
	}
}

num sum(Vector2 a, pos size) {
	num toReturn = 0;
	for (pos i = 0; i < size; i++) {
		toReturn += a[i];
	}
	return toReturn;
}
num mean(Vector2 a, pos size) {
	return sum(a, size) / num(size);
}
Vector2 mean(Matrix2 toMean, pos n, pos p) {
	Vector2 toReturn = new num[p];
	for (int i = 0; i < p; i++) {
		Vector2 temp = new num[n];
		for (int j = 0; j < n; j++) {
			temp[j] = toMean[j][i];
		}
		toReturn[i] = mean(temp, n);
		delete[] temp;
	}
	return toReturn;
}
num Variance(Vector2 a, pos size) {
	num meanA = mean(a, size);
	num total = 0.0;
	for (int i = 0; i < size; i++) {
		total += pow((a[i] - meanA), 2);
	}
	if (size - 1 == 0) { return 0; }
	return total / num(size - 1);
}
num std_dev(Vector2 a, pos size) {
	return sqrt(Variance(a, size));
}

Vector2 std_dev(Matrix2 a, pos n, pos p) {
	Vector2 toReturn = new num[p];
	Vector2 temp = new num[n];
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < n; j++) { temp[j] = a[j][i]; }
		toReturn[i] = std_dev(temp, n);
	}
	delete[] temp;
	return toReturn;
}
Matrix2 scale(Matrix2 X, pos n, pos p, bool center, Vector2 s, pos lenS) {
	Matrix2 input = new num * [n];
	for (pos i = 0; i < n; i++) {
		input[i] = new num[p];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p; j++) {
			input[i][j] = X[i][j];
		}
	}
	if (center) {
		Vector2 meanX = mean(X, n, p);
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
Vector2 mean(Matrix2 toMean, pos rowToMean, pos colToMean, int rowCol) {
	Vector2 toReturn;
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
void swap(Vector2 a, Vector2 b) {
	num t = *a;
	*a = *b;
	*b = t;
}
int partition(Vector2 arr, int low, int high) {
	num pivot = arr[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++) {
		if (arr[j] < pivot) {
			i++;
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}
void quickSort(Vector2 arr, int low, int high) {
	if (low < high) {
		int pi = partition(arr, low, high);
		quickSort(arr, low, pi - 1);
		quickSort(arr, pi + 1, high);
	}
}
void merge(Vector2 arr, int l, int m, int r) {
	int i, j, k;
	int n1 = m - l + 1;
	int n2 = r - m;

	/* create temp arrays */
	Vector2 L = new num[n1];
	Vector2 R = new num[n2];

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++)
		L[i] = arr[l + i];
	for (j = 0; j < n2; j++)
		R[j] = arr[m + 1 + j];

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0;
	j = 0;
	k = l;
	while (i < n1 && j < n2) {
		if (L[i] <= R[j]) {
			arr[k] = L[i];
			i++;
		}
		else {
			arr[k] = R[j];
			j++;
		}
		k++;
	}

	/* Copy the remaining elements of L[], if there are any */
	while (i < n1)
	{
		arr[k] = L[i];
		i++;
		k++;
	}

	/* Copy the remaining elements of R[], if there are any */
	while (j < n2)
	{
		arr[k] = R[j];
		j++;
		k++;
	}
	delete[] L;
	delete[] R;
}
void mergeSort(Vector2 arr, int n) {
	int curr_size;
	int left_start;

	for (curr_size = 1; curr_size <= n - 1; curr_size = 2 * curr_size) {
		for (left_start = 0; left_start < n - 1; left_start += 2 * curr_size) {
			int mid = min(left_start + curr_size - 1, n - 1);
			int right_end = min(left_start + 2 * curr_size - 1, n - 1);
			merge(arr, left_start, mid, right_end);
		}
	}
}
Vector2 arrayRankTransform(Vector2 arr, pos size) {
	unordered_map<int, int> m;
	Vector2 unsorted_array = new num[size];
	for (int i = 0; i < size; i++) { unsorted_array[i] = -abs(arr[i]); }
	mergeSort(unsorted_array, size);
	int rank = 1;

	for (int i = 0; i < size; i++) {
		if (m[unsorted_array[i]] == 0) {
			m[unsorted_array[i]] = rank;
			rank++;
		}

	}
	Vector2 res = new num[size];
	for (int i = 0; i < size; i++) {
		res[i] = (m[arr[i]]);
	}
	delete[] unsorted_array;
	return res;
}
void buildMaxHeap(Vector2 arr, int n)
{
	for (int i = 1; i < n; i++)
	{
		// if child is bigger than parent 
		if (arr[i] > arr[(i - 1) / 2])
		{
			int j = i;

			// swap child and parent until 
			// parent is smaller 
			while (arr[j] > arr[(j - 1) / 2])
			{
				swap(arr[j], arr[(j - 1) / 2]);
				j = (j - 1) / 2;
			}
		}
	}
}

void heapSort(Vector2 arr, int n)
{
	buildMaxHeap(arr, n);

	for (int i = n - 1; i > 0; i--)
	{
		// swap value of first indexed  
		// with last indexed  
		swap(arr[0], arr[i]);

		// maintaining heap property 
		// after each swapping 
		int j = 0, index;

		do
		{
			index = (2 * j + 1);

			// if left child is smaller than  
			// right child point index variable  
			// to right child 
			if (arr[index] < arr[index + 1] &&
				index < (i - 1))
				index++;

			// if parent is smaller than child  
			// then swapping parent with child  
			// having higher value 
			if (arr[j] < arr[index] && index < i)
				swap(arr[j], arr[index]);

			j = index;

		} while (index < i);
	}
}

Vector2 Neg_Abs_R_Rank(Vector2 toRank, pos size) {
	Vector2 sorted = new num[size];
	Vector2 toReturn = new num[size];
	for (int i = 0; i < size; i++) { sorted[i] = -abs(toRank[i]); }
	//return arrayRankTransform(sorted, size);
	//mergeSort(sorted, size);
	//quickSort(sorted, 0, size-1);
	heapSort(sorted, size);
	//print(sorted, size);
	for (int i = 0; i < size; i++) {
		int value = i + 1;
		int repeats = 1;
		while (i + 1 < size && sorted[i] == sorted[i + 1]) {
			i++;
			value += i + 1;
			repeats++;
		}
		for (int j = 0; j < size; j++) {
			int num_inserts = 0;
			if (sorted[i] == -abs(toRank[j])) { toReturn[j] = num(value) / num(repeats); num_inserts++; }
			if (num_inserts == repeats) { j = size + 1; }
		}
	}
	delete[] sorted;
	return toReturn;
}

Matrix2 convert_rowwise(Matrix2 input, pos n, pos p) {
	Matrix2 newMatrix = new num * [p];
	for (int i = 0; i < p; i++) {
		newMatrix[i] = new num[n];
		for (int j = 0; j < n; j++) {
			newMatrix[i][j] = input[j][i];
		}
	}
	for (int i = 0; i < n; i++) {
		delete[] input[i];
	}
	delete[] input;
	return newMatrix;
}

//OPTIONAL FOR TESTING
Matrix readMatrix_csv(string filename, bool header) {
	Matrix toReturn;
	Vector toReturnV;
	fstream file;
	file.open(filename.c_str(), ios::in);
	string newLine;
	unsigned lastcomma;
	string subs;
	pos len;
	if (!file.is_open()) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Invalid file name."); }
	if (header) {
		getline(file, newLine);
	}
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
		if (toReturnV.size() != 0) {
			toReturn.push_back(toReturnV);
		}
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