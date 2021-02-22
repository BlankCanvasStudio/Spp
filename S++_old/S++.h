//R.hVector
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <chrono>
#include <stack>
#include <map>
#define MIN_NUM -pow(1.17, 38)
#define num double
#define pos unsigned
#define Matrix num**
#define Vector num*
using namespace std;

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

//Official Vector functions
void print(Matrix toPrint, pos n, pos p) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p; j++) {
			cout << toPrint[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
void print(Vector toPrint, pos n) {
	for (int i = 0; i < n; i++) {
		cout << toPrint[i]<<endl;
	}
}
num sum(Vector a, pos size) {
	num toReturn = 0;
	for (pos i = 0; i < size; i++) {
		toReturn += a[i];
	}
	return toReturn;
}
num mean(Vector a, pos size) {
	return sum(a, size) / num(size);
}
Vector mean(Matrix toMean, pos n, pos p) {
	Vector toReturn = new num[p];
	for (int i = 0; i < p; i++) {
		Vector temp = new num[n];
		for (int j = 0; j < n; j++) {
			temp[j] = toMean[j][i];
		}
		toReturn[i] = mean(temp, n);
		delete[] temp;
	}
	return toReturn;
}
num Variance(Vector a, pos size) {
	num meanA = mean(a, size);
	num total = 0.0;
	for (int i = 0; i < size; i++) {
		total += pow((a[i] - meanA), 2);
	}
	if (size - 1 == 0) { return 0; }
	return total / num(size - 1);
}
num std_dev(Vector a, pos size) {
	return sqrt(Variance(a, size));
}

Vector std_dev(Matrix a, pos n, pos p) {
	Vector toReturn = new num[p];
	Vector temp = new num[n];
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < n; j++) { temp[j] = a[j][i]; }
		toReturn[i] = std_dev(temp, n);
	}
	delete[] temp;
	return toReturn;
}
Matrix scale(Matrix X, pos n, pos p, bool center, Vector s, pos lenS) {
	Matrix input = new num * [n];
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
void swap(Vector a, Vector b) {
	num t = *a;
	*a = *b;
	*b = t;
}
int partition(Vector arr, int low, int high) {
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
void quickSort(Vector arr, int low, int high) {
	if (low < high) {
		int pi = partition(arr, low, high);
		quickSort(arr, low, pi - 1);
		quickSort(arr, pi + 1, high);
	}
}
void merge(Vector arr, int l, int m, int r) {
	int i, j, k;
	int n1 = m - l + 1;
	int n2 = r - m;

	/* create temp arrays */
	Vector L = new num[n1];
	Vector R = new num[n2];

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
void mergeSort(Vector arr, int n) {
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
void buildMaxHeap(Vector arr, int n)
{
	for (int i = 1; i < n; i++){
		if (arr[i] > arr[(i - 1) / 2]){
			int j = i;
 
			while (arr[j] > arr[(j - 1) / 2])	{
				swap(arr[j], arr[(j - 1) / 2]);
				j = (j - 1) / 2;
			}
		}
	}
}

void heapSort(Vector arr, int n){
	buildMaxHeap(arr, n);

	for (int i = n - 1; i > 0; i--){
		swap(arr[0], arr[i]);

		int j = 0, index;

		do{
			index = (2 * j + 1);

			if (index < (i - 1) && arr[index] < arr[index + 1])
				index++;

			if (index < i && arr[j] < arr[index])
				swap(arr[j], arr[index]);

			j = index;

		} while (index < i);
	}
}
Vector Neg_Abs_R_Rank(Vector toRank, pos size) {
	Vector sorted = new num[size];
	multimap<num, pos> locations;
	for (int i = 0; i < size; i++) {
		locations.insert(pair<num, pos>(-abs(toRank[i]), i));
		sorted[i] = -abs(toRank[i]); 
	}
	Vector toReturn = new num[size];
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
		multimap<num, pos>::const_iterator Values =locations.find(sorted[i]);
		int Number = locations.count(sorted[i]);
		for (int i = 0; i < Number; i++){
			toReturn[Values->second] = num(value) / num(repeats);
			++Values;
		}
	}
	delete[] sorted;
	return toReturn;
}


int binarySearch(Vector arr, int l, int r, num x)
{
	while (l <= r) {
		int m = l + (r - l) / 2;

		// Check if x is present at mid 
		if (arr[m] == x)
			return m;

		// If x greater, ignore left half 
		if (arr[m] < x)
			l = m + 1;

		// If x is smaller, ignore right half 
		else
			r = m - 1;
	}

	// if we reach here, then element was 
	// not present 
	return -1;
}
Matrix convert_rowwise(Matrix input, pos n, pos p) {
	Matrix newMatrix = new num * [p];
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

Matrix readMatrix_csv(string filename, pos n, pos p, bool header) {
	Matrix toReturn = new num * [n];
	for (int i = 0; i < n; i++) {
		toReturn[i] = new num[p];
	}
	pos index_row = 0;
	pos index_col = 0;
	fstream file;
	file.open(filename, ios::in);
	string newLine;
	unsigned lastcomma;
	string subs;
	pos len;
	if (!file.is_open()) { cout << "ERROR ENCOUNTERED" << endl; throw invalid_argument("Invalid file name."); }
	if (header) {
		getline(file, newLine);
	}
	while (!file.eof() && index_row < n) {
		getline(file, newLine);
		lastcomma = 0;
		len = newLine.length();
		for (unsigned i = 0; i < len; i++) {
			if (newLine[i] == ',' || i == len - 1) {
				toReturn[index_row][index_col] = stod(newLine.substr(lastcomma, i));
				index_col = (index_col + 1) % p;
				lastcomma = i + 1;
			}
		}
		index_col = 0;
		index_row++;
	}
	return toReturn;
}
Matrix matrix_full(num toFill, pos n, pos p) {
	Matrix toReturn = new num * [n];
	for (int i = 0; i < n; i++) {
		toReturn[i] = new num[p];
		for (int j = 0; j < p; j++) {
			toReturn[i][j] = toFill;
		}
	}
	return toReturn;
}
void diag(Matrix toFill, pos n, num val) {
	for (int i = 0; i < n; i++) {
		toFill[i][i] = val;
	}
}
Vector multiply(Matrix a, pos rowA, pos colA, Vector b, pos lenB) {
	Vector toReturn = new num[rowA];
	for (pos i = 0; i < rowA; i++) {
		toReturn[i] = 0.0;
		for (pos j = 0; j < colA; j++) {
			toReturn[i] += a[i][j] * b[j];
		}
	}
	return toReturn;
}
Matrix subtract(Matrix a, Matrix b, pos rowA, pos colA, pos rowB, pos colB) {
	if (rowA != rowB && colA != colB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("SIze of A and B must be the same."); }
	Matrix toReturn = new num * [rowA];
	for (int i = 0; i < rowA; i++) {
		toReturn[i] = new num[colA];
		for (int j = 0; j < colA; j++) {
			toReturn[i][j] = a[i][j] - b[i][j];
		}
	}
	return toReturn;
}
Matrix square(Matrix toSquare, pos n, pos p) {
	Matrix toReturn = new num * [n];
	for (int i = 0; i < n; i++) {
		toReturn[i] = new num[p];
		for (int j = 0; j < p; j++) {
			toReturn[i][j] = toSquare[i][j] * toSquare[i][j];
		}
	}
	return toReturn;
}
num arithmeticMean(Matrix toMean, pos rowToMean, pos colToMean) {
	num mean = 0;
	for (pos i = 0; i < rowToMean; i++) {
		for (pos j = 0; j < colToMean; j++) {
			mean += toMean[i][j];
		}
	}
	return mean / num(rowToMean * colToMean);
}
num mad(Vector toMad, pos len) {
	num absSum = 0;
	num meanVal = mean(toMad, len);
	for (int i = 0; i < len; i++) {
		absSum += abs(toMad[i] - meanVal);
	}
	return absSum / len;
}
Vector col_mad(Matrix toMad, pos row, pos col) {
	Vector toReturn = new num[col];
	Vector temp = new num[row];
	num mean;
	for (int i = 0; i < col; i++) {
		mean = 0;
		for (int j = 0; j < row; j++) {
			temp[j] = toMad[j][i];
		}
		toReturn[i] = mad(temp, row);
	}
	delete[] temp;
	return toReturn;
}
Matrix multiply(Matrix a, Matrix b, pos rowA, pos colA, pos rowB, pos colB) {
	if (colA != rowB) { cout << "AN ERROR OCCURED." << endl; throw invalid_argument("The column size of A and row size of B must be the same."); }
	Matrix result = new num * [rowA];
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
int* min_index(Matrix a, pos row, pos col) {
	num min_num = DBL_MAX;
	int* min_location = new int[2];
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
Matrix mean_over_3D(Matrix* toMean, pos TD, pos row, pos col) {
	Matrix toReturn = new num * [row];
	for (int i = 0; i < row; i++) {
		toReturn[i] = new num[col];
		for (int j = 0; j < col; j++) {
			toReturn[i][j] = 0.0;
			for (int k = 0; k < TD; k++) {
				toReturn[i][j] += toMean[k][i][j];
			}
			toReturn[i][j] /= TD;
		}
	}
	return toReturn;
}
Vector seq(num lower, num upper, num steps) {
	Vector toReturn = new num[steps];
	for (int i = 0; i < steps; i++) {
		toReturn[i] = lower + (i * ((upper - lower) / steps));
	}
	return toReturn;
}
num max(Vector input, pos len) {
	num max = input[0];
	for (int i = 1; i < len; i++) {
		if (max < input[i]) { max = input[1]; }
	}
	return max;
}
Matrix t(Matrix toTranspose, pos n, pos q) {
	Matrix toReturn = new num * [q];
	for (int i = 0; i < q; i++) {
		toReturn[i] = new num[n];
		for (int j = 0; j < n; j++) {
			toReturn[i][j] = toTranspose[j][i];
		}
	}
	return toReturn;
}
void saveAsCSV(Matrix Input, pos rowInput, pos columnInput, string filename) {
	if (filename.substr(filename.length() - 5) != ".csv") {
		filename.append(".csv");
	}
	string line;
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

void save_CV_return(cv_returnList out, string filename) {
	if (filename.substr(filename.length() - 5) != ".csv") {
		filename.append(".csv");
	}
	string line;
	ofstream outputFile;
	outputFile.open(filename);
	string toAdd = "";
	for (int i = 0; i < out.nlam; i++) {
		toAdd += " " + to_string(out.lambda[i]);
	}
	outputFile << "lambda:" + toAdd + "\n";

	toAdd = "";
	for (int i = 0; i < out.nr; i++) {
		toAdd += " " + to_string(out.vr[i]);
	}
	outputFile << "vr:" + toAdd + "\n";


	outputFile << "r_min: " + to_string(out.r_opt) + "\n";
	outputFile << "lambda_min: " + to_string(out.lambda_opt) + "\n";

	outputFile << "optloc:\n";
	for (int i = 0; i < out.nfold + 1; i++) {
		toAdd = "";
		for (int j = 0; j < 2; j++) {
			toAdd.append(to_string(out.optloc[i][j]));
			toAdd.append(",");
		}
		toAdd.append("\n");
		outputFile << toAdd;
	}
	outputFile << "cvm:\n";
	for (int i = 0; i < out.nlam; i++) {
		toAdd = "";
		for (int j = 0; j < out.nr; j++) {
			toAdd.append(to_string(out.cvm[i][j]));
			toAdd.append(",");
		}
		toAdd.append("\n");
		outputFile << toAdd;
	}
	outputFile.close();
}