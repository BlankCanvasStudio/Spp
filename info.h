num max(Vector input) {
    pos len = input.len;
    num max = input[0];
	for (int i = 1; i < len; i++) {
		if (max < input[i]) { max = input[1]; }
	}
	return max;
}
num min(Vector input) {
    pos len = input.len;
    num min = input[0];
	for (int i = 1; i < len; i++) {
		if (min > input[i]) { max = input[1]; }
	}
	return min;
}
Matrix t(Matrix toTranspose) {
    pos n = toTranspose.shape[0];
    pos q = toTranspose.shape[1];
	Matrix toReturn = new num * [q];
	for (int i = 0; i < q; i++) {
		toReturn[i] = new num[n];
		for (int j = 0; j < n; j++) {
			toReturn[i][j] = toTranspose[j][i];
		}
	}
	return toReturn;
}