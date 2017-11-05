#include<iostream>
#include<cstdio>
#include<functional>
#include<omp.h>

//начало блока для SLT и UT матрицы.
inline size_t getBlShift(size_t ib, size_t jb, size_t bsize, size_t blockmSize) {
	return bsize * bsize * ((2 * blockmSize - jb)*(jb + 1) / 2 - blockmSize + ib);
}

template <typename T>
void genSimAsLowerTrMatrix(const char * fname,  size_t msize, size_t bsize, std::default_random_engine &re) {
	size_t blockmSize = msize / bsize;	//блочный размер
	size_t blockCount = (blockmSize * blockmSize + blockmSize) / 2; //общее число блоков
	int blockElemCount = bsize * bsize; // объем блока
	size_t arrSize = blockCount * blockElemCount; //msize bsize datasize data
	FILE *matrixFile;
	if ((matrixFile = fopen(fname, "wb")) == NULL) {
		printf("Cannot open file.\n");
	}
	T *matrixAsArray = new T[arrSize];
	std::uniform_int_distribution<int> unifIntPart(-999, 999);
	std::uniform_int_distribution<int> unifFracPart(-100, 100);
	for (size_t jb = 0; jb < blockmSize - 1; ++jb) //блоки на диагонали в другом цикле
		for (size_t ib = jb + 1; ib < blockmSize; ++ib) {
			size_t blShift = getBlShift(ib, jb, bsize, blockmSize);
			T *curBlock = matrixAsArray + blShift;
			for (size_t i = 0; i < bsize; i++)
				for (size_t j = 0; j < bsize; j++)
					curBlock[i * bsize + j] = T(unifIntPart(re) + unifFracPart(re) / 100.0);
		}

	for (size_t ijb = 0; ijb < blockmSize; ++ijb) { //диагональ
		size_t blShift = getBlShift(ijb, ijb, bsize, blockmSize);
		T *curBlock = matrixAsArray + blShift;
		for (size_t i = 0; i < bsize; i++)
			for (size_t j = 0; j < bsize; j++)
				if (i > j)
					curBlock[i * bsize + j] = curBlock[j * bsize + i];
				else curBlock[i * bsize + j] = T(unifIntPart(re) + unifFracPart(re) / 100.0);
	}


	size_t sizes[] = { msize, bsize, arrSize };
	fwrite(&sizes, sizeof(size_t), 3, matrixFile);
	fwrite(matrixAsArray, sizeof(T), arrSize, matrixFile);
	fclose(matrixFile);
	delete[] matrixAsArray;
}

template <typename T>
void genUpperTrMatrix(const char * fname, size_t msize, size_t bsize, std::default_random_engine &re) {
	size_t blockmSize = msize / bsize;	//блочный размер
	size_t blockCount = (blockmSize * blockmSize + blockmSize) / 2; //общее число блоков
	int blockElemCount = bsize * bsize; // объем блока
	size_t arrSize = blockCount * blockElemCount; //msize bsize datasize data
	FILE *matrixFile;
	if ((matrixFile = fopen(fname, "wb")) == NULL) {
		printf("Cannot open file.\n");
	}
	T *matrixAsArray = new T[arrSize];
	std::uniform_int_distribution<int> unifIntPart(-999, 999);
	std::uniform_int_distribution<int> unifFracPart(-100, 100);
	for (size_t ib = 0; ib < blockmSize - 1; ++ib) //блоки на диагонали в другом цикле
		for (size_t jb = ib + 1; jb < blockmSize; ++jb) {
			size_t blShift = getBlShift(jb, ib, bsize, blockmSize);
			T *curBlock = matrixAsArray + blShift;
			for (size_t i = 0; i < bsize; i++)
				for (size_t j = 0; j < bsize; j++)
					curBlock[i * bsize + j] = T(unifIntPart(re) + unifFracPart(re) / 100.0);
		}

	for (size_t ijb = 0; ijb < blockmSize; ++ijb) { //диагональ
		size_t blShift = getBlShift(ijb, ijb, bsize, blockmSize);
		T *curBlock = matrixAsArray + blShift;
		for (size_t i = 0; i < bsize; i++)
			for (size_t j = 0; j < bsize; j++)
				if (i > j)
					curBlock[i * bsize + j] = 0;
				else curBlock[i * bsize + j] = T(unifIntPart(re) + unifFracPart(re) / 100.0);
	}


	size_t sizes[] = { msize, bsize, arrSize };
	fwrite(&sizes, sizeof(size_t), 3, matrixFile);
	fwrite(matrixAsArray, sizeof(T), arrSize, matrixFile);
	fclose(matrixFile);
	delete[] matrixAsArray;

}

template <typename T>
class Matrix {
public:
	T* data;
	size_t msize, bsize, blVolume, rowVolume;
	Matrix(size_t msize, size_t bsize) : msize(msize) {
		setBsize(bsize);
		data = (T*) calloc(msize * msize, sizeof(T));
	}

	void setBsize(size_t size) {
		bsize = size;
		blVolume = bsize * bsize;
		rowVolume = bsize * msize;
	} 
	

	void reset(size_t bsize) {
		free(data);
		setBsize(bsize);
	}

	T* getBlock(int i, int j) {
		return data + (rowVolume * i + j * blVolume);
	}

	T* getData() {
		return data;
	}

	void printBlocks() {
		for (int i = 0; i < msize*msize; ++i)
			if (i % blVolume)
				if (i % bsize)
					std::cout << data[i] << ' ';
				else
					std::cout << '\n' << data[i] << ' ';
			else
				std::cout << "\n\n" << data[i] << ' ';
	}

	~Matrix() {
		free(data);
	}
};

template <typename T>
class TriangleMatrix {
public:
	T* data;
	size_t msize, bsize, blockmSize, arrSize;
	TriangleMatrix(const char * fname) {
		data = nullptr;
		readFromFile(fname);
	}
	void readFromFile(const char * fname) {
		delete[] data;
		FILE *matrixFile;
		if (fopen_s(&matrixFile, fname, "rb") != 0) {
			printf("Cannot open file.\n"); //debug
		}
		size_t sizes[3];
		fread(&sizes, sizeof(size_t), 3, matrixFile);
		msize = sizes[0];
		bsize = sizes[1];
		blockmSize = msize / bsize;
		arrSize = sizes[2];
		if (arrSize > 0)
			data = new T[arrSize];
		fread(data, sizeof(T), arrSize, matrixFile);
	}

	//для A(SLT) - getBlock(i, j)
	//для B(UT) - getBlock(j, i)
	T* getBlock(int i, int j) {
		return  i <= j ? data + getBlShift(j, i, bsize, blockmSize) : data + getBlShift(i, j, bsize, blockmSize);
	}

	void printBlocks() {
		for (int i = 0; i < arrSize; ++i)
			if (i % bsize*bsize)
				if (i % bsize)
					std::cout << data[i] << ' ';
				else
					std::cout << '\n' << data[i] << ' ';
			else
				std::cout << "\n\n" << data[i] << ' ';
	}

	~TriangleMatrix() {
		delete[] data;
	}
};

template <typename T>
void multBlocks(T *a, T *b, T *c,  int &bi, int &bk, int &bj, int &bsize) {
	for (int i = 0; i < bsize; ++i)
		for (int j = 0; j < bsize; ++j) 
			for (int k = 0; k < bsize; ++k)
				c[i * bsize + j] += a[bi < bk ? k*bsize + i : i*bsize + k] * b[k*bsize + j];
}

template <typename T>
void multBlocksPar(T *a, T *b, T *c,  int &bi, int &bk, int &bj, int &bsize) {
#pragma omp parallel for
	 for (int i = 0; i < bsize; ++i)
		for (int j = 0; j < bsize; ++j) 
			for (int k = 0; k < bsize; ++k)
				c[i * bsize + j] += a[bi < bk ? k*bsize + i : i*bsize + k] * b[k*bsize + j];
				//c[i * bsize + j] += a[indA(i, k)] * b[k*bsize + j];//b[indB(k, j)];
}


template <typename T>
void multMatrix(TriangleMatrix<T> &a, TriangleMatrix<T> &b, Matrix<T> &c) {
	int blockmSize = a.blockmSize;
	int bSize = a.bsize;
	for (int i = 0; i < blockmSize; ++i)
		for (int j = 0; j < blockmSize; ++j)
			for (int k = 0; k < blockmSize; ++k){
				if (k > j)
					continue;
				multBlocks<T>(a.getBlock(i, k), b.getBlock(k, j), c.getBlock(i, j), i, k, j, bSize);
			}
}

template <typename T>
void multMatrixParInternal(TriangleMatrix<T> &a, TriangleMatrix<T> &b, Matrix<T> &c) {
	int blockmSize = a.blockmSize;
	int bSize = a.bsize;
	for (int i = 0; i < blockmSize; ++i)
		for (int j = 0; j < blockmSize; ++j)
			for (int k = 0; k < blockmSize; ++k){
				if (k > j)
					continue;
				multBlocksPar<T>(a.getBlock(i, k), b.getBlock(k, j), c.getBlock(i, j), i, k, j, bSize);
			}
}

template <typename T>
void multMatrixParExternal(TriangleMatrix<T> &a, TriangleMatrix<T> &b, Matrix<T> &c) {
	int blockmSize = a.blockmSize;
	int bSize = a.bsize;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < blockmSize; ++i)
		for (int j = 0; j < blockmSize; ++j)
			for (int k = 0; k < blockmSize; ++k){
				if (k > j)
					continue;
				multBlocks<T>(a.getBlock(i, k), b.getBlock(k, j), c.getBlock(i, j), i, k, j, bSize);
			}
}

template <typename T>
bool checkMatrices(T *a, T *b, T *c, int size) {
	std::cout << "\nCheck:";
	for (int i = 0; i < size; ++i)
		if (std::abs(a[i] - b[i]) + std::abs(b[i] - c[i]) > 2 * DBL_EPSILON)
			return false;
		return true;
}

bool flEqual(float A, float B,
	float maxRelDiff = 0.0001f) {
	float diff = fabs(A - B);
	A = fabs(A);
	B = fabs(B);
	float largest = (B > A) ? B : A;

	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}
bool checkFloatMatrices(float *a, float *b, float *c, int size) {
	std::cout << "\nCheck:";
	for (int i = 0; i < size; ++i)
		if (!flEqual(a[i], b[i]) || !flEqual(b[i], c[i]))
			return false;
	return true;
}



