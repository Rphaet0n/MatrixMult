#include<ctime>
#include<random>
#include<limits>
#include<string>
#include <chrono>
#include <ratio>
#include "mm.h"

#define N 2880

int main() {
	using namespace std::chrono;
	std::default_random_engine re;
	size_t sizes[21] = { 1, 6, 10, 15, 20, 24, 30, 36, 40, 60, 72, 80, 96, 120, 144, 160, 180, 240, 360, 480, 720 };	
	for (int i = 0; i < 21; ++i){
		genSimAsLowerTrMatrix<double>(("dA" + std::to_string(sizes[i])).c_str(), N, sizes[i], re);
		genUpperTrMatrix<double>(("dB" + std::to_string(sizes[i])).c_str(), N, sizes[i], re);
		genSimAsLowerTrMatrix<float>(("fA" + std::to_string(sizes[i])).c_str(), N, sizes[i], re);
		genUpperTrMatrix<float>(("fB" + std::to_string(sizes[i])).c_str(), N, sizes[i], re);
	}
	omp_set_num_threads(4);
	std::cout << "\n##################### Double #####################\n";
	for (int i = 0; i < 21; ++i) {

		TriangleMatrix<double> dA(("dA" + std::to_string(sizes[i])).c_str());
		TriangleMatrix<double> dB(("dB" + std::to_string(sizes[i])).c_str());
		Matrix<double> cLin(N, sizes[i]), cParInner(N, sizes[i]), cParOuter(N, sizes[i]);
		steady_clock::time_point t1 = steady_clock::now();
		multMatrix<double>(dA, dB, cLin);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> linTime = duration_cast<duration<double>>(t2 - t1);
		std::cout << "linTime: " << linTime.count() * 1000 << " ms ";

		t1 = steady_clock::now();
		multMatrixParInternal<double>(dA, dB, cParInner);
		t2 = steady_clock::now();
		duration<double> inParTime = duration_cast<duration<double>>(t2 - t1);
		std::cout << "inParTime: " << inParTime.count() * 1000 << " ms ";

		t1 = steady_clock::now();
		multMatrixParExternal<double>(dA, dB, cParOuter);
		t2 = steady_clock::now();
		duration<double> outParTime = duration_cast<duration<double>>(t2 - t1);
		std::cout << "outParTime: " << outParTime.count() * 1000 << " ms ";

		bool correct = checkMatrices<double>(cLin.getData(), cParInner.getData(), cParOuter.getData(), N*N);
		std::cout << "bsize: " << sizes[i] << " correct: " << (correct? "True":"False") << " Double\n";

	}

	std::cout << "\n##################### Float #####################\n";
	for (int i = 0; i < 21; ++i) {
		TriangleMatrix<float> fA(("fA" + std::to_string(sizes[i])).c_str());
		TriangleMatrix<float> fB(("fB" + std::to_string(sizes[i])).c_str());
		Matrix<float> cLin(N, sizes[i]), cParInner(N, sizes[i]), cParOuter(N, sizes[i]);

		steady_clock::time_point t1 = steady_clock::now();
		multMatrix<float>(fA, fB, cLin);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> linTime = duration_cast<duration<double>>(t2 - t1);
		std::cout << "linTime: " << linTime.count() * 1000 << " ms ";

		t1 = steady_clock::now();
		multMatrixParInternal<float>(fA, fB, cParInner);
		t2 = steady_clock::now();
		duration<double> inParTime = duration_cast<duration<double>>(t2 - t1);
		std::cout << "inParTime: " << inParTime.count() * 1000 << " ms ";

		t1 = steady_clock::now();
		multMatrixParExternal<float>(fA, fB, cParOuter);
		t2 = steady_clock::now();
		duration<double> outParTime = duration_cast<duration<double>>(t2 - t1);
		std::cout << "outParTime: " << outParTime.count() * 1000 << " ms ";

		bool correct = checkFloatMatrices(cLin.getData(), cParInner.getData(), cParOuter.getData(), N*N);
		std::cout << "bsize: " << sizes[i] << " correct: " << (correct ? "True" : "False") << " Float\n";
	}

}

