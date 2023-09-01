#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>


typedef std::complex<double> base;


// sources
// https://en.wikipedia.org/wiki/Fast_Fourier_transform
// https://algorithmica.org/ru/fft


// FFT
// Time Complexity: O(nlogn)
// Space Complexity: O(n) (can be improved to O(1))
class fft_converter{
public:
    void fft(std::vector<base> & a){
        int size = a.size();
        a.resize(to_power(size));
        fft(a, false);
        a.resize(size);
    }
    void inverted_fft(std::vector<base> &a){
        int size = a.size();
        a.resize(to_power(size));
        fft(a,true);
        a.resize(size);
    }

private:
    void fft (std::vector<base> & a, bool invert) {
        int n = static_cast<int>(a.size());
        if (n == 1)  return;
        std::vector<base> a0 (n/2),  a1 (n/2);
        for (int i=0, j=0; i<n; i+=2, ++j) {
            a0[j] = a[i];
            a1[j] = a[i+1];
        }
        fft (a0, invert);
        fft (a1, invert);

        double ang = 2*M_PI/n * (invert ? -1 : 1);
        base w (1),  wn (cos(ang), sin(ang));
        for (int i=0; i<n/2; ++i) {
            a[i] = a0[i] + w * a1[i];
            a[i+n/2] = a0[i] - w * a1[i];
            if (invert)
                a[i] /= 2,  a[i+n/2] /= 2;
            w *= wn;
        }
    }

    size_t to_power(int size){
        size_t n = 1;
        while (n < size) n <<= 1;
        return n;
    }
};




int main() {

    // Generating random positive numbers
    std::mt19937 rand;
    rand.seed(std::random_device()());

    // From 1 to 100
    std::uniform_int_distribution<int> distribution(1, 100);

    const int AMOUNT = 100;

    std::vector<double> error_amounts(AMOUNT);

    const int SIZE = 8; // FFT size

    fft_converter fft;


    for(int i = 0; i < AMOUNT;++i){
        std::vector<base> input_data;
        input_data.reserve(SIZE);

        for (int i = 0; i < SIZE; ++i) {
            input_data.emplace_back(distribution(rand), distribution(rand));
        }

        std::vector<base> temp(input_data);

        fft.fft(temp);
        fft.inverted_fft(temp);

        // Calculating error
        double error = 0.0;
        for (int i = 0; i < SIZE; ++i) {
            error += std::abs(input_data[i].real() - temp[i].real()) +
                     std::abs(input_data[i].imag() - temp[i].imag());
        }
        error /= SIZE;
        error_amounts[i] = error;
    }

    // Print results


    std::cout << "Average error on each iteration:";
    double average = 0;
    for(int i = 0; i < AMOUNT;++i){
        average += error_amounts[i];
        std::cout << "Iteration " << i << ": " << error_amounts[i] << "\n";
    }

    average /= AMOUNT;

    std::cout << "Average error from all iterations: " << average;
}
