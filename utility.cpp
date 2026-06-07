#include <complex>
#include <iostream>
#include <vector>
#include <time.h>       /* time */
#include<numeric>
#include <math.h>

// #include <iomanip>
// #include <ctime>
#include <Rcpp.h>


using Rcpp::Rcout;


const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d__%X", &tstruct);

    return buf;
}

const std::string currentDateTime(double time_offset) {
     time_t now = time( NULL);
    struct tm now_tm = *localtime( &now);
    struct tm then_tm = now_tm;
    then_tm.tm_sec += time_offset;   // add 50 seconds to the time
    mktime( &then_tm);
    char       buf[80];
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d__%X", &then_tm);
    return buf;
}

double calculate_stddev(std::vector<double>& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stddev = std::sqrt(sq_sum / v.size() - mean * mean);
    return stddev;
}

// 
// NumericMatrix calculate_amplitude(vector<std::complex<float>> &x, double &dt, double &fcenter){
//     
// }



// int calculate_FFT()
// {
//     const int nfft=32;
//     kiss_fft_cfg fwd = kiss_fft_alloc(nfft,0,NULL,NULL);
//     //kiss_fft_cfg inv = kiss_fft_alloc(nfft,1,NULL,NULL);

//     //vector<std::complex<float>> x(nfft, 0.0);
//     //vector<std::complex<float>> fx(nfft, 0.0);

//     x[0] = 1;
//     x[1] = std::complex<float>(0,3);

//     kiss_fft(fwd,(kiss_fft_cpx*)&x[0],(kiss_fft_cpx*)&fx[0]);
//     for (int k=0;k<nfft;++k) {
//         fx[k] = std::norm(fx[k]);
//         fx[k] *= 1./nfft;
//     }
//     //kiss_fft(inv,(kiss_fft_cpx*)&fx[0],(kiss_fft_cpx*)&x[0]);
//     //cout << "the circular correlation of [1, 3i, 0 0 ....] with itself = ";
//     cout
//         << fx[0] << ","
//         << fx[1] << ","
//         << fx[2] << ","
//         << fx[3] << " ... " << endl;
//     double mysum = 0;
//     for (int k=0;k<nfft;++k) {
//         mysum += real(fx[k]);
//     }

//     cout << "mysum: " << mysum << "\n";

//     kiss_fft_free(fwd);
//     kiss_fft_free(inv);
//     return 0;
// }