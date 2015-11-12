#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>

#include <boost/program_options.hpp>

#include <ippcore.h>
#include <ipps.h>

namespace       po = boost::program_options;
using namespace std::chrono;

const char*   _data_file_name = "fft-data.txt";
const char*   _fft_file_name  = "fft-forward.txt";
const char*   _bak_file_name  = "fft-backward.txt";

size_t        _fft_size       = 8192;
size_t        _batch_size     = 1000;
bool          _use_periodic   = false;
double        _mean           = 0.5;
double        _std            = 0.2;
long          _iterations     = 100000;

void populate_periodic(Ipp32f* buf) {
  for (int i = 0; i < _fft_size; ++i) {
    float t = i * .002;
    float amp = sin(M_PI * t);
    amp += sin(2 * M_PI * t);
    amp += sin(3 * M_PI * t);
    buf[i] = amp;
  }
}

void populate_random(Ipp32f* buf) {
  std::default_random_engine       generator(std::random_device{}());
  std::normal_distribution<float> distribution(_mean, _std);
  for (int i = 0; i < _fft_size; ++i) {
    buf[i] = distribution(generator);
  }
}

void populate(Ipp32f* buf) {
  if (_use_periodic)
    populate_periodic(buf);
  else
    populate_random(buf);
}

void write_real(std::string filename, Ipp32f* buf) {
  std::ofstream ofs;
  ofs.open(filename);
  ofs.precision(10);

  for (int i = 0; i < _fft_size; ++i) {
    ofs << buf[i] << std::endl;
  }

  ofs.close();
}

void write_pack(std::string filename, Ipp32f* buf) {
  std::ofstream ofs;
  ofs.open(filename);

  for (int i = 1; i < _fft_size / 2; i += 2) {
    float real = buf[i];
    float imag = buf[i+1];
    float amp = sqrt(pow(real, 2) + pow(imag, 2));
    float phase = atan2(imag, real);
    ofs << amp << ", " << phase << std::endl;
  }

  ofs.close();
}

std::string check_result(Ipp32f* pSrc, Ipp32f* pInv) {
  int OK = 1;
  for (int i = 0; i < _fft_size; i++) {
    if (0.001 < (abs(pSrc[i] - pInv[i]))) {
      OK = 0;
      break;
    }
  }
  return 1 == OK ? "OK" : "Fail";
}

double sqer(Ipp32f* pSrc, Ipp32f* pInv) {
  // signal power
  double sp = 0;
  for (int i = 0; i < _fft_size; i++) {
    sp += pow(pSrc[i], 2);
  }

  // quant error energy
  double qe = 0;
  for (int i = 0; i < _fft_size; i++) {
    qe += pow(pSrc[i] - pInv[i], 2);
  }

  // sqer
  return 10.0 * log(sp/qe);
}

void report_settings(int count) {

  const int order = (int)(log((double)_fft_size) / log(2.0));

  if (0 != count)
    std::cout << "Iterations:       " << count << std::endl;
  std::cout << "Data size:        " << _fft_size << std::endl;
  std::cout << "Data order:       " << order << std::endl;
  std::cout << "Data type:        " << (_use_periodic ? "Periodic" : "Random") << std::endl;
  if (!_use_periodic) {
    std::cout << "Mean:             " << _mean << std::endl;
    std::cout << "STD:              " << _std << std::endl;
  }
}

void test_fft() {
  //Set the size
  const int order = (int)(log((double)_fft_size) / log(2.0));

  // Spec and working buffers
  IppsFFTSpec_R_32f*  pFFTSpec = 0;
  Ipp8u*              pFFTSpecBuf;
  Ipp8u*              pFFTInitBuf;
  Ipp8u*              pFFTWorkBuf;

  // status
  IppStatus status;

  // Allocate buffers
  Ipp32f* pSrc = ippsMalloc_32f(_fft_size); // real
  Ipp32f* pDst = ippsMalloc_32f(_fft_size); // real, pack format
  Ipp32f* pInv = ippsMalloc_32f(_fft_size); // real, inverted data

  // populate the buffers
  populate(pSrc);
  write_real(_data_file_name, pSrc);

  // Query to get buffer sizes
  int sizeFFTSpec,sizeFFTInitBuf,sizeFFTWorkBuf;
  ippsFFTGetSize_R_32f(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintAccurate, &sizeFFTSpec, &sizeFFTInitBuf, &sizeFFTWorkBuf);

  // Alloc FFT buffers
  pFFTSpecBuf = ippsMalloc_8u(sizeFFTSpec);
  pFFTInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
  pFFTWorkBuf = ippsMalloc_8u(sizeFFTWorkBuf);

  // Initialize FFT
  status = ippsFFTInit_R_32f(&pFFTSpec, order, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate, pFFTSpecBuf, pFFTInitBuf);
  if (0 != status)
    std::cout << "FFT init status " << status << std::endl;
  if (pFFTInitBuf)
    ippFree(pFFTInitBuf);

  // Do the FFT
  status = ippsFFTFwd_RToPack_32f(pSrc, pDst, pFFTSpec, pFFTWorkBuf);
  if (0 != status)
    std::cout << "FFT status " << status << std::endl;

  write_pack(_fft_file_name, pDst);

  // Invert the FFT
  status = ippsFFTInv_PackToR_32f(pDst, pInv, pFFTSpec, pFFTWorkBuf);
  if (0 != status)
    std::cout << "iFFT status " << status << std::endl;

  // Write results
  write_real(_bak_file_name, pInv);

  //check results
  report_settings(0);
  std::cout << "FFT:              " << check_result(pSrc, pInv) << std::endl;
  std::cout << "SQER:             " << sqer(pSrc, pInv) << std::endl;

  // free
  if (pFFTWorkBuf)
    ippFree(pFFTWorkBuf);
  if (pFFTSpecBuf)
    ippFree(pFFTSpecBuf);

  ippFree(pSrc);
  ippFree(pDst);
  ippFree(pInv);
}

void time_fft() {
  //Set the size
  const int order = (int)(log((double)_fft_size) / log(2.0));

  // Spec and working buffers
  IppsFFTSpec_R_32f*  pFFTSpec = 0;
  Ipp8u*              pFFTSpecBuf;
  Ipp8u*              pFFTInitBuf;
  Ipp8u*              pFFTWorkBuf;

  // status
  IppStatus status;
  double    total_sqer = 0;

  // Allocate buffers
  std::vector<Ipp32f*> src = std::vector<Ipp32f*>(_batch_size);
  std::vector<Ipp32f*> dst = std::vector<Ipp32f*>(_batch_size);
  std::vector<Ipp32f*> inv = std::vector<Ipp32f*>(_batch_size);

  for (int i = 0; i < _batch_size; ++i) {
    src[i] = ippsMalloc_32f(_fft_size); // real
    dst[i] = ippsMalloc_32f(_fft_size); // real, pack format
    inv[i] = ippsMalloc_32f(_fft_size); // real, inverted data
  }

  // Query to get buffer sizes
  int sizeFFTSpec,sizeFFTInitBuf,sizeFFTWorkBuf;
  ippsFFTGetSize_R_32f(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintAccurate, &sizeFFTSpec, &sizeFFTInitBuf, &sizeFFTWorkBuf);

  // Alloc FFT buffers
  pFFTSpecBuf = ippsMalloc_8u(sizeFFTSpec);
  pFFTInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
  pFFTWorkBuf = ippsMalloc_8u(sizeFFTWorkBuf);

  // Initialize FFT
  status = ippsFFTInit_R_32f(&pFFTSpec, order, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate, pFFTSpecBuf, pFFTInitBuf);
  if (0 != status)
    std::cout << "FFT init status " << status << std::endl;
  if (pFFTInitBuf)
    ippFree(pFFTInitBuf);

  int last_percent = -1;
  nanoseconds duration(0);

  for (int l = 0; l < _iterations; ++l) {

    // Populate data
    for (int d = 0; d < _batch_size; ++d) {
      populate(src[d]);
    }

    // Start clock
    high_resolution_clock::time_point start = high_resolution_clock::now();

    // Do the batch of FFTs
    for (int i = 0; i < _batch_size; ++i) {
      // Do the FFT
      ippsFFTFwd_RToPack_32f(src[i], dst[i], pFFTSpec, pFFTWorkBuf);

      // Invert the FFT
      ippsFFTInv_PackToR_32f(dst[i], inv[i], pFFTSpec, pFFTWorkBuf);
    }

    // Stop clock
    high_resolution_clock::time_point finish = high_resolution_clock::now();
    duration += duration_cast<nanoseconds>(finish - start);

    // Compute the error
    for (int i = 0; i < _batch_size; ++i) {
      total_sqer += sqer(src[i], inv[i]);
    }

    // Report progress
    int percent = (int) ((double) l / (double) _iterations * 100.0);
    if (percent != last_percent) {
      std::cout << "\r" << percent << "%";
      std::cout.flush();
      last_percent = percent;
    }
  }

  // Report results
  int count = _batch_size * _iterations;
  double ave_dur = duration.count() / (double) (count * 2); // we do 2 FFTs per loop

  std::cout << "\r";
  report_settings(count);
  std::cout << "Total duration:   " << duration.count()  << " ns" << std::endl;
  std::cout << "Total SQER:       " << total_sqer << std::endl;
  std::cout << "Average duration: " << ave_dur << " ns (" << (ave_dur / 1000.0) << " μs)" << std::endl;
  std::cout << "Ave SQER:         " << (total_sqer / (double) count) << std::endl;

  // free
  if (pFFTWorkBuf)
    ippFree(pFFTWorkBuf);
  if (pFFTSpecBuf)
    ippFree(pFFTSpecBuf);

  for (int i = 0; i < _batch_size; ++i) {
    ippFree(src[i]);
    ippFree(dst[i]);
    ippFree(inv[i]);
  }
}

int main(int ac, char* av[])
{
  int  ret = 0;
  bool time = false;

  try {
    po::options_description desc("Allowed options");

    desc.add_options()
    ("help,h",         "Produce help message")

    ("size,s",         po::value<long>(), "Set the size of the buffer")

    ("periodic,p",     "Use a periodic data set")
    ("random,r",       "Use a gaussian distributed random data set")
    ("mean,m",         po::value<double>(), "Mean for random data")
    ("deviation,d",    po::value<double>(), "Standard deviation for random data")

    ("time,t",         "Time the FFT")
    ("iterations,i",   po::value<long>(), "Set the number of iterations to perform");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }

    if (vm.count("size")) {
      _fft_size = vm["size"].as<long>();
    }

    if (vm.count("periodic")) {
      _use_periodic = true;
    }

    if (vm.count("random")) {
    }

    if (vm.count("mean")) {
      _mean = vm["mean"].as<double>();
    }

    if (vm.count("deviation")) {
      _std = vm["deviation"].as<double>();
    }

    if (vm.count("time")) {
      time = true;
    }

    if (vm.count("iterations")) {
      _iterations = vm["iterations"].as<long>();
    }
    _iterations /= _batch_size;

  } catch (std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown error" << std::endl;
    return 1;
  }

  if (time)
    time_fft();
  else
    test_fft();

  return ret;
}
