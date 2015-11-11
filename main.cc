#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

#include <boost/program_options.hpp>

#include <ippcore.h>
#include <ipps.h>

namespace     po = boost::program_options;

const char*   _data_file_name = "fft-data.txt";
const char*   _fft_file_name  = "fft-forward.txt";
const char*   _bak_file_name  = "fft-backward.txt";

size_t        _fft_size       = 8192;
bool          _use_periodic   = false;
double        _mean           = 0.5;
double        _std            = 0.2;
long          _iterations     = 1000000;

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
  Ipp32f *pSrc = ippsMalloc_32f(_fft_size);  // real
  Ipp32f *pDst = ippsMalloc_32f(_fft_size);  // real, pack format

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
  status = ippsFFTInit_R_32f(&pFFTSpec, order, IPP_FFT_NODIV_BY_ANY, ippAlgHintAccurate, pFFTSpecBuf, pFFTInitBuf);
  std::cout << "FFT init status " << status << std::endl;
  if (pFFTInitBuf)
    ippFree(pFFTInitBuf);

  // Do the FFT
  status = ippsFFTFwd_RToPack_32f(pSrc, pDst, pFFTSpec, pFFTWorkBuf);
  std::cout << "FFT status " << status << std::endl;

  write_pack(_fft_file_name, pDst);

/*
  //check results
  ippsFFTInv_CToC_32fc(pDst, pDst, pFFTSpec, pFFTWorkBuf);
  int OK = 1;
  for (int i = 0; i < _fft_size; i++) {
     pDst[i].re /= (Ipp32f) _fft_size;
     pDst[i].im /= (Ipp32f) _fft_size;
     if ((abs(pSrc[i].re - pDst[i].re) > .001) ||
         (abs(pSrc[i].im - pDst[i].im) > .001) )
    {
      OK = 0;
      break;
    }
  }

  std::cout << "FFT " << (1 == OK ? "OK" : "Fail") << std::endl;
*/
  if (pFFTWorkBuf)
    ippFree(pFFTWorkBuf);
  if (pFFTSpecBuf)
    ippFree(pFFTSpecBuf);

  ippFree(pSrc);
  ippFree(pDst);

}

void time_fft() {

}

int main(int ac, char* av[])
{
  int  ret = 0;
  bool time = false;

  try {
    po::options_description desc("Allowed options");

    desc.add_options()
    ("help,h",         "Produce help message")

    ("size,s",         po::value<int>(), "Set the size of the buffer [8192]")

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
      _fft_size = vm["size"].as<int>();
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
