#include <iostream>
#include <cmath>

#include <ippcore.h>
#include <ipps.h>

int main(int ac, char* av[]) {
  //Set the size
  const int N = 128;
  const int order = (int)(log((double)N) / log(2.0));

  // Spec and working buffers
  IppsFFTSpec_C_32fc *pFFTSpec=0;
  Ipp8u *pFFTSpecBuf, *pFFTInitBuf, *pFFTWorkBuf;

  // Allocate complex buffers
  Ipp32fc *pSrc=ippsMalloc_32fc(N);
  Ipp32fc *pDst=ippsMalloc_32fc(N);

  // Query to get buffer sizes
  int sizeFFTSpec,sizeFFTInitBuf,sizeFFTWorkBuf;
  ippsFFTGetSize_C_32fc(order, IPP_FFT_NODIV_BY_ANY,
    ippAlgHintAccurate, &sizeFFTSpec, &sizeFFTInitBuf, &sizeFFTWorkBuf);

  // Alloc FFT buffers
  pFFTSpecBuf = ippsMalloc_8u(sizeFFTSpec);
  pFFTInitBuf = ippsMalloc_8u(sizeFFTInitBuf);
  pFFTWorkBuf = ippsMalloc_8u(sizeFFTWorkBuf);

  // Initialize FFT
  ippsFFTInit_C_32fc(&pFFTSpec, order, IPP_FFT_NODIV_BY_ANY,
    ippAlgHintAccurate, pFFTSpecBuf, pFFTInitBuf);
  if (pFFTInitBuf)
    ippFree(pFFTInitBuf);

  // Do the FFT
  ippsFFTFwd_CToC_32fc(pSrc,pDst,pFFTSpec,pFFTWorkBuf);

  //check results
  ippsFFTInv_CToC_32fc(pDst,pDst,pFFTSpec,pFFTWorkBuf);
  int OK = 1;
  for (int i = 0; i < N; i++) {
     pDst[i].re /= (Ipp32f)N;
     pDst[i].im /= (Ipp32f)N;
     if ((abs(pSrc[i].re - pDst[i].re) > .001) ||
         (abs(pSrc[i].im - pDst[i].im) > .001) )
    {
      OK = 0;
      break;
    }
  }

  std::cout << "FFT " << (1 == OK ? "OK" : "Fail") << std::endl;

  if (pFFTWorkBuf)
    ippFree(pFFTWorkBuf);
  if (pFFTSpecBuf)
    ippFree(pFFTSpecBuf);

  ippFree(pSrc);
  ippFree(pDst);

  return 0;
}
