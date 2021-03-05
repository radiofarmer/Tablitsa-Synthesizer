#pragma once

#include "vectorclass.h"

#include <initializer_list>

// Development macros
#define ALG2

template <typename T>
class ChebyshevBL
{
public:
  template <typename T>
  class Biquad
  {
  public:
    Biquad(std::initializer_list<double> coefs) : mCoefs(coefs)
    {
      mB = Vec4d(mCoefs[0], mCoefs[1], mCoefs[2], 0.);
      mA = Vec4d(0., 1., -1. * mCoefs[4], -1. * mCoefs[5]);

      // Recursive coefficient expressions
      T a1 = -1. * mCoefs[4];
      T a2 = -1. * mCoefs[5];
      T a_rec1[4]{
        1.,
        a1,
        a1 * a1 + a2,
        a1 * a2
      };
      T a_rec2[4]{
        a1,
        a1 * a1 + a2,
        a1 * a1 * a1 + 2 * a1 * a2,
        a1 * a1 * a2 + a2 * a2
      };
      T a_rec3[4]{
        a1 * a1 + a2,
        a1 * a1 * a1 + 2 * a1 * a2,
        a1 * a1 * a1 * a1 + 3 * a1 * a1 * a2 + a2 * a2,
        a1 * a1 * a1 * a2 + 2 * a1 * a2 * a2
      };
#if defined ALG1
      // Recursive a-coefficients (ALG1)
      mAr1.load(a_rec1);
      mAr2.load(a_rec2);
      mAr3.load(a_rec3);
#elif defined ALG2
      // Column vectors of coefficients (ALG2)
      mCol1 = Vec4d(a_rec3[0], a_rec2[0], a_rec1[0], mA[0]);
      mCol2 = Vec4d(a_rec3[1], a_rec2[1], a_rec1[1], mA[1]);
      mCol3 = Vec4d(a_rec3[2], a_rec2[2], a_rec1[2], mA[2]);
      mCol4 = Vec4d(a_rec3[3], a_rec2[3], a_rec1[3], mA[3]);

      // Other terms
      mA1 = Vec4d(-mCoefs[4], 1., 0., 0.);
      mA2 = Vec4d(1., 0., 0., 0.);
#endif
    }

    inline T Process(T s)
    {
      // Direct Form 1
      T sum = mCoefs[0] * s + mCoefs[1] * mX[0] + mCoefs[2] * mX[1] - mCoefs[4] * mY[0] - mCoefs[5] * mY[1];
      mX.push(s);
      mY.push(sum);
      return sum;
    }

    inline T ProcessDFII(T s)
    {
      // Direct Form 2 - summed last two terms in each expression are the values stored in mAddends
      T sum = s - mCoefs[4] * mZ[0] - mCoefs[5] * mZ[1]; // SUM1
      T out = mCoefs[0] * sum + mCoefs[1] * mZ[0] + mCoefs[2] * mZ[1]; // SUM2
      mZ.push(sum);
      return out;
    }

    inline T ProcessPrecomp(T s)
    {
      // Direct Form 2
      T sum = s - mAddends[0]; // SUM1 
      T out = mCoefs[0] * sum + mAddends[1]; // SUM2
      mMatrixDelay[1] = mMatrixDelay[0];
      mMatrixDelay[0] = sum;
      return out;
    }

    inline Vec4d __vectorcall ProcessRecursive(const Vec4d& x)
    {
      // x = (x[n], x[n+1], x[n+2], x[n+3]), i.e. earliest -> latest
      mDelayVector = blend4<1, 0, 4, 5>(x, mDelayVector); // x[n+1], x[n], w[n-1], w[n-2], i.e. latest -> earliest

#ifdef ALG1
      T w3 = x[3] - mCoefs[4] * x[2] + horizontal_add(mAr3 * mDelayVector); // w[n+3]
      T w2 = x[2] + horizontal_add(mAr2 * mDelayVector); // w[n+2]
      T w1 = horizontal_add(mAr1 * mDelayVector); // w[n+1]
      T w0 = horizontal_add(mA * mDelayVector); // w[n]
      Vec4d y_in = Vec4d(w3, w2, w1, w0); // latest -> earliest
#elif defined ALG2

      // All delay terms:
      Vec4d y_in = mA2 * x[3] + mA1 * x[2] + mCol1 * mDelayVector[0] + mCol2 * mDelayVector[1] + mCol3 * mDelayVector[2] + mCol4 * mDelayVector[3];
#endif

      const Vec4d y_out = Vec4d(
        horizontal_add(blend4<3, 6, 7, V_DC>(y_in, mDelayVector) * mB), // y[n] = b0 * w[n] + b1 * w[n-1] + b2 * w[n-2]
        horizontal_add(blend4<2, 3, 6, V_DC>(y_in, mDelayVector) * mB), // y[n+1] = b0 * w[n+1] + b1 * w[n] + b2 * w[n-1]
        horizontal_add(permute4<1, 2, 3, V_DC>(y_in) * mB),
        horizontal_add(y_in * mB) // y[n+3] = b0 * w[n+3] + b1 * w[n+2] + b2 * w[n+1]
      ); // earliest -> latest

      mDelayVector = permute4<0, 1, V_DC, V_DC>(y_in); // w[n+3] -> w[n-1]; w[n+2] -> w[n-2]; Order is latest -> earliest
      return y_out;
    }

    void SetAddendPtr(T* ptr)
    {
      mAddends = ptr;
    }

    void SetDelayPtr(T* ptr)
    {
      mMatrixDelay = ptr;
    }

  protected:
    const std::vector<double> mCoefs;
    DelayLine<T, 2> mZ;
    Vec4d mB;
    Vec4d mA;
#if defined ALG1
    Vec4d mAr3; // Recursive coefficients, third iteration
    Vec4d mAr2; // Recursive coefficients, second iteration
    Vec4d mAr1; // Recursive coefficients, first
#elif defined ALG2
    Vec4d mA1;
    Vec4d mA2;
    Vec4d mCol1;
    Vec4d mCol2;
    Vec4d mCol3;
    Vec4d mCol4;
#endif

    Vec4d mDelayVector = Vec4d(0.);

    T* mAddends = nullptr;
    T* mMatrixDelay = nullptr;

    friend class ChebyshevBL<T>;
  };

public:
  ChebyshevBL()
  {
    T z1[12]{};
    T z2[12]{};

    for (int i{ 0 }; i < 6; ++i)
    {
      // Provide each second-order filter with a pointer to the addend terms
      mSOS[i].SetAddendPtr(&mAddendPtrs[2 * i]);
      mSOS[i].SetDelayPtr(mDelayPtrs[i]);
      // Store the coefficients in this parent filter object
      z1[2 * i] = mSOS[i].mCoefs[4];
      z1[2 * i + 1] = mSOS[i].mCoefs[1];
      z2[2 * i] = mSOS[i].mCoefs[5];
      z2[2 * i + 1] = mSOS[i].mCoefs[2];
    }

    mZ1_Coefs1.load(z1);
    mZ1_Coefs2.load(z1 + 4);
    mZ1_Coefs3.load(z1 + 8);

    mZ2_Coefs1.load(z2);
    mZ2_Coefs2.load(z2 + 4);
    mZ2_Coefs3.load(z2 + 8);
  }

  inline void Precompute()
  {
    // Get values from the 6x2 delay matrix
    Vec4d flt01_taps = Vec4d().load(mSOS[0].mMatrixDelay);
    Vec4d flt23_taps = Vec4d().load(mSOS[2].mMatrixDelay);
    Vec4d flt45_taps = Vec4d().load(mSOS[4].mMatrixDelay);
    // Alternatively:
//  Vec4d z1_0 = blend<0, 0, 4, 4>(mSOS[0].mDelayVector, mSOS[1].mDelayVector);
    // etc...

    // Tap 1
    Vec4d z1_0 = permute4<0, 0, 2, 2>(flt01_taps);
    Vec4d z1_1 = permute4<0, 0, 2, 2>(flt23_taps);
    Vec4d z1_2 = permute4<0, 0, 2, 2>(flt45_taps);
    // Tap 2
    Vec4d z2_0 = permute4<1, 1, 3, 3>(flt01_taps);
    Vec4d z2_1 = permute4<1, 1, 3, 3>(flt23_taps);
    Vec4d z2_2 = permute4<1, 1, 3, 3>(flt45_taps);

    // First 4 terms
    Vec4d addends1 = (z1_0 * mZ1_Coefs1) + (z2_0 * mZ2_Coefs1);
    addends1.store(mAddendPtrs);
    // Next 4 terms
    Vec4d addends2 = (z1_1 * mZ1_Coefs2) + (z2_1 * mZ2_Coefs2);
    addends2.store(mAddendPtrs + 4);
    // Last 4 terms
    Vec4d addends3 = (z1_2 * mZ1_Coefs3) + (z2_2 * mZ2_Coefs3);
    addends3.store(mAddendPtrs + 8);
  }

  inline Vec4d __vectorcall ProcessAndDownsample_Recursive(Vec4d& x)
  {
    for (int i{ 0 }; i < 6; ++i)
      x = mSOS[i].ProcessRecursive(x);
    Vec4d out = permute4<1, 3, V_DC, V_DC>(x);
    return out;
  }

  inline Vec4d __vectorcall ProcessAndDownsample_Recursive(Vec4d& x1, Vec4d& x2)
  {
    for (int i{ 0 }; i < 6; ++i)
      x1 = mSOS[i].ProcessRecursive(x1);
    for (int i{ 0 }; i < 6; ++i)
      x2 = mSOS[i].ProcessRecursive(x2);
    Vec4d out = blend4<1, 3, 5, 7>(x1, x2);
    return out;
  }

  inline T ProcessAndDownsample(T* s)
  {
    Process(s[0]); // Process and throw away sample
    return Process(s[1]); // Process and return sample
  }

  inline T Process(T s)
  {
    T sum = s;
    Precompute();
    sum = mSOS[0].ProcessPrecomp(sum);
    sum = mSOS[1].ProcessPrecomp(sum);
    sum = mSOS[2].ProcessPrecomp(sum);
    sum = mSOS[3].ProcessPrecomp(sum);
    sum = mSOS[4].ProcessPrecomp(sum);
    sum = mSOS[5].ProcessPrecomp(sum);
    return sum;
  }

private:
  Biquad<T> mSOS[6]{ { 1.0194978892532574e-05, 2.0389957785065147e-05, 1.0194978892532574e-05, 1.0, -1.6759049412919762, 0.7386853939104704 },
    { 1.0, 2.0, 1.0, 1.0, -1.4134265459338102, 0.7727803769380925 },
    { 1.0, 2.0, 1.0, 1.0, -1.0347899844575996, 0.8248787327348859 },
    { 1.0, 2.0, 1.0, 1.0, -0.6862027551047557, 0.8794409054482396 },
    { 1.0, 2.0, 1.0, 1.0, -0.4401516515832734, 0.9299883719435977 },
    { 1.0, 2.0, 1.0, 1.0, -0.317836352796945, 0.9768802444563486 } };

  Vec4d mZ1_Coefs1, mZ1_Coefs2, mZ1_Coefs3;
  Vec4d mZ2_Coefs1, mZ2_Coefs2, mZ2_Coefs3;

  T mAddendPtrs[12]{ 0. };
  T mDelayPtrs[6][2]{ 0. };
};