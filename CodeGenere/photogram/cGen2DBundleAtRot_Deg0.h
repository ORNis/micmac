// File Automatically generated by eLiSe
#include "StdAfx.h"


class cGen2DBundleAtRot_Deg0: public cElCompiledFonc
{
   public :

      cGen2DBundleAtRot_Deg0();
      void ComputeVal();
      void ComputeValDeriv();
      void ComputeValDerivHessian();
      double * AdrVarLocFromString(const std::string &);
      void SetDepR1_x(double);
      void SetDepR1_y(double);
      void SetDepR2_x(double);
      void SetDepR2_y(double);
      void SetDepR3_x(double);
      void SetDepR3_y(double);


      static cAutoAddEntry  mTheAuto;
      static cElCompiledFonc *  Alloc();
   private :

      double mLocDepR1_x;
      double mLocDepR1_y;
      double mLocDepR2_x;
      double mLocDepR2_y;
      double mLocDepR3_x;
      double mLocDepR3_y;
};
