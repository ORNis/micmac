#ifndef ARSENIC_H
#define ARSENIC_H
#include "StdAfx.h"

class ArsenicImage
{
public:
	ArsenicImage(){}
	~ArsenicImage(){}
	Im2D_REAL4 RChan;
	Im2D_REAL4 GChan;
	Im2D_REAL4 BChan;
	cElNuage3DMaille* info3D;
	Pt2di SZ;
};

class GrpVodka
{
public:
    GrpVodka(double diaph, double foc, bool isComputed){this->foc=foc;this->diaph=diaph;this->isComputed=isComputed;}
    ~GrpVodka(){}
    double diaph,foc;
    bool isComputed;
    vector<string> aListIm;
    vector<double> ExpTime, ISO;
    int size(){return (int)this->aListIm.size();}
};


class Param3Chan
{
public:
    Param3Chan(){}
    ~Param3Chan(){}
    std::vector<double> parRed, parBlue, parGreen;
    int size(){return (int)this->parRed.size();}
private:

};

class PtsHom
{
public:
    PtsHom(){}
    ~PtsHom(){}
    std::vector<double> Gr1, Gr2, Dist1, Dist2;
    Pt2di SZ;
    size_t size()
    {
      return Gr1.size();
    }
};

class cl_PtsRadio
{
public:
    cl_PtsRadio(){}
    ~cl_PtsRadio(){}
    vector<double> kR, kG, kB;
    std::vector<Pt2dr> Pts;
    std::vector<size_t> OtherIm;
    Pt2di SZ;
    size_t size()
    {
      return Pts.size();
    }
};

class cl_MatPtsHom
{
public:
    cl_MatPtsHom(){}
    ~cl_MatPtsHom(){}
    vector<cl_PtsRadio> aMat;

    size_t nbTotalPts()
    {
      int sum = 0;
      for(size_t i = 0; i < aMat.size(); ++i)
      {
        sum += aMat[i].size();
      }
      return sum;
    }
};

class cAppliCorrColor :  public cAppliWithSetImage
{
public:
	cAppliCorrColor(int argc, char ** argv);
  void rescaleImagesAndClouds();
	void doVignetteCorrection();
  void doRescaleData();
  void loadGrpImagesMMByP();
  void readPtsHom3D(const double &TPA);
  void computeNeighborsFromChImSec();
  void egalFieldCorrect(int nbIte, double aThresh);
  bool haveVignetteEqual() {return mVignetteDir != "";};

  std::string fullPathOrigRGBFile(size_t idx)
  {
    std::string postFix = "";
    if(haveVignetteEqual()) {
      postFix = "_Vodka.tif";
      return Dir() + ELISE_STR_DIR + mVignetteDir + ELISE_STR_DIR + mVectIm->at(idx) + postFix;
    } else {
      return Dir() + ELISE_STR_DIR + mVectIm->at(idx);
    }
  }
  std::string fullPathOrigCloud(size_t idx)
  {
    return Dir() + ELISE_STR_DIR + mMMIN->NameFileXml(eTMIN_Depth, mVectIm->at(idx));
  }

  std::string fullPathArsenicRGBFile(size_t idx)
  {
    if(mAreImagesDS) return Dir() + ELISE_STR_DIR + NameTmpDirArsenic() + ELISE_STR_DIR + mVectIm->at(idx) + ".tif";
    else fullPathOrigRGBFile(idx);
  }

  std::string fullPathArsenicCloud(size_t idx)
  {
    if(mAreCloudsDS) return Dir() + ELISE_STR_DIR + NameTmpDirArsenic()
                            + ELISE_STR_DIR + cElFilename(mMMIN->NameFileXml(eTMIN_Depth, mVectIm->at(idx))).m_basename;
    else fullPathOrigCloud(idx);
  }

  std::string fullPathCorrectedRGBFile(size_t idx)
  {
    Dir() + ELISE_STR_DIR + mOutDir + ELISE_STR_DIR + mVectIm->at(idx) + "_egal.tif";
  }

  static const std::string NameTmpDirArsenic(){return "Tmp-Arsenic";};

private:
	eTypeMMByP mModeMMByP;
  cInterfChantierNameManipulateur * mICNM;
  const std::vector<std::string> * mVectIm;
	cMMByImNM * mMMIN;
  int mDS;
  bool mAreImagesDS;
  bool mAreCloudsDS;
  std::string mVignetteDir;
  std::string mOutDir;
  cl_MatPtsHom mMatPtsHom;
  std::map<size_t, std::list<size_t> > mMapNeighborsIm;
	std::vector<ArsenicImage> mVectArsenicImage;
};

#endif // ARSENIC_H
