/*Header-MicMac-eLiSe-25/06/2007

    MicMac : Multi Image Correspondances par Methodes Automatiques de Correlation
    eLiSe  : ELements of an Image Software Environnement

    www.micmac.ign.fr

   
    Copyright : Institut Geographique National
    Author : Marc Pierrot Deseilligny
    Contributors : Gregoire Maillet, Didier Boldo.

[1] M. Pierrot-Deseilligny, N. Paparoditis.
    "A multiresolution and optimization-based image matching approach:
    An application to surface reconstruction from SPOT5-HRS stereo imagery."
    In IAPRS vol XXXVI-1/W41 in ISPRS Workshop On Topographic Mapping From Space
    (With Special Emphasis on Small Satellites), Ankara, Turquie, 02-2006.

[2] M. Pierrot-Deseilligny, "MicMac, un lociel de mise en correspondance
    d'images, adapte au contexte geograhique" to appears in 
    Bulletin d'information de l'Institut Geographique National, 2007.

Francais :

   MicMac est un logiciel de mise en correspondance d'image adapte 
   au contexte de recherche en information geographique. Il s'appuie sur
   la bibliotheque de manipulation d'image eLiSe. Il est distibue sous la
   licences Cecill-B.  Voir en bas de fichier et  http://www.cecill.info.


English :

    MicMac is an open source software specialized in image matching
    for research in geographic information. MicMac is built on the
    eLiSe image library. MicMac is governed by the  "Cecill-B licence".
    See below and http://www.cecill.info.

Header-MicMac-eLiSe-25/06/2007*/


#include "NewRechPH.h"

double TimeAndReset(ElTimer &aChrono )
{
   double aRes = aChrono.uval();
   aChrono.reinit();

   return aRes;
}

std::string NameFileNewPCarac(const std::string & aNameGlob,bool Bin,const std::string & anExt)
{
    std::string aDirGlob = DirOfFile(aNameGlob);
    std::string aDirLoc= "NewPH" + anExt + "/";
    ELISE_fp::MkDirSvp(aDirGlob+aDirLoc);
    return aDirGlob+aDirLoc + NameWithoutDir(aNameGlob) + (Bin ? ".dmp" : ".xml");
}


void  cAppli_NewRechPH::Clik()
{
   if (mW1) mW1->clik_in();
}

void cAppli_NewRechPH::AddScale(cOneScaleImRechPH * aI1,cOneScaleImRechPH * aI2)
{
    // std::cout << "---------PPppp   " << aI1->PowDecim() << " " << aI2->PowDecim() << "\n";
    mVI1.push_back(aI1);
    mVImHS.push_back(aI2);
}

void cAppli_NewRechPH::AddBrin(cBrinPtRemark * aBr)
{
   mVecB.push_back(aBr);
}



cPtSc CreatePtSc(const Pt2dr & aP,double aSc)
{
    cPtSc aRes;
    aRes.Pt()    = aP;
    aRes.Scale() = aSc;
    return aRes;
}

double Som_K_ApK(int aN,double a)
{
    double aSom = 0;
    for (int aK=1 ; aK<= aN ; aK++)
        aSom +=  pow(a,aK) * aK;
    return aSom;
}

/*
void TestSom(int aN,double a)
{
    // Verifie Som k=1,N  { a^k * k}

     double aCheck = Som_K_ApK(aN,a);
     // double aFormul = ((1-pow(a,aN+1)
}
*/
/*
double Sigma2FromFactExp(double a);
double FactExpFromSigma2(double aS2);


void TestSigma2(double a)
{
   int aNb = 100 + 100/(1-a);
   double aSom=0;
   double aSomK2=0;
   for (int aK=-aNb ; aK<= aNb ; aK++)
   {
       double aPds = pow(a,ElAbs(aK));
       aSom += aPds;
       aSomK2 += aPds * ElSquare(aK);
   }

   double aSigmaExp =  aSomK2 / aSom;
   // double aSigmaTh = (2*a) /((1+a) * ElSquare(1-a));
   // double aSigmaTh = (2*a) /(ElSquare(1-a));
   double aSigmaTh = Sigma2FromFactExp(a);

   double aATh = FactExpFromSigma2(aSigmaTh);

   std::cout << "TestSigma2 " << aSigmaExp << " " << aSigmaTh/aSigmaExp - 1 << " aaa " << a << " " << aATh << "\n";
}
*/

double cAppli_NewRechPH::IndScalDecim(const double & aS ) const
{
   return log(aS/mEch0Decim) / log(2.0);
}


cAppli_NewRechPH::cAppli_NewRechPH(int argc,char ** argv,bool ModeTest) :
    mNbByOct     (5.0),
    mPowS        (pow(2.0,1/mNbByOct)),
    mEch0Decim   (pow(2.0,1.5)),
    mNbS         (30),
    mISF         (-1,1e10),
    mStepSR      (1.0),
    mRollNorm    (false),
    mNbSR2Use    (10),
    mDeltaSR     (1),
    mNbTetaIm    (16),
    mMulNbTetaInv (4),
    mNbTetaInv   (mNbTetaIm*mMulNbTetaInv),
    mS0          (1.0),
    mScaleStab   (1.0),
    mSeuilAC     (0.95),
    mSeuilCR     (0.6),
    mScaleCorr   (false),
    mW1          (0),
    mModeTest    (ModeTest),
    mDistMinMax  (3.0),
    mDoMin       (true),
    mDoMax       (true),
    mDoPly       (false),
    mPlyC        (0),
    mDeZoomLn2   (7),
    mHistLong    (1000,0),
    mHistN0      (1000,0),
    mExtSave     ("Std"),
    mBasic       (false),
    mAddModeSift (true),
    mAddModeTopo (true),
    mLapMS       (false),
    mTestDirac   (false),
    mSaveFileLapl (false),
    // mPropCtrsIm0 (0.1),
    mNbSpace           (0),
    mNbScaleSpace      (0),
    mNbScaleSpaceCstr  (0),
    mDistAttenContr    (100.0),
    mPropContrAbs      (0.3),
    mSzContrast        (2),
    mPropCalcContr     (0.05),
    mIm0               (1,1),
    mTIm0              (mIm0),
    mImContrast        (1,1),
    mTImContrast       (mImContrast)
{
   cSinCardApodInterpol1D * aSinC = new cSinCardApodInterpol1D(cSinCardApodInterpol1D::eTukeyApod,5.0,5.0,1e-4,false);
   mInterp = new cTabIM2D_FromIm2D<tElNewRechPH>(aSinC,1000,false);

/*
   TestSigma2(0.1);
   TestSigma2(0.5);
   TestSigma2(0.9);
   TestSigma2(0.95);
   TestSigma2(0.99);
   getchar();
*/

   double aSeuilPersist = 1.0;

   MMD_InitArgcArgv(argc,argv);
   ElInitArgMain
   (
         argc,argv,
         LArgMain()  << EAMC(mName, "Name Image",  eSAM_IsPatFile),
         LArgMain()   << EAM(mPowS, "PowS",true,"Scale Pow")
                      << EAM(mNbS,  "NbS",true,"Number of level")
                      << EAM(mS0,   "S0",true,"ScaleInit, Def=1")
                      << EAM(mDoPly, "DoPly",true,"Generate ply file, for didactic purpose")
                      << EAM(mBoxCalc, "Box",true,"Box for computation")
                      << EAM(mModeTest, "Test",true,"if true add W")
                      << EAM(aSeuilPersist, "SP",true,"Threshold persistance")
                      << EAM(mBasic, "Basic",true,"Basic")
                      << EAM(mAddModeSift, "Sift",true,"Add SIFT Mode")
                      << EAM(mAddModeTopo, "Topo",true,"Add Topo Mode")
                      << EAM(mLapMS, "LapMS",true,"MulScale in Laplacian, def=false")
                      << EAM(mTestDirac, "Dirac",true,"Test with dirac image")
                      << EAM(mSaveFileLapl, "SaveLapl",true,"Save Laplacian file, def=false")
                      << EAM(mScaleStab, "SS",true,"Scale of Stability")
                      << EAM(mExtSave, "Save",true,"Extension for save")
                      << EAM(mScaleCorr, "ScCor",true,"Scale by correl")
                      << EAM(mISF, "ISF",true,"Interval scale forced")
   );


   ElTimer  aChrono;

    mNbSR2Calc = mRollNorm ?  (mNbSR2Use * 2 - 1) : mNbSR2Use;
    mMaxLevR     = mNbS - (mNbSR2Calc-1) * mDeltaSR;


   if (! EAMIsInit(&mExtSave))
   {
      mExtSave  = mBasic ? "Basic" : "Std";
   }
   if (! EAMIsInit(&mNbS))
   {
      if (mBasic) 
          mNbS = 1;
   }
   mNbInOct = log(2) / log(mPowS);

   if (mDoPly)
   {
      mPlyC = new cPlyCloud;
   }

   mP0Calc = Pt2di(0,0);
   mP1Calc = mTestDirac ? Pt2di(1000,1000) : Pt2di(-1,-1);
   if (EAMIsInit(&mBoxCalc))
   {
       mP0Calc = mBoxCalc._p0;
       mP1Calc = mBoxCalc._p1;
   }
   // Create top scale
   cOneScaleImRechPH * aIm0 = cOneScaleImRechPH::FromFile(*this,mS0,mName,mP0Calc,mP1Calc);
   AddScale(aIm0,aIm0);

   // Create matr of link, will have do it much less memory consuming (tiling of list ?)
   mIm0         = mVI1.back()->Im();
   mTIm0        = tTImNRPH(mIm0);
   mSzIm = mIm0.sz();
   mImContrast  = tImNRPH(mSzIm.x,mSzIm.y);
   mTImContrast = tTImNRPH(mImContrast);
   ComputeContrast();

   mSzLn2 = (mSzIm + Pt2di(mDeZoomLn2,mDeZoomLn2)) / mDeZoomLn2;
   mBufLnk2  = std::vector<std::vector<tLPtBuf> >(mSzLn2.y,std::vector<tLPtBuf>(mSzLn2.x,tLPtBuf() ));
   

   double aScaleMax = mS0*pow(mPowS,mNbS);
   // Used for  nearest point researh
   mVoisLnk = SortedVoisinDisk(-1,aScaleMax+4,true);


   if (! EAMIsInit(&mDZPlyLay))
   {
      mDZPlyLay = ElMin(mSzIm.x,mSzIm.y)/ double(mNbS);
   }
   if (mModeTest)
   {
      mW1 = Video_Win::PtrWStd(mSzIm);
   }
   // mVI1.back()->Show(mW1);

   int aPowDecim =1;
   // bool DecimBegun = false;  
   // Variables utilisee pour traiter
   int  KLastDecim=-11000;

   double aTimeLecture = TimeAndReset(aChrono);
  
   for (int aK=0 ; aK<mNbS ; aK++)
   {
        // Init from low resol
        if (aK!=0)
        {
           double aScale = mS0*pow(mPowS,aK);
           double ISD = IndScalDecim(aScale);
           double ISDPrec = IndScalDecim(mVI1.back()->ScaleAbs());
           bool  Decim = (ISD>0) && (round_down(ISD)!=round_down(ISDPrec));
           if (Decim) 
           {
              aPowDecim *= 2;
              KLastDecim = aK;
           }
           bool ReplikScale = (aK<= (KLastDecim+2));
           if (0)
           {
               std::cout << " IIII " << IndScalDecim(aScale) << " " << ISDPrec
                         << " RRR " << round_down(ISD) << " " << round_down(ISDPrec) 
                         << " " << aPowDecim 
                         <<  (Decim  ? "*" : " ") 
                         <<  (ReplikScale ? "#" : " ")
                         << "\n";
           }
           
           cOneScaleImRechPH * aImLS =  cOneScaleImRechPH::FromScale(*this,*mVI1.back(),aScale,aPowDecim,Decim);
           cOneScaleImRechPH * aImHS = aImLS;
           if (ReplikScale)
           {
               cOneScaleImRechPH * aBackHS = mVImHS.back();
               aImHS = cOneScaleImRechPH::FromScale(*this,*mVImHS.back(),aScale,aBackHS->PowDecim(),false);
           }
           AddScale(aImLS,aImHS);
           

        }
        std::cout << "DONE SCALE " << aK << " on " << mNbS  << "\n";
   }
   double aTimeGauss = TimeAndReset(aChrono);
   cSetPCarac aSPC;
   if (mAddModeSift)
   {
       for (int aK=0 ; aK<mNbS-1 ; aK++)
       {
           // std::cout << "SFUFT MAKE DIFF " << aK << "\n";
           bool aDoDifLS =  mVI1[aK]->SameDecim(*mVI1[aK+1]);
           if (aDoDifLS) 
              mVI1[aK]->SiftMakeDif(mVI1[aK+1]);

           bool aDoDifHS =      mVImHS[aK]->SameDecim(*mVImHS[aK+1])  // Il faut que les Hs soient compatible
                           && ( (mVI1[aK]!= mVImHS[aK]) || (mVI1[aK+1]!= mVImHS[aK+1]) ); // et que cela apport qqch
           if (aDoDifHS)
              mVImHS[aK]->SiftMakeDif(mVImHS[aK+1]);

           ELISE_ASSERT(aDoDifLS||aDoDifHS,"None diff possible");
           // std::cout << "HHHhhh " << aDoDifLS << " " << aDoDifHS << "\n";
       }
       for (int aK=1 ; aK<mNbS-2 ; aK++)
       {
            // bool
            cOneScaleImRechPH *  aILsUp   = mVI1[aK]->GetSameDecimSiftMade(mVI1[aK-1],mVImHS[aK-1]);
            cOneScaleImRechPH *  aILsDown = mVI1[aK]->GetSameDecimSiftMade(mVI1[aK+1],mVImHS[aK+1]);

            if ((aILsUp!=0) && (aILsDown!=0))
            {
                mVI1[aK]->SiftMaxLoc(aILsUp,aILsDown,aSPC);
            }
            else
            {
                cOneScaleImRechPH *  aIHsUp   = mVImHS[aK]->GetSameDecimSiftMade(mVI1[aK-1],mVImHS[aK-1]);
                cOneScaleImRechPH *  aIHsDown = mVImHS[aK]->GetSameDecimSiftMade(mVI1[aK+1],mVImHS[aK+1]);
                if ((aIHsUp!=0) && (aIHsDown!=0))
                {
                    mVImHS[aK]->SiftMaxLoc(aIHsUp,aIHsDown,aSPC);
                }
                else
                {
                    ELISE_ASSERT(false,"SIFT, scale assertion failed");
                }
            }
       }
   }
   double aTimeSift = TimeAndReset(aChrono);
   if (mAddModeTopo)
   {
       for (int aK=0 ; aK<mNbS ; aK++)
       {
            // Compute point of scale
            mVI1[aK]->CalcPtsCarac(mBasic);
            mVI1[aK]->Show(mW1);
          
            // Links the point at different scale
            if (aK!=0)
            {
               mVI1[aK]->CreateLink(*(mVI1[aK-1]));
            }
            std::cout << "DONE CARAC " << aK << " on " << mNbS << "\n";
       }


       Clik();

       for (int aK=0 ; aK<mNbS ; aK++)
       {
           mVI1[aK]->Export(aSPC,mPlyC);
       }

       if (mPlyC)
       {
           mPlyC->PutFile("NewH.ply");
       }

   }
   else
   {
   }


   double aTimeTopo = TimeAndReset(aChrono);

   for (auto & aPt : aSPC.OnePCarac())
       aPt.OK() = true;

   std::vector<cOnePCarac> aNewL;
   for (auto & aPt : aSPC.OnePCarac())
   {
       if (aPt.OK())  
          ComputeContrastePt(aPt);
       // 
       if (aPt.OK())
       {
          mVI1[aPt.NivScale()]->ComputeDirAC(aPt);
       }
       //  ComputeDirAC(cBrinPtRemark &);

       // C'est la que la decimation est invalidee
       if (aPt.OK())
       {
          mVI1[aPt.NivScale()]->AffinePosition(aPt);
       }
       mVI1[aPt.NivScale()]->NbPOfLab(int(aPt.Kind())) ++;
  }

  double aTimeAffine = TimeAndReset(aChrono);


  for (auto & aPt : aSPC.OnePCarac())
  {
       if (aPt.OK())
       {
          CalvInvariantRot(aPt);
       }

       // Put in global coord
       aPt.Pt() =  aPt.Pt() + Pt2dr(mP0Calc);
       if (aPt.OK())
          aNewL.push_back(aPt);
  }
  aSPC.OnePCarac() = aNewL;

  double aTimeCalcImage = TimeAndReset(aChrono);


  // MakeFileXML(aSPC,NameFileNewPCarac(mName,true,mExtSave));
  MakeFileXML(aSPC,NameFileNewPCarac(mName,true,mExtSave));

  double aTimeXml = TimeAndReset(aChrono);


  std::cout << "================== Time ========= \n";
  std::cout << "   Lecture =" << aTimeLecture << "\n";
  std::cout << "   Gaussian=" << aTimeGauss << "\n";
  std::cout << "   Laplacia=" << aTimeSift << "\n";
  std::cout << "   Topologi=" << aTimeTopo << "\n";
  std::cout << "   Affinage=" << aTimeAffine << "\n";
  std::cout << "   CalcInva=" << aTimeCalcImage << "\n";
  std::cout << "   SauvXml =" << aTimeXml << "\n";
}

bool  cAppli_NewRechPH::ComputeContrastePt(cOnePCarac & aPt)
{
   
   Symb_FNum aFI0  (mIm0.in_proj());
   Symb_FNum aFPds (1.0);
   double aS0,aS1,aS2;
   ELISE_COPY
   (
       disc(aPt.Pt(),aPt.Scale()),
       Virgule(aFPds,aFPds*aFI0,aFPds*Square(aFI0)),
       Virgule(sigma(aS0),sigma(aS1),sigma(aS2))
   );

   aS1 /= aS0;
   aS2 /= aS0;
   aS2 -= ElSquare(aS1);
   aS2 = sqrt(ElMax(aS2,1e-10)) * (aS0 / (aS0-1.0));
   aPt.Contraste() = aS2;
   aPt.ContrasteRel() = aS2 / mTImContrast.getproj(round_ni(aPt.Pt()));
    
   aPt.OK() = aPt.ContrasteRel() > mSeuilCR;

   return aPt.OK();
}

cOneScaleImRechPH * cAppli_NewRechPH::GetImOfNiv(int aNiv)
{
   return mVI1.at(aNiv);
}



void cAppli_NewRechPH::ComputeContrast()
{
   Symb_FNum aFI0  (mIm0.in_proj());
   Symb_FNum aFPds (1.0);
   Symb_FNum aFSom (rect_som(Virgule(aFPds,aFPds*aFI0,aFPds*Square(aFI0)),mSzContrast));

   Symb_FNum aS0 (aFSom.v0());
   Symb_FNum aS1 (aFSom.v1()/aS0);
   Symb_FNum aS2 (Max(0.0,aFSom.v2()/aS0 -Square(aS1)));

   tImNRPH aImC0  (mSzIm.x,mSzIm.y);
   double aNbVois = ElSquare(1+2*mSzContrast);
   // compute variance of image
   ELISE_COPY
   (
        mIm0.all_pts(),
       // ect_max(mIm0.in_proj(),mSzContrast)-rect_min(mIm0.in_proj(),mSzContrast),
        sqrt(aS2) * (aNbVois/(aNbVois-1.0)),
        mImContrast.out() | aImC0.out()
   );
   

   // Calcul d'une valeur  moyenne robuste
   std::vector<double> aVC;
   int aStep = 2*mSzContrast+1;
   for (int aX0=0 ; aX0<mSzIm.x ; aX0+= aStep)
   {
      for (int aY0=0 ; aY0<mSzIm.y ; aY0+= aStep)
      {
          int aX1 = ElMin(aX0+aStep,mSzIm.x);
          int aY1 = ElMin(aY0+aStep,mSzIm.y);
          // Calcul de la moyenne par carre
          double aSom = 0.0;
          double aSomF = 0.0;
          for (int aX=aX0 ; aX<aX1; aX++)
          {
              for (int aY=aY0 ; aY<aY1; aY++)
              {
                  aSom++;
                  aSomF += mTImContrast.get(Pt2di(aX,aY));
              }
          }
          aVC.push_back(aSomF/aSom);
      }
   }
   double aV0 = KthValProp(aVC,mPropCalcContr);
   double aV1 = KthValProp(aVC,1.0-mPropCalcContr);

   double aSom=0.0;
   double aSomF=0.0;
   for (const auto & aV : aVC)
   {
       if ((aV>=aV0) && (aV<=aV1))
       {
           aSom  ++;
           aSomF += aV;
       }
   }
   double aMoy = aSomF / aSom;
   double aFact = 1.0-1.0/mDistAttenContr;
   FilterExp(mImContrast,aFact);
   Im2D_REAL4 aIP1(mSzIm.x,mSzIm.y,1.0);
   FilterExp(aIP1,aFact);
   ELISE_COPY(mImContrast.all_pts(),mImContrast.in()/aIP1.in(),mImContrast.out());

   ELISE_COPY
   (
      mImContrast.all_pts(),
      mImContrast.in()*(1-mPropContrAbs)+mPropContrAbs*aMoy,
      mImContrast.out()
   );
   std::cout << "MOY = " << aMoy << "\n";


   // Pb si lance en //, si necessaire le mettre en option
   if (0)
   {
      Tiff_Im::CreateFromIm(aImC0,"ImC0.tif");
      Tiff_Im::CreateFromIm(mImContrast,"ImSeuilContraste.tif");
      Tiff_Im::CreateFromFonc("ImCRatio.tif",mSzIm,aImC0.in()/Max(1e-10,mImContrast.in()),GenIm::real4);
   }
}


void cAppli_NewRechPH::AdaptScaleValide(cOnePCarac & aPC)
{
    // Pour le faire bien il faut suivre le point, on le fait direct
    // pour les point max-min topo, et pour sift on verra
    ELISE_ASSERT(false,"cAppli_NewRechPH::AdaptScaleValide"); 
/*
   aPC.ScaleNature() = aPC.Scale();

   if (ScaleIsValid(aPC.Scale())) 
      return;

   for (const auto & aIm :  mVI1)
   {
      if (OkNivLapl(aIm->Niv()))
      {
      }

   }
*/
}



bool cAppli_NewRechPH::OkNivStab(int aNiv)
{
   return mVI1.at(aNiv)->ScaleAbs() >= mScaleStab;
}

bool cAppli_NewRechPH::Inside(const Pt2di & aP) const
{
    return (aP.x>=0) && (aP.y>=0) && (aP.x<mSzIm.x) && (aP.y<mSzIm.y);
}


double  cAppli_NewRechPH::DistMinMax(bool Basic) const  
{
   if (Basic)
   {
       return  60;
   }
   return mDistMinMax;
}

///   ============== Buf2 ==================

bool   cAppli_NewRechPH::InsideBuf2(const Pt2di & aP)
{
    return (aP.x>=0) && (aP.y>=0) && (aP.x<mSzLn2.x) && (aP.y<mSzLn2.y);
}

Pt2di cAppli_NewRechPH::PIndex2(const Pt2di & aP) const
{
   return Pt2di((aP.x+mDeZoomLn2/2)/mDeZoomLn2,(aP.y+mDeZoomLn2/2)/mDeZoomLn2);
}

tLPtBuf & cAppli_NewRechPH::LPOfBuf2(const Pt2dr & aPG)
{
    Pt2di aPInd = PIndex2(round_ni(aPG));

    return mBufLnk2.at(aPInd.y).at(aPInd.x);
}

void cAppli_NewRechPH::ClearBuff2(const tPtrPtRemark & aP)
{
   LPOfBuf2(aP->Pt()).clear();
}


void cAppli_NewRechPH::AddBuf2(const tPtrPtRemark & aP)
{
    LPOfBuf2(aP->Pt()).push_front(aP);
}

tPtrPtRemark  cAppli_NewRechPH::NearestPoint2(const Pt2di & aPG,const double & aDistFull)
{
    Pt2di aPInd = PIndex2(aPG);
    int aDistRed = round_up(aDistFull/mDeZoomLn2+1e-5);
    tPtrPtRemark aResult = nullptr;
    double aDMin= ElSquare(aDistFull);

    int aX0 = ElMax(0,aPInd.x-aDistRed);
    int aX1 = ElMin(mSzLn2.x-1,aPInd.x+aDistRed);
    int aY0 = ElMax(0,aPInd.y-aDistRed);
    int aY1 = ElMin(mSzLn2.y-1,aPInd.y+aDistRed);

    for (int aX=aX0 ; aX<=aX1; aX++)
    {
        for (int aY=aY0 ; aY<=aY1; aY++)
        {
             for (auto & itP : mBufLnk2[aY][aX])
             {
                 double aDist = square_euclid(Pt2dr(aPG)-(itP)->Pt());
                 if (aDist<=aDMin)
                 {
                      aDMin = aDist;
                      aResult = itP;
                 }
             }
        }
    }

    return aResult;
}

/*
        tPtrPtRemark  NearestPoint2(const Pt2di &,const double & aDist);
        void AddBuf2(const tPtrPtRemark &);
        void ClearBuff2(const tPtrPtRemark &);
*/



const Pt2di & cAppli_NewRechPH::SzIm() const  {return mSzIm;}

double cAppli_NewRechPH::ScaleAbsOfNiv(const int & aNiv) const
{
   return mVI1.at(aNiv)->ScaleAbs();
}

bool cAppli_NewRechPH::OkNivLapl(int aNiv)
{
   return (aNiv < int(mVI1.size())-2) ;
}

double cAppli_NewRechPH::GetLapl(int aNiv,const Pt2di & aP,bool &Ok)
{
   Ok = false;
   if (! OkNivLapl(aNiv))
      return 0;
   double aV1 = mVI1.at(aNiv)->GetVal(aP,Ok);
   if (!Ok)  return 0;
   double aV2 = mVI1.at(aNiv+1)->GetVal(aP,Ok);
   if (!Ok)  return 0;
   return aV1 - aV2;
}




int Test_NewRechPH(int argc,char ** argv)
{
   cAppli_NewRechPH anAppli(argc,argv,false);


   return EXIT_SUCCESS;

}

/*
int Generate_ImagSift(int argc,char ** argv)
{
     Pt2di aSz(1000,1000);
     Im2D_REAL4 aIm(aSz.x,aSz.y);

     for (int aKx=0 ; aKx<10 ; aKx++)
     {
         for (int aKy=0 ; aKy<10 ; aKy++)
         {
             Pt2di aP0(aKx*100,aKy*100);
             Pt2di aP1((aKx+1)*100,(aKy+1)*100);
             Pt2dr aMil = Pt2dr(aP0+aP1) / 2.0;

             double aSigma = (0.5*aKx+1.5*aKy+1);
             double aSign = ((aKx+aKy) % 2)   ? 1 : -1;

             ELISE_COPY
             (
                  rectangle(aP0,aP1),
                  128 * (1+aSign * exp(-  ( Square(FX-aMil.x) + Square(FY-aMil.y)) / Square(aSigma))),
                  aIm.out()
             );

         }
     }
     Tiff_Im::CreateFromIm(aIm,"TestSift.tif");
}
*/


/*Footer-MicMac-eLiSe-25/06/2007

Ce logiciel est un programme informatique servant à la mise en
correspondances d'images pour la reconstruction du relief.

Ce logiciel est régi par la licence CeCILL-B soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL-B telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL-B, et que vous en avez accepté les
termes.
aooter-MicMac-eLiSe-25/06/2007*/
