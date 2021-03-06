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


#ifndef _NewRechMATCH_IMAGE_H_
#define _NewRechMATCH_IMAGE_H_



ElSimilitude SimilRobustInit(const ElPackHomologue & aPackFull,double aPropRan);

class cIndexCodeBinaire;
class cAFM_Im;
class cAFM_Im_Master;
class cAFM_Im_Sec ;
class cAppli_FitsMatch1Im;

class cSetOPC;

//============================================================

class cIndexCodeBinaire
{
    public :
         cIndexCodeBinaire (const cCompCB &,bool Overlap);
         const std::vector<cCompileOPC *> & VectVois(const cCompileOPC & aPC);
         void Add(cSetOPC &,const cFitsOneLabel &);
    private :
         void Add(cCompileOPC & anOpc);
         int                      mNBBTot;
         int                      mNBBVois;
         const std::vector<int> * mFlagV;
         std::vector<std::vector<cCompileOPC *> > mVTabIndex;
         bool                                     mOverlap;
};

class cCdtCplHom
{
    public :
       cCdtCplHom(cCompileOPC * aPM,cCompileOPC * aPS,double aCorr,int aShift) :
           mPM    (aPM),
           mPS    (aPS),
           mCorr  (aCorr),
           mShift (aShift),
           mOk    (true),
           mDistS (0)
       {
       }

       Pt2dr PM() const {return mPM->mOPC.Pt();}
       Pt2dr PS() const {return mPS->mOPC.Pt();}
       cCompileOPC * mPM;
       cCompileOPC * mPS;
       double        mCorr;
       int           mShift;
       bool          mOk;
       double        mDistS;
};

class cPtFromPCC
{
   public :
       Pt2dr operator() (cCdtCplHom * aCC) { return aCC->PM(); }
};

typedef ElQT<cCdtCplHom*,Pt2dr,cPtFromPCC> tQtCC ;

class cSetOPC
{
    public :
       void FiltrageFromHighestScale(const cSetOPC & aSet,int aNb,bool Show);

       void InitLabel(const cFitsOneLabel &,const cSeuilFitsParam&,bool DoIndex,bool Overlap);
       cIndexCodeBinaire & Ind();
       const cFitsOneLabel &  FOL() const;
       const cSeuilFitsParam &  Seuil() const;
       cSetOPC();
       const std::vector<cCompileOPC> &  VOpc() const;
       std::vector<cCompileOPC> &  VOpc() ;
       void Add(const cCompileOPC&);
       cCompileOPC& At(int aK);
        

    private :
       std::vector<cCompileOPC>  mVOpc;
       cIndexCodeBinaire  *   mIndexCB;
       const cFitsOneLabel *  mFOL;
       const cSeuilFitsParam *  mSeuil;
};


class cAFM_Im
{
     public :
         friend class cAFM_Im_Master;
         friend class cAFM_Im_Sec;
         friend class cAppli_FitsMatch1Im;

         cAFM_Im (const std::string  &,cAppli_FitsMatch1Im &);
         ~cAFM_Im ();
         void SetFlagVSetCC(bool Index);

     protected :
         cAppli_FitsMatch1Im & mAppli;
         std::string mNameIm;
         cMetaDataPhoto mMTD;
         Pt2di       mSzIm;
         cSetPCarac                             mSetPC;
         std::vector<cSetOPC >                  mVSetCC;
         cSetOPC                                mSetInd0; // decision rapide sur l'overlap
         cSetOPC                                mSetInd1;
         // std::vector<cCompileOPC> &             mVIndex;

};


class cAFM_Im_Master : public  cAFM_Im
{
     public :
         cAFM_Im_Master (const std::string  &,cAppli_FitsMatch1Im &);
         void MatchOne(bool OverLap,cAFM_Im_Sec & , cSetOPC & ,cSetOPC & ,std::vector<cCdtCplHom> & ,int aNbMin);

         bool             MatchGlob(cAFM_Im_Sec &);

         void FiltrageSpatialGlob(std::vector<cCdtCplHom> & aVCpl,int aNbMin);


         // std::vector<std::vector<cCompileOPC *> > mVTabIndex;
         void FilterVoisCplCt(std::vector<cCdtCplHom> & aV);
         void RemoveCpleQdt(cCdtCplHom &);
         cPtFromPCC  mArgQt;
         tQtCC   mQt;
};


class cAFM_Im_Sec : public  cAFM_Im
{
     public :
         cAFM_Im_Sec (const std::string  &,cAppli_FitsMatch1Im &);
};

class cAppli_FitsMatch1Im
{
     public :
          cAppli_FitsMatch1Im(int argc,char ** argv);
          const std::string &   ExtNewH () const {return    mExtNewH;}
          const cFitsParam & FitsPm() const {return mFitsPm;}
          std::string NameCple(const std::string & aN1,const std::string & aN2) const;
          int NbBIndex() const;
          int ThreshBIndex() const;
          Pt2di  NbMaxS0() const;
          bool  ShowDet() const;
          eTypePtRemark    LabOL() const;

     private :
          cFitsParam         mFitsPm;
          std::string        mNameMaster;
          std::string        mPatIm;
          cElemAppliSetFile  mEASF;
          cAFM_Im_Master *   mImMast;
          cAFM_Im_Sec *      mCurImSec;
          std::string        mNameXmlFits;
          std::string        mExtNewH;
          std::string        mSH;
          std::string        mPostHom;
          bool               mExpTxt;
          int                mNbBIndex;
          int                mThreshBIndex;
          bool               mOneWay;
          bool               mShowDet;
          bool               mCallBack;
          Pt2di              mNbMaxS0;  // Nb max en presel x=> pour overlap en point a analyser, y=> pour modele 3D, y en point voulu
          eTypePtRemark      mLabOL;
};






#endif //  _NewRechMATCH_IMAGE_H_


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
