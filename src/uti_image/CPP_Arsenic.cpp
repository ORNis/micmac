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

#include "../../include/StdAfx.h"
#include "../src/uti_image/Arsenic.h"
#include "hassan/reechantillonnage.h"
#include <algorithm>
#include <functional>
#include <numeric>
#include <math.h>

void Arsenic_Banniere()
{
    std::cout <<  "\n";
    std::cout <<  " **********************************\n";
    std::cout <<  " *     A-utomated                 *\n";
    std::cout <<  " *     R-adiometric               *\n";
    std::cout <<  " *     S-hift                     *\n";
    std::cout <<  " *     E-qualization and          *\n";
    std::cout <<  " *     N-ormalization for         *\n";
    std::cout <<  " *     I-nter-image               *\n";
    std::cout <<  " *     C-orrection                *\n";
    std::cout <<  " **********************************\n\n";
}

inline double Dist3d(const Pt3d<double> & aP1, const Pt3d<double> & aP2) //TODO: check for euclid() for Pt3d<Type>
{
    return std::sqrt(std::pow(aP1.x-aP2.x,2)+std::pow(aP1.y-aP2.y,2)+std::pow(aP1.z-aP2.z,2));
}

std::vector<ArsenicImage> LoadGrpImagesMMByP(const cMMByImNM * aMMIN, int aDs, string InVig)
{

    cInterfChantierNameManipulateur * aICNM = cInterfChantierNameManipulateur::BasicAlloc(aDir);
    const std::vector<std::string> * aVectIm = aICNM->Get(aPatIm);

    std::list<string> ListConvert, ListVig;
    std::vector<std::string> VectImSc, VectMasq;
    size_t nbIm = aVectIm->size();

    // If a vignette correction folder was provided
    string postfix = "";
    if(InVig != "")
    {
        std::string cmdVig = MMDir() + "bin/mm3d Vodka \"" + aPatIm + "\" DoCor=1 Out=" + InVig + " InCal=" + InVig;
        postfix = "_Vodka.tif";
        ListVig.push_back(cmdVig);
        cEl_GPAO::DoComInParal(ListVig, aDir + "MkVig");
    }


    // Read images and Masks
    for (size_t aK1=0 ; aK1<nbIm ; ++aK1)
    {
        string cmdConv = MMDir() + "bin/ScaleIm " + InVig + (*aVectIm)[aK1] + postfix + " F8B=1 Out=" + (*aVectIm)[aK1] + "_Scaled.tif";
        ListConvert.push_back(cmdConv);

        //TODO: VectMasq.push_back("MM-Malt-Img-" + StdPrefix((aVectIm)[aK1]) + "/AutoMask_STD-MALT_Num_" + aEtape + ".tif");
        VectImSc.push_back((*aVectIm)[aK1] + std::string("_Scaled.tif"));
    }

    cEl_GPAO::DoComInParal(ListConvert, aDir + "MkScale");

    //Read the data
    std::vector<ArsenicImage> aGrIm;

    for (size_t aK1=0 ; aK1<nbIm ; ++aK1)
    {
        ArsenicImage aIm;
        //reading 3D info
        //TODO cElNuage3DMaille * info3D1 = cElNuage3DMaille::FromFileIm("MM-Malt-Img-" + StdPrefix((*aVectIm)[aK1]) + "/NuageImProf_STD-MALT_Etape_" + aEtape + ".xml");

        //TODO: aIm.info3D=info3D1;

        Tiff_Im aTF1= Tiff_Im::StdConvGen(aDir + VectImSc[aK1],3,false);
        Tiff_Im aTFM= Tiff_Im::StdConvGen(aDir + VectMasq[aK1],1,false);
        Pt2di aSz = aTF1.sz();
        Im2D_REAL4  aIm1R(aSz.x,aSz.y);
        Im2D_REAL4  aIm1G(aSz.x,aSz.y);
        Im2D_REAL4  aIm1B(aSz.x,aSz.y);
        Im2D_INT1  aMasq(aSz.x,aSz.y);
        ELISE_COPY
            (
                aTF1.all_pts(),
                aTF1.in(),
                Virgule(aIm1R.out(),aIm1G.out(),aIm1B.out())
            );

        ELISE_COPY
            (
                aTFM.all_pts(),
                aTFM.in(),
                aMasq.out()
            );

        aIm.Mask=aMasq;
        aIm.RChan=aIm1R;
        aIm.GChan=aIm1G;
        aIm.BChan=aIm1B;
        aIm.SZ=aSz;
        aGrIm.push_back(aIm);
    }

    return aGrIm;
}


cl_MatPtsHom ReadPtsHom3D(std::vector<ArsenicImage> & aGrIm,
                          const std::vector<std::pair<size_t, size_t>> & somePairs,
                          const int & ResolModel, const double & TPA)
{

    cl_MatPtsHom aMatPtsHomol;
    int nbIm = static_cast<int>(aGrIm.size());
    //going through each pair of different images
    for (int aK1=0 ; aK1<nbIm; ++aK1) //TODO: use somePairs there
    {
        //Creating the tie points vector
        cl_PtsRadio aPtsRadio;
        aMatPtsHomol.aMat.push_back(aPtsRadio);
        std::cout << "Extracting point from image "<< aK1+1<<" out of "<< nbIm << endl;
        //going through each point of image
        for (int aY=0 ; aY<aGrIm[aK1].SZ.y  ; aY++)
        {
            for (int aX=0 ; aX<aGrIm[aK1].SZ.x  ; aX++)
            {
                Pt2dr pos2DPtIm1;pos2DPtIm1.x=aX;pos2DPtIm1.y=aY;
                if(aGrIm[aK1].Mask.data()[aY][aX]==0){continue;}
                //If pts in masq, go look for 3D position
                    Pt3d<double> pos3DPtIm1=aGrIm[aK1].info3D->PreciseCapteur2Terrain(pos2DPtIm1);
                        //Testing the position of the point in other images
                        vector<double> distances(nbIm,10000000); //distances between original 3D point and reprojection from other images (init to "a lot")
                        vector<Pt2dr> pos2DOtherIm(nbIm);
                        for (int aK2=0 ; aK2<nbIm ; aK2++)
                        {
                            if (aK1!=aK2)// && aGrIm[aK2].info3D->PIsVisibleInImage(pos3DPtIm1))
                            {
                            Pt2dr pos2DPtIm2=aGrIm[aK2].info3D->Ter2Capteur(pos3DPtIm1);
                            //if pt in image and in masq, go look for 2D position, then 3D position
                            if(pos2DPtIm2.x>0 && pos2DPtIm2.x<aGrIm[aK2].SZ.x && pos2DPtIm2.y>0 && pos2DPtIm2.y<aGrIm[aK2].SZ.y){
                                        if(aGrIm[aK2].Mask.data()[int(pos2DPtIm2.y)][int(pos2DPtIm2.x)]){
                                            pos2DOtherIm[aK2]=pos2DPtIm2;
                                            Pt3d<double> pos3DPtIm2=aGrIm[aK2].info3D->PreciseCapteur2Terrain(pos2DPtIm2);
                                            //Compute Distance between the 2 3D points to check if they are the same ones (occlusion, beware!)
                                            distances[aK2]=Dist3d(pos3DPtIm1,pos3DPtIm2);
                                        }
                                    }
                            }else{aMatPtsHomol.aMat[aK1].SZ=aGrIm[aK1].SZ;}
                        }
                        for (int aK2=0 ; aK2<int(nbIm) ; aK2++){
                            if(distances[aK2]<(aGrIm[aK1].info3D->ResolSolOfPt(pos3DPtIm1))/TPA){//id pos3DPtIm1~=pos3DPtIm2 -->pt is considered homologous,it is added to PtsHom (Gr1, R1, G1, B1, X1, Y1, idem 2, NbPtsCouple++)
                                //Go looking for grey value of the point for each chan (Reechantillonnage/interpolation because pos2DPtIm not always integer)
                                //double Red1   =Reechantillonnage::biline(aGrIm[aK1].RChan.data(), aGrIm[aK1].SZ.x, aGrIm[aK1].SZ.y, pos2DPtIm1);
                                //double Green1 =Reechantillonnage::biline(aGrIm[aK1].GChan.data(), aGrIm[aK1].SZ.x, aGrIm[aK1].SZ.y, pos2DPtIm1);
                                //double Blue1  =Reechantillonnage::biline(aGrIm[aK1].BChan.data(), aGrIm[aK1].SZ.x, aGrIm[aK1].SZ.y, pos2DPtIm1);
                                double Red1   = aGrIm[aK1].RChan.data()[int(pos2DPtIm1.y)][int(pos2DPtIm1.x)];
                                double Green1 = aGrIm[aK1].GChan.data()[int(pos2DPtIm1.y)][int(pos2DPtIm1.x)];
                                double Blue1  = aGrIm[aK1].BChan.data()[int(pos2DPtIm1.y)][int(pos2DPtIm1.x)];
                                double Red2   = Reechantillonnage::biline(aGrIm[aK2].RChan.data(), aGrIm[aK2].SZ.x, aGrIm[aK2].SZ.y, pos2DOtherIm[aK2]);
                                double Green2 = Reechantillonnage::biline(aGrIm[aK2].GChan.data(), aGrIm[aK2].SZ.x, aGrIm[aK2].SZ.y, pos2DOtherIm[aK2]);
                                double Blue2  = Reechantillonnage::biline(aGrIm[aK2].BChan.data(), aGrIm[aK2].SZ.x, aGrIm[aK2].SZ.y, pos2DOtherIm[aK2]);
                                aMatPtsHomol.aMat[aK1].Pts.push_back(pos2DPtIm1.mul(ResolModel));
                                if(Red1>0){aMatPtsHomol.aMat[aK1].kR.push_back((1 + Red2/Red1 )/2);}else{aMatPtsHomol.aMat[aK1].kR.push_back((1 + Red2)/2);}
                                if(Green1>0){aMatPtsHomol.aMat[aK1].kG.push_back((1 + Green2/Green1 )/2);}else{aMatPtsHomol.aMat[aK1].kG.push_back((1 + Green2)/2);}
                                if(Blue1>0){aMatPtsHomol.aMat[aK1].kB.push_back((1 + Blue2/Blue1 )/2);}else{aMatPtsHomol.aMat[aK1].kB.push_back((1 + Blue2)/2);}
                                aMatPtsHomol.aMat[aK1].OtherIm.push_back(aK2);
                                aMatPtsHomol.aMat[aK1].SZ=aGrIm[aK1].SZ;
                                //file_out <<(Red2+Blue2+Green2)/(Red1+Blue1+Green1)<<endl;
                            }
                        }

            }
        }
    }
    return aMatPtsHomol;

}

cl_MatPtsHom TiePtsFilter(cl_MatPtsHom aMatPtsHomol, double aThresh)
{
    //Computing the mean of the correction factor
    int nbIm = (int)aMatPtsHomol.aMat.size();
    for (int numIm=0 ; numIm<nbIm ; numIm++)
    {
        std::cout<< "[Im "<< numIm <<"] NbPts before filter : "<< aMatPtsHomol.aMat[numIm].size()<<endl;
        double meanR=0,meanG=0,meanB=0;
        for(int i=0 ; i<aMatPtsHomol.aMat[numIm].size() ; i++)
        {
            meanR += aMatPtsHomol.aMat[numIm].kR[i];
            meanG += aMatPtsHomol.aMat[numIm].kG[i];
            meanB += aMatPtsHomol.aMat[numIm].kB[i];
        }
        meanR = meanR/aMatPtsHomol.aMat[numIm].size();
        std::cout << "Mean R = " << meanR << std::endl;
        meanG = meanG/aMatPtsHomol.aMat[numIm].size();
        std::cout << "Mean G = " << meanG << std::endl;
        meanB = meanB/aMatPtsHomol.aMat[numIm].size();
        std::cout << "Mean B = " << meanB << std::endl;

        //If a factor is different by more than "seuil" from the mean, the point is considered an outlier
        for(int i=aMatPtsHomol.aMat[numIm].size()-1 ; i>=0 ; i--)
        {
            if(aMatPtsHomol.aMat[numIm].kR[i]>meanR*aThresh || aMatPtsHomol.aMat[numIm].kR[i]<meanR/aThresh ||
               aMatPtsHomol.aMat[numIm].kG[i]>meanG*aThresh || aMatPtsHomol.aMat[numIm].kG[i]<meanG/aThresh ||
               aMatPtsHomol.aMat[numIm].kB[i]>meanB*aThresh || aMatPtsHomol.aMat[numIm].kB[i]<meanB/aThresh)
            {
                aMatPtsHomol.aMat[numIm].kR.erase(aMatPtsHomol.aMat[numIm].kR.begin() + i);
                aMatPtsHomol.aMat[numIm].kG.erase(aMatPtsHomol.aMat[numIm].kG.begin() + i);
                aMatPtsHomol.aMat[numIm].kB.erase(aMatPtsHomol.aMat[numIm].kB.begin() + i);
                aMatPtsHomol.aMat[numIm].Pts.erase(aMatPtsHomol.aMat[numIm].Pts.begin() + i);
                aMatPtsHomol.aMat[numIm].OtherIm.erase(aMatPtsHomol.aMat[numIm].OtherIm.begin() + i);
            }
        }
        std::cout << "[Im " << numIm << "] NbPts after filter: " << aMatPtsHomol.aMat[numIm].size() << std::endl;
    }
    return aMatPtsHomol;
}

void Egal_field_correct_ite(string aDir,std::vector<std::string> * aSetIm, cl_MatPtsHom aMatPtsHomol , string aDirOut, string InVig, int ResolModel, int nbIm, int nbIte, double aThresh)
{
for(int iter=0;iter<nbIte;iter++){
    std::cout << "Pass "<<iter+1<<" out of "<< nbIte<<endl;

    //Filtering the tie points
    aMatPtsHomol = TiePtsFilter(aMatPtsHomol, aThresh);

//Correcting the tie points

//#pragma omp parallel for

    for(int numImage1=0;numImage1<nbIm;numImage1++)
    {
        vector<int> cpt(nbIm,0);
        std::cout << "Computing factors for Im " << numImage1 << std::endl;

        //For each tie point point, compute correction value (distance-ponderated mean value of all the tie points)
        for(int k = 0; k<aMatPtsHomol.aMat[numImage1].size() ; k++){//go through each tie point
            double aCorR=0.0,aCorG=0.0,aCorB=0.0;
            double aSumDist=0;
            Pt2dr aPt(aMatPtsHomol.aMat[numImage1].Pts[k].x/ResolModel,aMatPtsHomol.aMat[numImage1].Pts[k].x/ResolModel);
            for(int numPt = 0; numPt<int(aMatPtsHomol.aMat[numImage1].size()) ; numPt++){//go through each tie point
                Pt2dr aPtIn(aMatPtsHomol.aMat[numImage1].Pts[numPt].x/ResolModel,aMatPtsHomol.aMat[numImage1].Pts[numPt].y/ResolModel);
                double aDist=euclid(aPtIn, aPt);
                if(aDist<1){aDist=1;}
                aSumDist=aSumDist+1/(aDist);
                aCorR = aCorR + aMatPtsHomol.aMat[numImage1].kR[numPt]/(aDist);
                aCorG = aCorG + aMatPtsHomol.aMat[numImage1].kG[numPt]/(aDist);
                aCorB = aCorB + aMatPtsHomol.aMat[numImage1].kB[numPt]/(aDist);
            }
            //Normalize
            aCorR = aCorR/aSumDist;
            aCorG = aCorG/aSumDist;
            aCorB = aCorB/aSumDist;

            //correcting Tie points color with computed surface
            aMatPtsHomol.aMat[numImage1].kR[k]=aCorR;
            aMatPtsHomol.aMat[numImage1].kG[k]=aCorG;
            aMatPtsHomol.aMat[numImage1].kB[k]=aCorB;
            }
    }
}

//Filtering the tie points
aMatPtsHomol = TiePtsFilter(aMatPtsHomol, aThresh);

cout<<"Factors were computed"<<endl;
//end truc à iterer--------------------------------------------------------------------------------------------------------------------------------------



    //Applying the correction to the images
    //Bulding the output file system
    ELISE_fp::MkDirRec(aDir + aDirOut);
    //Reading input files
    string suffix="";if(InVig!=""){suffix="_Vodka.tif";}


#ifdef USE_OPEN_MP
#pragma omp parallel for
#endif
    for(int i=0;i<nbIm;i++)
    {
        string aNameIm=InVig + (*aSetIm)[i] + suffix;//if vignette is used, change the name of input file to read
        std::cout<<"Correcting " << aNameIm <<" (with "<< aMatPtsHomol.aMat[i].size()<<" data points)"<<endl;
        string aNameOut=aDir + aDirOut + (*aSetIm)[i] +"_egal.tif";

        Pt2di aSzMod=aMatPtsHomol.aMat[i].SZ;//Size of the correction surface, taken from the size of the scaled image
        Im2D_REAL4  aImCorR(aSzMod.x,aSzMod.y,0.0);
        Im2D_REAL4  aImCorG(aSzMod.x,aSzMod.y,0.0);
        Im2D_REAL4  aImCorB(aSzMod.x,aSzMod.y,0.0);
        REAL4 ** aCorR = aImCorR.data();
        REAL4 ** aCorG = aImCorG.data();
        REAL4 ** aCorB = aImCorB.data();

        // For each point of the surface, compute correction value (distance - weighted mean value of all the tie points)
        long start=time(NULL);
        for (int aY=0 ; aY<aSzMod.y  ; aY++)
            {
                for (int aX=0 ; aX<aSzMod.x  ; aX++)
                {
                    double aCorPtR=0,aCorPtG=0,aCorPtB=0;
                    double aSumDist=0;
                    Pt2dr aPt(aX,aY);
                    for(int j = 0; j<aMatPtsHomol.aMat[i].size() ; j++){//go through each tie point
                        Pt2dr aPtIn(aMatPtsHomol.aMat[i].Pts[j].x/ResolModel,aMatPtsHomol.aMat[i].Pts[j].y/ResolModel);
                        float aDist = static_cast<float>(euclid(aPtIn, aPt));
                        if(aDist<1){aDist=1;}
                        aSumDist = aSumDist+1/(aDist);
                        aCorPtR = aCorPtR + aMatPtsHomol.aMat[i].kR[j]/(aDist);
                        aCorPtG = aCorPtG + aMatPtsHomol.aMat[i].kG[j]/(aDist);
                        aCorPtB = aCorPtB + aMatPtsHomol.aMat[i].kB[j]/(aDist);
                    }
                    //Normalize
                    aCorR[aY][aX] = static_cast<float>(aCorPtR/aSumDist);
                    aCorG[aY][aX] = static_cast<float>(aCorPtG/aSumDist);
                    aCorB[aY][aX] = static_cast<float>(aCorPtB/aSumDist);
                }
            }

        long end = time(NULL);
        cout<<"Correction field computed in "<<end-start<<" sec, applying..."<<endl;

        //Reading the image and creating the objects to be manipulated
        Tiff_Im aTF= Tiff_Im::StdConvGen(aDir + aNameIm,3,false);
        Pt2di aSz = aTF.sz();

        Im2D_U_INT1  aImR(aSz.x,aSz.y);
        Im2D_U_INT1  aImG(aSz.x,aSz.y);
        Im2D_U_INT1  aImB(aSz.x,aSz.y);

        ELISE_COPY
        (
           aTF.all_pts(),
           aTF.in(),
           Virgule(aImR.out(),aImG.out(),aImB.out())
        );

        U_INT1 ** aDataR = aImR.data();
        U_INT1 ** aDataG = aImG.data();
        U_INT1 ** aDataB = aImB.data();

        for (int aY=0 ; aY<aSz.y  ; aY++)
            {
                for (int aX=0 ; aX<aSz.x  ; aX++)
                {
                    Pt2dr aPt(double(aX/ResolModel),double(aY/ResolModel));
                    //To be able to correct the edges
                    if(aPt.x>aSzMod.x-2){aPt.x=aSzMod.x-2;}
                    if(aPt.y>aSzMod.y-2){aPt.y=aSzMod.y-2;}
                    //Bilinear interpolation from the scaled surface to the full scale image
                    double R = aDataR[aY][aX]*Reechantillonnage::biline(aCorR, aSzMod.x, aSzMod.y, aPt);
                    double G = aDataG[aY][aX]*Reechantillonnage::biline(aCorG, aSzMod.x, aSzMod.y, aPt);
                    double B = aDataB[aY][aX]*Reechantillonnage::biline(aCorB, aSzMod.x, aSzMod.y, aPt);
                    //Underflow and Overflow handling:
                    if(R>255){aDataR[aY][aX]=255;}else if(R<0){aDataR[aY][aX]=0;}else{aDataR[aY][aX]=R;}
                    if(G>255){aDataG[aY][aX]=255;}else if(G<0){aDataG[aY][aX]=0;}else{aDataG[aY][aX]=G;}
                    if(B>255){aDataB[aY][aX]=255;}else if(B<0){aDataB[aY][aX]=0;}else{aDataB[aY][aX]=B;}
                }
        }


        //Writing ouput image
         Tiff_Im  aTOut
            (
                aNameOut.c_str(),
                aSz,
                GenIm::u_int1,
                Tiff_Im::No_Compr,
                Tiff_Im::RGB
            );


         ELISE_COPY
             (
                 aTOut.all_pts(),
                 Virgule(aImR.in(),aImG.in(),aImB.in()),
                 aTOut.out()
             );

    }
}

int Arsenic_main(int argc,char ** argv)
{

    std::string aFullPattern,aDirOut="Egal/",InVig="";
    int ResolModel=16;
    double TPA=16,aThresh=1.4;
    int nbIte=5;
        //Reading arguments
        ElInitArgMain
        (
            argc,argv,
            LArgMain()  << EAMC(aFullPattern,"Images Pattern", eSAM_IsPatFile),
            LArgMain()  << EAM(aDirOut,"Out",true,"Output folder (end with /) and/or prefix (end with another char)")
                        << EAM(InVig,"InVig",true,"Input vignette folder (for example : Vignette/ )", eSAM_IsDir)
                        << EAM(ResolModel,"ResolModel",true,"Resol of input model (Def=16)")
                        << EAM(TPA,"TPA",true,"Tie Point Accuracy (Higher is better, lower gives more points Def=16)")
                        << EAM(nbIte,"NbIte",true,"Number of iteration of the process (default=5)")
                        << EAM(aThresh,"ThreshDisp",true,"Disparity threshold between the tie points (Def=1.4 for 40%)")
        );

        if (!MMVisualMode)
        {
            std::string aDir,aPatIm;
            SplitDirAndFile(aDir,aPatIm,aFullPattern);

            cInterfChantierNameManipulateur * aICNM = cInterfChantierNameManipulateur::BasicAlloc(aDir);
            const std::vector<std::string> * aSetIm = aICNM->Get(aPatIm);

            std::vector<std::string> aVectIm=*aSetIm;
            int nbIm = (int)aVectIm.size();

            ELISE_ASSERT(nbIm>1,"Less than two images found with this pattern");

            //Computing homologous points
            std::cout << "Computing homologous points" << std::endl;
            cl_MatPtsHom aMatPtsHomol=ReadPtsHom3D(aDir, aPatIm, InVig, ResolModel, TPA);

            //Computing and applying the equalization surface
            std::cout << "Computing and applying the equalization surface" << std::endl;
            Egal_field_correct_ite(aDir, & aVectIm, aMatPtsHomol, aDirOut, InVig, ResolModel, nbIm, nbIte, aThresh);

            Arsenic_Banniere();
        }

        return 0;
}


cAppliCorrColor::cAppliCorrColor(int argc, char **argv):
        cAppliWithSetImage(argc-1,argv+1,0),
        mModeMMByP (eQuickMac)
{

    std::string aPat, anOri,
            aDirOut="Egal/", aVignDir="";
    std::string aModeMMByP;
    int aDs=1;
    double TPA=16,aThresh=1.4;
    int nbIte=5;

    ElInitArgMain
            (
                    argc,argv,
                    LArgMain()  << EAMC(aPat,"Full Directory (Dir+Pattern)",eSAM_IsPatFile)
                                << EAMC(anOri,"Orientation", eSAM_IsExistDirOri),
                    LArgMain()  << EAM(aDirOut,"Out",true,"Output folder (end with /) and/or prefix (end with another char)")
                                << EAM(aModeMMByP,"ModeMMByP",true,"Mode used in MMByP", eSAM_None,ListOfVal(eNbTypeMMByP,"e"))
                                << EAM(aVignDir,"InVig",true,"Input vignette folder (for example : Vignette/ )")
                                << EAM(aDs,"DownScale",true,"Down scale images reduce computation time (Def=1)")
                                << EAM(TPA,"TPA",true,"Tie Point Accuracy (Higher is better, lower gives more points Def=16)")
                                << EAM(nbIte,"NbIte",true,"Number of iterations of the process (default=5)")
                                << EAM(aThresh,"ThreshDisp",true,"color discrepancy threshold between tie points (Def=1.4 for 40%)")
            );

    if (EAMIsInit(&aModeMMByP))
    {
        bool aModeHelp;
        StdReadEnum(aModeHelp,mModeMMByP,aModeMMByP,eNbTypeMMByP);
    }

    cMMByImNM * aMMIN = cMMByImNM::ForGlobMerge(Dir(),aDs,aModeMMByP); //TODO: make this a class field

    std::cout << this->DirAndPatFileOfImSec() << std::endl;

    Arsenic_Banniere();
}


int PIMsCorrColor_main(int argc, char ** argv)
{
    cAppliCorrColor anAppli(argc, argv);
    return EXIT_SUCCESS;
}

/*Footer-MicMac-eLiSe-25/06/2007

Ce logiciel est un programme informatique servant � la mise en
correspondances d'images pour la reconstruction du relief.

Ce logiciel est r�gi par la licence CeCILL-B soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL-B telle que diffus�e par le CEA, le CNRS et l'INRIA
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe �
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement,
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�.

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez
pris connaissance de la licence CeCILL-B, et que vous en avez accept� les
termes.
Footer-MicMac-eLiSe-25/06/2007*/
