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

void Arsenic_Banniere() {
    std::cout << "\n";
    std::cout << " **********************************\n";
    std::cout << " *     A-utomated                 *\n";
    std::cout << " *     R-adiometric               *\n";
    std::cout << " *     S-hift                     *\n";
    std::cout << " *     E-qualization and          *\n";
    std::cout << " *     N-ormalization for         *\n";
    std::cout << " *     I-nter-image               *\n";
    std::cout << " *     C-orrection                *\n";
    std::cout << " **********************************\n\n";
}

inline double Dist3d(const Pt3d<double> &aP1, const Pt3d<double> &aP2) //TODO: check for euclid() for Pt3d<Type>
{
    return std::sqrt(std::pow(aP1.x - aP2.x, 2) + std::pow(aP1.y - aP2.y, 2) + std::pow(aP1.z - aP2.z, 2));
}

void TiePtsFilter(cl_MatPtsHom &aMatPtsHomol, double aThresh) {
    //Computing the mean of the correction factor
    int nbIm = (int) aMatPtsHomol.aMat.size();
    for (int numIm = 0; numIm < nbIm; numIm++) {
        std::cout << "[Im " << numIm << "] NbPts before filter : " << aMatPtsHomol.aMat[numIm].size() << std::endl;
        double meanR = 0, meanG = 0, meanB = 0;
        for (int i = 0; i < aMatPtsHomol.aMat[numIm].size(); i++) {
            meanR += aMatPtsHomol.aMat[numIm].kR[i];
            meanG += aMatPtsHomol.aMat[numIm].kG[i];
            meanB += aMatPtsHomol.aMat[numIm].kB[i];
        }
        meanR = meanR / aMatPtsHomol.aMat[numIm].size();
        std::cout << "Mean R = " << meanR << std::endl;
        meanG = meanG / aMatPtsHomol.aMat[numIm].size();
        std::cout << "Mean G = " << meanG << std::endl;
        meanB = meanB / aMatPtsHomol.aMat[numIm].size();
        std::cout << "Mean B = " << meanB << std::endl;

        //If a factor is different by more than "seuil" from the mean, the point is considered an outlier
        for (int i = aMatPtsHomol.aMat[numIm].size() - 1; i >= 0; i--) {
            if (aMatPtsHomol.aMat[numIm].kR[i] > meanR * aThresh || aMatPtsHomol.aMat[numIm].kR[i] < meanR / aThresh ||
                aMatPtsHomol.aMat[numIm].kG[i] > meanG * aThresh || aMatPtsHomol.aMat[numIm].kG[i] < meanG / aThresh ||
                aMatPtsHomol.aMat[numIm].kB[i] > meanB * aThresh || aMatPtsHomol.aMat[numIm].kB[i] < meanB / aThresh) {
                aMatPtsHomol.aMat[numIm].kR.erase(aMatPtsHomol.aMat[numIm].kR.begin() + i);
                aMatPtsHomol.aMat[numIm].kG.erase(aMatPtsHomol.aMat[numIm].kG.begin() + i);
                aMatPtsHomol.aMat[numIm].kB.erase(aMatPtsHomol.aMat[numIm].kB.begin() + i);
                aMatPtsHomol.aMat[numIm].Pts.erase(aMatPtsHomol.aMat[numIm].Pts.begin() + i);
                aMatPtsHomol.aMat[numIm].OtherIm.erase(aMatPtsHomol.aMat[numIm].OtherIm.begin() + i);
            }
        }
        std::cout << "[Im " << numIm << "] NbPts after filter: " << aMatPtsHomol.aMat[numIm].size() << std::endl;
    }
}

void Egal_field_correct_ite(
        std::string aDir,
        std::vector<std::string> *aSetIm,
        cl_MatPtsHom &aMatPtsHomol,
        const std::string &aDirOut,
        const std::string InVig,
        int ResolModel,
        int nbIm,
        int nbIte,
        double aThresh) {
    for (int iter = 0; iter < nbIte; iter++) {
        std::cout << "Pass " << iter + 1 << " out of " << nbIte << std::endl;
        //Filtering the tie points
        TiePtsFilter(aMatPtsHomol, aThresh);

        //Correcting the tie points
        //TODO: #pragma omp parallel for
        for (int numImage1 = 0; numImage1 < nbIm; ++numImage1) {
            std::cout << "Computing factors for Im " << numImage1 << std::endl;

            // For each tie point point, compute correction value (distance-ponderated mean value of all the tie points)
            for (int k = 0; k < aMatPtsHomol.aMat[numImage1].size(); k++) {//go through each tie point
                double aCorR = 0.0, aCorG = 0.0, aCorB = 0.0;
                double aSumDist = 0;
                Pt2dr aPt(aMatPtsHomol.aMat[numImage1].Pts[k].x / ResolModel,
                          aMatPtsHomol.aMat[numImage1].Pts[k].x / ResolModel);
                for (int numPt = 0;
                     numPt < int(aMatPtsHomol.aMat[numImage1].size()); numPt++) {//go through each tie point
                    Pt2dr aPtIn(aMatPtsHomol.aMat[numImage1].Pts[numPt].x / ResolModel,
                                aMatPtsHomol.aMat[numImage1].Pts[numPt].y / ResolModel);
                    double aDist = euclid(aPtIn, aPt);
                    if (aDist < 1) { aDist = 1; }
                    aSumDist = aSumDist + 1 / (aDist);
                    aCorR = aCorR + aMatPtsHomol.aMat[numImage1].kR[numPt] / (aDist);
                    aCorG = aCorG + aMatPtsHomol.aMat[numImage1].kG[numPt] / (aDist);
                    aCorB = aCorB + aMatPtsHomol.aMat[numImage1].kB[numPt] / (aDist);
                }
                //Normalize
                aCorR = aCorR / aSumDist;
                aCorG = aCorG / aSumDist;
                aCorB = aCorB / aSumDist;

                //correcting Tie points color with computed surface
                aMatPtsHomol.aMat[numImage1].kR[k] = aCorR;
                aMatPtsHomol.aMat[numImage1].kG[k] = aCorG;
                aMatPtsHomol.aMat[numImage1].kB[k] = aCorB;
            }
        }
    }

    // Filter the tie points
    TiePtsFilter(aMatPtsHomol, aThresh);

    std::cout << "Factors were computed" << std::endl;

    // Applying the correction to the images
    // Bulding the output file system
    ELISE_fp::MkDirRec(aDir + aDirOut);
    //Reading input files
    string suffix = "";
    if (InVig != "") { suffix = "_Vodka.tif"; }


#ifdef USE_OPEN_MP
#pragma omp parallel for
#endif
    for (int i = 0; i < nbIm; i++) {
        string aNameIm = InVig + (*aSetIm)[i] + suffix;//if vignette is used, change the name of input file to read
        std::cout << "Correcting " << aNameIm << " (with " << aMatPtsHomol.aMat[i].size() << " data points)" << endl;
        string aNameOut = aDir + aDirOut + (*aSetIm)[i] + "_egal.tif";

        Pt2di aSzMod = aMatPtsHomol.aMat[i].SZ;//Size of the correction surface, taken from the size of the scaled image
        Im2D_REAL4 aImCorR(aSzMod.x, aSzMod.y, 0.0);
        Im2D_REAL4 aImCorG(aSzMod.x, aSzMod.y, 0.0);
        Im2D_REAL4 aImCorB(aSzMod.x, aSzMod.y, 0.0);
        REAL4 **aCorR = aImCorR.data();
        REAL4 **aCorG = aImCorG.data();
        REAL4 **aCorB = aImCorB.data();

        // For each point of the surface, compute correction value (distance - weighted mean value of all the tie points)
        long start = time(NULL);
        for (int aY = 0; aY < aSzMod.y; aY++) {
            for (int aX = 0; aX < aSzMod.x; aX++) {
                double aCorPtR = 0, aCorPtG = 0, aCorPtB = 0;
                double aSumDist = 0;
                Pt2dr aPt(aX, aY);
                for (int j = 0; j < aMatPtsHomol.aMat[i].size(); j++) {//go through each tie point
                    Pt2dr aPtIn(aMatPtsHomol.aMat[i].Pts[j].x / ResolModel, aMatPtsHomol.aMat[i].Pts[j].y / ResolModel);
                    float aDist = static_cast<float>(euclid(aPtIn, aPt));
                    if (aDist < 1) { aDist = 1; }
                    aSumDist = aSumDist + 1 / (aDist);
                    aCorPtR = aCorPtR + aMatPtsHomol.aMat[i].kR[j] / (aDist);
                    aCorPtG = aCorPtG + aMatPtsHomol.aMat[i].kG[j] / (aDist);
                    aCorPtB = aCorPtB + aMatPtsHomol.aMat[i].kB[j] / (aDist);
                }
                //Normalize
                aCorR[aY][aX] = static_cast<float>(aCorPtR / aSumDist);
                aCorG[aY][aX] = static_cast<float>(aCorPtG / aSumDist);
                aCorB[aY][aX] = static_cast<float>(aCorPtB / aSumDist);
            }
        }

        long end = time(NULL);
        cout << "Correction field computed in " << end - start << " sec, applying..." << std::endl;

        //Reading the image and creating the objects to be manipulated
        Tiff_Im aTF = Tiff_Im::StdConvGen(aDir + aNameIm, 3, false);
        Pt2di aSz = aTF.sz();

        Im2D_U_INT1 aImR(aSz.x, aSz.y);
        Im2D_U_INT1 aImG(aSz.x, aSz.y);
        Im2D_U_INT1 aImB(aSz.x, aSz.y);

        ELISE_COPY
        (
                aTF.all_pts(),
                aTF.in(),
                Virgule(aImR.out(), aImG.out(), aImB.out())
        );

        U_INT1 **aDataR = aImR.data();
        U_INT1 **aDataG = aImG.data();
        U_INT1 **aDataB = aImB.data();

        for (int aY = 0; aY < aSz.y; aY++) {
            for (int aX = 0; aX < aSz.x; aX++) {
                Pt2dr aPt(double(aX / ResolModel), double(aY / ResolModel));
                //To be able to correct the edges
                if (aPt.x > aSzMod.x - 2) { aPt.x = aSzMod.x - 2; }
                if (aPt.y > aSzMod.y - 2) { aPt.y = aSzMod.y - 2; }
                //Bilinear interpolation from the scaled surface to the full scale image
                double R = aDataR[aY][aX] * Reechantillonnage::biline(aCorR, aSzMod.x, aSzMod.y, aPt);
                double G = aDataG[aY][aX] * Reechantillonnage::biline(aCorG, aSzMod.x, aSzMod.y, aPt);
                double B = aDataB[aY][aX] * Reechantillonnage::biline(aCorB, aSzMod.x, aSzMod.y, aPt);
                //Underflow and Overflow handling:
                if (R > 255) { aDataR[aY][aX] = 255; } else if (R <
                                                                0) { aDataR[aY][aX] = 0; } else { aDataR[aY][aX] = R; }
                if (G > 255) { aDataG[aY][aX] = 255; } else if (G <
                                                                0) { aDataG[aY][aX] = 0; } else { aDataG[aY][aX] = G; }
                if (B > 255) { aDataB[aY][aX] = 255; } else if (B <
                                                                0) { aDataB[aY][aX] = 0; } else { aDataB[aY][aX] = B; }
            }
        }


        //Writing ouput image
        Tiff_Im aTOut
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
                Virgule(aImR.in(), aImG.in(), aImB.in()),
                aTOut.out()
        );

    }
}

int Arsenic_main(int argc, char **argv) { } //FIXME Axe that

cAppliCorrColor::cAppliCorrColor(int argc, char **argv) :
        cAppliWithSetImage(argc - 1, argv + 1, 0),
        mModeMMByP(eQuickMac) {

    std::string aPat, anOri,
            aDirOut = "Egal/", aVignDir = "";
    std::string aModeMMByP;
    int aDs = 1;
    double TPA = 16, aThresh = 1.4;
    int nbIte = 5;

    ElInitArgMain
            (
                    argc, argv,
                    LArgMain() << EAMC(aPat, "Full Directory (Dir+Pattern)", eSAM_IsPatFile)
                    << EAMC(anOri, "Orientation", eSAM_IsExistDirOri),
                    LArgMain() <<
                    EAM(aDirOut, "Out", true, "Output folder (end with /) and/or prefix (end with another char)")
                    << EAM(aModeMMByP, "ModeMMByP", true, "Mode used in MMByP", eSAM_None, ListOfVal(eNbTypeMMByP, "e"))
                    << EAM(aVignDir, "InVig", true, "Input vignette folder (for example : Vignette/ )")
                    << EAM(aDs, "DownScale", true, "Down scale images reduce computation time (Def=1)")
                    << EAM(TPA, "TPA", true, "Tie Point Accuracy (Higher is better, lower gives more points Def=16)")
                    << EAM(nbIte, "NbIte", true, "Number of iterations of the process (default=5)")
                    <<
                    EAM(aThresh, "ThreshDisp", true, "color discrepancy threshold between tie points (Def=1.4 for 40%)")
            );

    if (EAMIsInit(&aModeMMByP)) {
        bool aModeHelp;
        StdReadEnum(aModeHelp, mModeMMByP, aModeMMByP, eNbTypeMMByP);
    }

    //mMMIN = cMMByImNM::ForGlobMerge(Dir(),aDs,aModeMMByP); //FIXME aDs != MMByP aDs

    ELISE_fp::MkDirSvp(NameTmpDirArsenic());
    mMMIN = cMMByImNM::ForGlobMerge(Dir(), 1, aModeMMByP);
    LoadGrpImagesMMByP(aVignDir, aDs);
    //TODO: read clouds and discard image without enough points

    //Computing homologous points
    std::cout << "Computing homologous points" << std::endl;
    ReadPtsHom3D(TPA);

    /*
    //Computing and applying the equalization surface
    std::cout << "Computing and applying the equalization surface" << std::endl;
    Egal_field_correct_ite(aDir, & aVectIm, aMatPtsHomol, aDirOut, InVig, ResolModel, nbIm, nbIte, aThresh);
    */

    Arsenic_Banniere();
}


void cAppliCorrColor::LoadGrpImagesMMByP(const std::string &InVig, int aDs) //TODO: make these class fields
{

    cInterfChantierNameManipulateur *aICNM = cInterfChantierNameManipulateur::BasicAlloc(Dir());
    const std::vector<std::string> *aVectIm = aICNM->Get(mMMIN->KeyFileLON());

    std::list<string> ListConvert, ListVig;
    std::vector<std::string> VectImSc, VectCloud;
    size_t nbIm = aVectIm->size();

    // If a vignette correction folder was provided
    std::string postfix = "";
    if (InVig != "") {
        std::string cmdVig = MM3dBinFile_quotes("Vodka")
                             + mMMIN->KeyFileLON()
                             + " DoCor=1 Out=" + InVig
                             + " InCal=" + InVig;

        postfix = "_Vodka.tif";
        ListVig.push_back(cmdVig);
        cEl_GPAO::DoComInParal(ListVig);
    }

    //TODO: use a mechanism "a la" MergeCloud

    // Read images, clouds
    for (size_t aK1 = 0; aK1 < nbIm; ++aK1) {
        std::string aRelPathIm = InVig + "." + ELISE_STR_DIR + (*aVectIm)[aK1] + postfix;
        std::string aRelPathCloud = mMMIN->NameFileXml(eTMIN_Depth,
                                                       (*aVectIm)[aK1]); //FIXME names should be computed by a method
        //TODO: Glob dir ? ^
        if (aDs != 1) {
            std::string aNewRelPathIm = NameTmpDirArsenic() + ELISE_STR_DIR + (*aVectIm)[aK1] + ".tif";
            std::string aNewRelPathCloud = NameTmpDirArsenic() + ELISE_STR_DIR +
                                           cElFilename(mMMIN->NameFileXml(eTMIN_Depth, (*aVectIm)[aK1])).m_basename;

            std::string cmdConvImg = MM3dBinFile_quotes("ScaleIm")
                                     + aRelPathIm + BLANK
                                     + ToString(aDs + 2) // FIXME: DOSomething with that
                                     + " F8B=1 Out=" + aNewRelPathIm;

            std::string cmdConvCloud = MM3dBinFile_quotes("ScaleNuage")
                                       + aRelPathCloud + BLANK
                                       + (aNewRelPathCloud.substr(0, aNewRelPathCloud.size() -
                                                                     4)) //TODO: is there something better in Elise?
                                       + BLANK
                                       + ToString(aDs) + BLANK
                                       + "InDirLoc=0";

            ListConvert.push_back(cmdConvImg);
            ListConvert.push_back(cmdConvCloud);

            std::swap(aRelPathIm, aNewRelPathIm);
            std::swap(aRelPathCloud, aNewRelPathCloud);
        }

        VectImSc.push_back(aRelPathIm);
        VectCloud.push_back(aRelPathCloud);
    }

    cEl_GPAO::DoComInParal(ListConvert);

    for (size_t aK1 = 0; aK1 < nbIm; ++aK1) {
        ArsenicImage aIm; //TODO: a true constructor or a struct

        // reading 3D info
        aIm.info3D = cElNuage3DMaille::FromFileIm(Dir() + VectCloud[aK1]);
        Tiff_Im aTF1 = Tiff_Im::StdConvGen(Dir() + VectImSc[aK1], 3, false);
        Pt2di aSz = aTF1.sz();
        Im2D_REAL4 aIm1R(aSz.x, aSz.y);
        Im2D_REAL4 aIm1G(aSz.x, aSz.y);
        Im2D_REAL4 aIm1B(aSz.x, aSz.y);
        ELISE_COPY(aTF1.all_pts(),
                   aTF1.in(),
                   Virgule(aIm1R.out(), aIm1G.out(), aIm1B.out()));
        aIm.RChan = aIm1R;
        aIm.GChan = aIm1G;
        aIm.BChan = aIm1B;
        aIm.SZ = aSz;
        mVectArsenicImage.push_back(aIm);
    }
}


void cAppliCorrColor::ReadPtsHom3D(const double &TPA) {

    size_t nbIm = mVectArsenicImage.size();
    //going through each pair of different images
    for (int aK1 = 0; aK1 < nbIm; ++aK1) //TODO: use Some Pairs there or a graph precomputed
    {
        // Create the vector of tie points
        cl_PtsRadio aPtsRadio;
        mMatPtsHom.aMat.push_back(aPtsRadio);
        std::cout << "Extracting points from image " << aK1 + 1 << " out of " << nbIm << std::endl;
        ArsenicImage & anIm1 = mVectArsenicImage[aK1];

        //going through each point of image

        for (int aY = 0; aY < anIm1.SZ.y; ++aY) {
            for (int aX = 0; aX < anIm1.SZ.x; ++aX) {
                Pt2dr pos2DPtIm1(aX, aY);

                if (anIm1.info3D->ImMask().get(aX, aY) == 0)
                    continue;

                std::cout << aX << " " << aY << std::endl;
                // If pt is in masq, go look for 3D position
                Pt3dr pos3DPtIm1 = anIm1.info3D->PreciseCapteur2Terrain(pos2DPtIm1);

                // Testing the position of the point in other images
                std::vector<double> distances(nbIm, std::numeric_limits<double>::max()); //distances between original 3D point and reprojection from other images (init to "a lot")
                std::vector<Pt2dr> pos2DOtherIm(nbIm);

                for (int aK2 = 0; aK2 < nbIm; ++aK2) {
                    if (aK1 != aK2)// && aGrIm[aK2].info3D->PIsVisibleInImage(pos3DPtIm1))
                    {
                        ArsenicImage & anIm2 = mVectArsenicImage[aK2];

                        Pt2dr pos2DPtIm2 = anIm2.info3D->Ter2Capteur(pos3DPtIm1);
                        // if the pt is in image and in masq, go look for 2D position, then 3D position
                        std::cout << pos2DPtIm2.x << " " << pos2DPtIm2.y << std::endl;
                        if (pos2DPtIm2.x > 0 && pos2DPtIm2.x < anIm2.SZ.x && pos2DPtIm2.y > 0 && pos2DPtIm2.y < anIm2.SZ.y)
                        {
                            if (anIm2.info3D->ImMask().get(round_ni(pos2DPtIm2.x), round_ni(pos2DPtIm2.y)) !=0) {
                                pos2DOtherIm[aK2] = pos2DPtIm2;
                                Pt3dr pos3DPtIm2 = anIm2.info3D->PreciseCapteur2Terrain(pos2DPtIm2);
                                // Compute Distance between the 2 3D points to check if they are the same ones (occlusion, beware!)
                                distances[aK2] = Dist3d(pos3DPtIm1, pos3DPtIm2);
                            }
                        }
                    }
                    else {
                        mMatPtsHom.aMat[aK1].SZ = anIm1.SZ;
                    }
                }

                for (int aK2 = 0; aK2 < nbIm; aK2++) {
                    if (distances[aK2] < (mVectArsenicImage[aK1].info3D->ResolSolOfPt(pos3DPtIm1)) /
                                         TPA) {//id pos3DPtIm1~=pos3DPtIm2 -->pt is considered homologous,it is added to PtsHom (Gr1, R1, G1, B1, X1, Y1, idem 2, NbPtsCouple++)
                        //Go looking for grey value of the point for each chan (Reechantillonnage/interpolation because pos2DPtIm not always integer)
                        double Red1 = mVectArsenicImage[aK1].RChan.data()[int(pos2DPtIm1.y)][int(pos2DPtIm1.x)];
                        double Green1 = mVectArsenicImage[aK1].GChan.data()[int(pos2DPtIm1.y)][int(pos2DPtIm1.x)];
                        double Blue1 = mVectArsenicImage[aK1].BChan.data()[int(pos2DPtIm1.y)][int(pos2DPtIm1.x)];

                        double Red2 = Reechantillonnage::biline(mVectArsenicImage[aK2].RChan.data(),
                                                                mVectArsenicImage[aK2].SZ.x,
                                                                mVectArsenicImage[aK2].SZ.y, pos2DOtherIm[aK2]);
                        double Green2 = Reechantillonnage::biline(mVectArsenicImage[aK2].GChan.data(),
                                                                  mVectArsenicImage[aK2].SZ.x,
                                                                  mVectArsenicImage[aK2].SZ.y, pos2DOtherIm[aK2]);
                        double Blue2 = Reechantillonnage::biline(mVectArsenicImage[aK2].BChan.data(),
                                                                 mVectArsenicImage[aK2].SZ.x,
                                                                 mVectArsenicImage[aK2].SZ.y, pos2DOtherIm[aK2]);
                        mMatPtsHom.aMat[aK1].Pts.push_back(
                                pos2DPtIm1.mul(mVectArsenicImage[aK1].info3D->ResolImRefFromCapteur()));
                        if (Red1 > 0) {
                            mMatPtsHom.aMat[aK1].kR.push_back((1 + Red2 / Red1) / 2);
                        } else { mMatPtsHom.aMat[aK1].kR.push_back((1 + Red2) / 2); }
                        if (Green1 > 0) {
                            mMatPtsHom.aMat[aK1].kG.push_back((1 + Green2 / Green1) / 2);
                        } else { mMatPtsHom.aMat[aK1].kG.push_back((1 + Green2) / 2); }
                        if (Blue1 > 0) {
                            mMatPtsHom.aMat[aK1].kB.push_back((1 + Blue2 / Blue1) / 2);
                        } else { mMatPtsHom.aMat[aK1].kB.push_back((1 + Blue2) / 2); }
                        mMatPtsHom.aMat[aK1].OtherIm.push_back(aK2);
                        mMatPtsHom.aMat[aK1].SZ = mVectArsenicImage[aK1].SZ;
                    }
                }

            }
        }
    }
}

int PIMsCorrColor_main(int argc, char **argv) {
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
