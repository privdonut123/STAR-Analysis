/********************************************************************
* .sl73_x8664_gcc485/obj/StRoot/StMuDSTMaker/StMuDSTMaker_Cint.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error .sl73_x8664_gcc485/obj/StRoot/StMuDSTMaker/StMuDSTMaker_Cint.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtableStMuDSTMaker_Cint();
extern void G__cpp_setup_inheritanceStMuDSTMaker_Cint();
extern void G__cpp_setup_typetableStMuDSTMaker_Cint();
extern void G__cpp_setup_memvarStMuDSTMaker_Cint();
extern void G__cpp_setup_globalStMuDSTMaker_Cint();
extern void G__cpp_setup_memfuncStMuDSTMaker_Cint();
extern void G__cpp_setup_funcStMuDSTMaker_Cint();
extern void G__set_cpp_environmentStMuDSTMaker_Cint();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "StMuETofCollection.h"
#include "StMuEmcCollection.h"
#include "StMuFcsCollection.h"
#include "StMuFmsCollection.h"
#include "StMuFstCollection.h"
#include "StMuFttCollection.h"
#include "StMuFwdTrackCollection.h"
#include "StMuMtdCollection.h"
#include "StMuPmdCollection.h"
#include "StMuRHICfCollection.h"
#include "StMuRpsCollection.h"
#include "StMuTriggerIdCollection.h"
#include "SchedulerExample.h"
#include "StAddRunInfoMaker.h"
#include "StMuArrays.h"
#include "StMuBTofHit.h"
#include "StMuBTofPidTraits.h"
#include "StMuBTofUtil.h"
#include "StMuChainMaker.h"
#include "StMuCut.h"
#include "StMuDbReader.h"
#include "StMuDebug.h"
#include "StMuDst.h"
#include "StMuDst2StEventMaker.h"
#include "StMuDstFilterMaker.h"
#include "StMuDstMaker.h"
#include "StMuETofDigi.h"
#include "StMuETofHeader.h"
#include "StMuETofHit.h"
#include "StMuETofPidTraits.h"
#include "StMuEmcCluster.h"
#include "StMuEmcHit.h"
#include "StMuEmcPoint.h"
#include "StMuEmcTowerData.h"
#include "StMuEmcUtil.h"
#include "StMuEpdHit.h"
#include "StMuEpdUtil.h"
#include "StMuEvent.h"
#include "StMuFcsCluster.h"
#include "StMuFcsHit.h"
#include "StMuFcsInfo.h"
#include "StMuFcsPoint.h"
#include "StMuFcsUtil.h"
#include "StMuFgtAdc.h"
#include "StMuFgtCluster.h"
#include "StMuFgtStrip.h"
#include "StMuFgtStripAssociation.h"
#include "StMuFilter.h"
#include "StMuFmsCluster.h"
#include "StMuFmsHit.h"
#include "StMuFmsInfo.h"
#include "StMuFmsPoint.h"
#include "StMuFmsUtil.h"
#include "StMuFstHit.h"
#include "StMuFstRawHit.h"
#include "StMuFstUtil.h"
#include "StMuFttCluster.h"
#include "StMuFttPoint.h"
#include "StMuFttRawHit.h"
#include "StMuFttUtil.h"
#include "StMuFwdTrack.h"
#include "StMuFwdTrackUtil.h"
#include "StMuHelix.h"
#include "StMuIOMaker.h"
#include "StMuL3EventSummary.h"
#include "StMuL3Filter.h"
#include "StMuMcTrack.h"
#include "StMuMcVertex.h"
#include "StMuMomentumShiftMaker.h"
#include "StMuMtdHeader.h"
#include "StMuMtdHit.h"
#include "StMuMtdPidTraits.h"
#include "StMuMtdRawHit.h"
#include "StMuPmdCluster.h"
#include "StMuPmdHit.h"
#include "StMuPmdUtil.h"
#include "StMuPrimaryTrackCovariance.h"
#include "StMuPrimaryVertex.h"
#include "StMuProbPidTraits.h"
#include "StMuRHICfHit.h"
#include "StMuRHICfPoint.h"
#include "StMuRHICfRawHit.h"
#include "StMuRHICfUtil.h"
#include "StMuRpsTrack.h"
#include "StMuRpsTrackPoint.h"
#include "StMuTimer.h"
#include "StMuTofHit.h"
#include "StMuTofUtil.h"
#include "StMuTrack.h"
#include "StuDraw3DMuEvent.h"
#include "EztEmcRawData.h"
#include "EztEventHeader.h"
#include "EztFpdBlob.h"
#include "EztTrigBlob.h"
#include "StMuEzTree.h"
#include "StTriggerDataMother.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__StMuDSTMaker_CintLN_TClass;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TBuffer;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMemberInspector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TObject;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TNamed;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TString;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEshortcOallocatorlEshortgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlElongcOallocatorlElonggRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEunsignedsPcharcOallocatorlEunsignedsPchargRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEdoublecOallocatorlEdoublegRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEboolcOallocatorlEboolgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_string;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TObjArray;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TClonesArray;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StETofCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuETofHeader;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuETofDigi;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuETofHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuETofCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuETofHeadercOallocatorlEStMuETofHeadergRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuETofHeadercOallocatorlEStMuETofHeadergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuETofDigicOallocatorlEStMuETofDigigRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuETofDigicOallocatorlEStMuETofDigigRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuETofHitcOallocatorlEStMuETofHitgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuETofHitcOallocatorlEStMuETofHitgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEmcHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TArrayS;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEmcCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEmcPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StDetectorId;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTrackType;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTrackFinderMethod;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTrackFittingMethod;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEmcCrateStatus;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StVertexFinderId;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEmcTowerData;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEmcTowerDatacLcLdA;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEmcCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFcsHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFcsCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFcsPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFcsInfo;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFcsCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TDataSet;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StObject;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFmsHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFmsCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFmsPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFmsInfo;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFmsCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFstRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFstHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFstCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFttRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFttCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFttPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFttCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFwdTrack;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFwdTrackCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMtdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMtdHeader;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMtdHeader;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMtdHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMtdRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMtdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuMtdHeadercOallocatorlEStMuMtdHeadergRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuMtdHeadercOallocatorlEStMuMtdHeadergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuMtdHitcOallocatorlEStMuMtdHitgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuMtdHitcOallocatorlEStMuMtdHitgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuMtdRawHitcOallocatorlEStMuMtdRawHitgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuMtdRawHitcOallocatorlEStMuMtdRawHitgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuPmdCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuPmdHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuPmdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRHICfRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRHICfHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRHICfPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRHICfCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_pairlEstringcOlonggR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_pairlEdoublecOdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TElementActionTlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TElementPosActionTlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TElementActionTlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TElementPosActionTlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTRow_constlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTRowlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTDiag_constlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTColumn_constlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTFlat_constlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSub_constlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseRow_constlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseDiag_constlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTColumnlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTDiaglEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTFlatlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSublEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseRowlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseDiaglEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TVector3;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRpsTrackPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsTrackPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsTrackPointcLcLStMuRpsTrackPointQuality;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsTrackPointcLcLdA;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TFile;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TRef;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsTrack;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsTrackcLcLStMuRpsTrackType;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsTrackcLcLStMuRpsAngles;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsTrackcLcLdA;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRpsCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRpsCollectioncLcLdA;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuRpsTrackPointmUcOallocatorlEStMuRpsTrackPointmUgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuRpsTrackPointmUcOallocatorlEStMuRpsTrackPointmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuRpsTrackmUcOallocatorlEStMuRpsTrackmUgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuRpsTrackmUcOallocatorlEStMuRpsTrackmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEStRpsTrackPointmUcOStMuRpsTrackPointmUcOlesslEStRpsTrackPointmUgRcOallocatorlEpairlEStRpsTrackPointmUsPconstcOStMuRpsTrackPointmUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTriggerId;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTriggerIdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuTriggerIdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TArrayC;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TArrayI;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TArrayF;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TH1D;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TChain;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TTree;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDstMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDst;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TH3D;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_SchedulerExample;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEvent;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StAddRunInfoMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_emcTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_fgtTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_fmsTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_rhicfTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_MCTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_muDstTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_pmdTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_tofTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_btofTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_etofTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_mtdTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_epdTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_eztTypes;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_NARRAYS;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuArrays;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StThreeVectorlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StThreeVectorlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StBTofHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEvent;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuPrimaryVertex;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuTrack;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRichSpectra;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StDetectorState;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StL3AlgorithmInfo;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StStrangeEvMuDst;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StV0MuDst;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StXiMuDst;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StKinkMuDst;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StV0Mc;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StXiMc;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StKinkMc;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StStrangeAssoc;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TCut;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTriggerData;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTrack;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTrackGeometry;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEmcCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFmsCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRHICfCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuTofHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTofData;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTofRawData;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StBTofCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuBTofHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StBTofRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StBTofHeader;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEpdHitCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEpdHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_EztEventHeader;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_EztTrigBlob;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_EztFpdBlob;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_EztEmcRawData;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StDcaGeometry;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuPrimaryTrackCovariance;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StPhysicalHelix;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDebug;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuIOMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuL3EventSummary;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEventInfo;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRunInfo;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEventSummary;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StVpdTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMtdTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StCtbTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StZdcTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StBbcTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEmcTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFpdTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFmsTriggerDetector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFpdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StL0Trigger;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuCut;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTofCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMomentumShiftMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuHelix;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuProbPidTraits;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StBTofPidTraits;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuBTofPidTraits;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StETofPidTraits;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuETofPidTraits;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMtdPidTraits;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMtdPidTraits;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StLorentzVectorlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMatrixlEfloatgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTrackTopologyMap;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTRow_constlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTRowlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTDiag_constlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTColumn_constlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTFlat_constlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSub_constlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseRow_constlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseDiag_constlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTColumnlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTDiaglEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTFlatlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSublEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseRowlEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TMatrixTSparseDiaglEdoublegR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlETHelixKFitterAuxcOallocatorlETHelixKFitterAuxgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlETHelixKFitterAuxcOallocatorlETHelixKFitterAuxgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StVertex;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEmcGeom;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StuProbabilityPidAlgorithm;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDstFilterMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuBTofHitCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuBTofUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEpairlEstringcOlonggRcOallocatorlEpairlEstringcOlonggRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEpairlEstringcOlonggRcOallocatorlEpairlEstringcOlonggRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDbReader;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuChainMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StV0Vertex;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StXiVertex;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StKinkVertex;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEpairlEstringcOintgRcOallocatorlEpairlEstringcOintgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEpairlEstringcOintgRcOallocatorlEpairlEstringcOintgRsPgRsPgRcLcLiterator;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEpairlEstringcOintgRcOallocatorlEpairlEstringcOintgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDst2StEventMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRTSBaseMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlETStringcOallocatorlETStringgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TVirtualPad;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StUKey;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StIOInterFace;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_listlEunsignedsPshortcOallocatorlEunsignedsPshortgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFilter;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_BetheBloch;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuL3Filter;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTrackNode;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StIOMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTreeMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StStrangeMuDstMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEmcUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFmsUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuRHICfUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFcsUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFttUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFstUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFwdTrackUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEpdHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuEpdHitmUcOallocatorlEStMuEpdHitmUgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuEpdHitmUcOallocatorlEStMuEpdHitmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEpdUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuPmdUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuTofHitCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuTofUtil;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEzTree;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TEventList;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDstMakercLcLioMode;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuDstMakercLcLioNameMode;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StETofDigi;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StETofHeader;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEunsignedsPintcOunsignedsPlongcOlesslEunsignedsPintgRcOallocatorlEpairlEconstsPunsignedsPintcOunsignedsPlonggRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEunsignedsPintcOunsignedsPlongsPlongcOlesslEunsignedsPintgRcOallocatorlEpairlEconstsPunsignedsPintcOunsignedsPlongsPlonggRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StETofHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEpdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TRefArray;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_TLorentzVector;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFcsCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEconstsPStFcsClustermUcOStMuFcsClustermUcOlesslEconstsPStFcsClustermUgRcOallocatorlEpairlEconstsPStFcsClustermUsPconstcOStMuFcsClustermUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEconstsPStFcsHitmUcOStMuFcsHitmUcOlesslEconstsPStFcsHitmUgRcOallocatorlEpairlEconstsPStFcsHitmUsPconstcOStMuFcsHitmUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEconstsPStFcsPointmUcOStMuFcsPointmUcOlesslEconstsPStFcsPointmUgRcOallocatorlEpairlEconstsPStFcsPointmUsPconstcOStMuFcsPointmUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFgtAdc;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFgtHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFgtCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFgtStrip;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFgtStrip;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFgtStripAssociation;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFmsCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFmsPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFmsDbMaker;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFstHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFstRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFstHitCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFstEvtCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFttCluster;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFttPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFttRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFttCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEconstsPStFttRawHitmUcOTObjectmUcOlesslEconstsPStFttRawHitmUgRcOallocatorlEpairlEconstsPStFttRawHitmUsPconstcOTObjectmUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEconstsPStFttClustermUcOTObjectmUcOlesslEconstsPStFttClustermUgRcOallocatorlEpairlEconstsPStFttClustermUsPconstcOTObjectmUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEconstsPStFttPointmUcOTObjectmUcOlesslEconstsPStFttPointmUgRcOallocatorlEpairlEconstsPStFttPointmUsPconstcOTObjectmUgRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFwdTrack;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFwdTrackProjection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuFwdTrackSeedPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuFwdTrackProjectioncOallocatorlEStMuFwdTrackProjectiongRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuFwdTrackProjectioncOallocatorlEStMuFwdTrackProjectiongRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_vectorlEStMuFwdTrackSeedPointcOallocatorlEStMuFwdTrackSeedPointgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_reverse_iteratorlEvectorlEStMuFwdTrackSeedPointcOallocatorlEStMuFwdTrackSeedPointgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StFwdTrackCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_g2t_track_st;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMcTrack;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMcTrackcLcLEHIT;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_g2t_vertex_st;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuMcVertex;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMtdHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMtdRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StPhmdCollection;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StPrimaryVertex;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRHICfRawHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRHICfHit;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StRHICfPoint;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuTimer;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_EDraw3DStyle;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StDraw3D;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_maplEEDraw3DStylecOStDraw3DStylecOlesslEEDraw3DStylegRcOallocatorlEpairlEconstsPEDraw3DStylecOStDraw3DStylegRsPgRsPgR;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_EEmcGeomSimple;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StuDraw3DMuEvent;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_EztEmcRawDatacLcLdA;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StEmcRawData;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StMuEzTreecLcLdA;
extern G__linked_taginfo G__StMuDSTMaker_CintLN_StTriggerDataMother;

/* STUB derived class for protected member access */