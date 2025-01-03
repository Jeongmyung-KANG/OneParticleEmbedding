#ifndef __JMJETVEC_H__ 
#define __JMJETVEC_H__ 

#include <iostream>
#include "TMath.h"
#include "TObject.h"
#include <cassert>
class jmJetVec : public TObject {
    private : 
        Double_t __px; 
        Double_t __py;
        Double_t __pz;
        Double_t __pt;
        Double_t __e;
        Double_t __eta; 
        Double_t __phi;
        Double_t __A;
        Double_t __rho;
        Double_t __rhoM;
		Double_t __m;
		
		Double_t __trigPx;
		Double_t __trigPy; 
		Double_t __trigPz;
		Double_t __trigE; 

        Double_t __trigPt;
        Double_t __trigPhi;
        Double_t __trigEta;

        Double_t __trackSum;
		Bool_t __IsTrueJet;	 
		int __jetIndex;  
		int __constiSize;
		std::vector<double> __constiVectorPhi; 
		std::vector<double> __constiVectorEta;  
		std::vector<double> __constiVectorPt; 
    public :

        jmJetVec() : __px(0), __py(0), __pz(0), __pt(0), __e(0), __m(0), __eta(0), __phi(0), __A(0), __rho(0), __rhoM(0), __trigPt(0), __trigPhi(0), __trigEta(0), __trackSum(0), __trigPx(0), __trigPy(0), __trigPz(0), __trigE(0), __IsTrueJet(0), __jetIndex(0), __constiSize(0), __constiVectorPhi(0), __constiVectorEta(0), __constiVectorPt(0){};
        jmJetVec(Double_t _px, Double_t _py, Double_t _pz, Double_t _e, Double_t _m, 
                 Double_t _phi, Double_t _eta, 
                 Double_t _A,
                 Double_t _rho,
                 Double_t _rhoM,
				 Double_t _trigPx,
				 Double_t _trigPy,
				 Double_t _trigPz,
				 Double_t _trigE,
                 Double_t _trigPt, Double_t _trigPhi, Double_t _trigEta,
                 Double_t _trackSum, Bool_t _IsTrueJet, int _jetIndex, int _constiSize, std::vector<double> _constiVectorPhi, std::vector<double> _constiVectorEta, std::vector<double> _constiVectorPt):
				__px(_px), __py(_py), __pz(_pz), __pt(TMath::Sqrt(_px*_px + _py*_py)), __e(_e), __m(_m),
                 __eta(_eta), __phi(_phi), __A(_A),
                __rho(_rho), __rhoM(_rhoM), 
                __trigPt(_trigPt), __trigPhi(_trigPhi), __trigEta(_trigEta),__trackSum(_trackSum), __trigPx(_trigPx), __trigPy(_trigPy), __trigPz(_trigPz), __trigE(_trigE), __IsTrueJet(_IsTrueJet), __jetIndex(_jetIndex), __constiSize(_constiSize), __constiVectorPhi(_constiVectorPhi), __constiVectorEta(_constiVectorEta), __constiVectorPt(_constiVectorPt){};

        ~jmJetVec() {};
        
        void setJetKin(Double_t _px, Double_t _py, Double_t _pz, Double_t _e, Double_t _m);
        void setJetPhiEta(Double_t _phi, Double_t _eta);
        void setA(Double_t _A);
        void setRho(Double_t _rho);
        void setRhoM(Double_t _rhoM);
        void setTriggerPtPhiEta(Double_t _pt, Double_t _phi, Double_t _eta);
		void setTriggerPxPyPzE(Double_t _px, Double_t _py, Double_t _pz, Double_t _e);
		void setTrackSum(Double_t _trackSum);  
		void setIsTrueJet(Bool_t _IsTrueJet); 
		void setJetIndex(int _jetIndex);
		void setConstiSize(int _ConstiSize); 
		void setConstiVectorPhi(std::vector<double> _constiVectorPhi);
		void setConstiVectorEta(std::vector<double> _constiVectorEta); 
		void setConstiVectorPt(std::vector<double> _constiVectorPt); 

		Double_t px() const { return __px; } 
        Double_t py() const { return __py; }
        Double_t pz() const { return __pz; }
        Double_t pt() const { return __pt; }
        Double_t e() const { return __e; }
		Double_t m() const {return __m;}
        Double_t eta() const { return __eta; }
        Double_t phi() const { return __phi; }
        Double_t A() const { return __A; }
        Double_t rho() const { return __rho; }
        Double_t rhoM() const { return __rhoM; }
        Double_t trigPt() const {return __trigPt;}
        Double_t trigPhi() const {return __trigPhi;}
        Double_t trigEta() const {return __trigEta;}
        Double_t trackSum() const {return __trackSum;}
		Double_t trigPx() const {return __trigPx;}
		Double_t trigPy() const {return __trigPy;}
		Double_t trigPz() const {return __trigPz;}
		Double_t trigE() const {return __trigE;}
		Bool_t IsTrueJet() const{return __IsTrueJet;} 
		int jetIndex() const{return __jetIndex;}  
		int constiSize() const{return __constiSize;} 
		std::vector<double> constiVectorPhi() const{return __constiVectorPhi;} 
		std::vector<double> constiVectorEta() const{return __constiVectorEta;} 
		std::vector<double> constiVectorPt() const{return __constiVectorPt;};

};

inline void jmJetVec::setJetKin(Double_t _px, Double_t _py, Double_t _pz, Double_t _e, Double_t _m){
    __px = _px;
    __py = _py;
    __pz = _pz;
    __e  = _e;
    __pt = TMath::Sqrt(_px*_px + _py*_py); 
	__m = _m;
}

inline void jmJetVec::setJetPhiEta(Double_t _phi, Double_t _eta){
    __phi = _phi;
    __eta = _eta;
}

inline void jmJetVec::setA(Double_t _A) { 
    __A = _A; 
}
inline void jmJetVec::setRho(Double_t _rho){
    __rho = _rho;
}
inline void jmJetVec::setRhoM(Double_t _rhoM){
    __rhoM = _rhoM;
}

inline void jmJetVec::setTriggerPxPyPzE(Double_t px, Double_t py, Double_t pz, Double_t e){
	__trigPx = px; 
	__trigPy = py; 
	__trigPz = pz;
	__trigE = e; 
}
inline void jmJetVec::setTriggerPtPhiEta(Double_t _pt, Double_t _phi, Double_t _eta){
    __trigPt = _pt;
    __trigPhi = _phi; 
    __trigEta = _eta;
}

inline void jmJetVec::setTrackSum(Double_t _trackSum){
	__trackSum = _trackSum;
} 

inline void jmJetVec::setIsTrueJet(Bool_t _IsTrueJet){
	__IsTrueJet = _IsTrueJet;
}

inline void jmJetVec::setJetIndex(int _jetIndex){
	__jetIndex = _jetIndex; 
}

inline void jmJetVec::setConstiSize(int _constiSize){
	__constiSize = _constiSize;
}

inline void jmJetVec::setConstiVectorPhi(std::vector<double> _constiVectorPhi){
	__constiVectorPhi = _constiVectorPhi; 
}

inline void jmJetVec::setConstiVectorEta(std::vector<double> _constiVectorEta){
	__constiVectorEta = _constiVectorEta; 
}
inline void jmJetVec::setConstiVectorPt(std::vector<double> _constiVectorPt){ 
	__constiVectorPt = _constiVectorPt; 
}
#endif // __JMJETVEC_H__
