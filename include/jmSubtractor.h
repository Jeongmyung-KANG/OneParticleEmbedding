#ifndef __JMSUBTRACTOR_H__
#define __JMSUBTRACTOR_H__

#include "TMath.h" 

#include "jmJetVec.h"

class jmSubtractor {

    private : 
    jmJetVec *inputJet; 
    jmJetVec *subtractedJet; 
    jmJetVec *subtractingArea;

    Double_t _rho; 
    Double_t _rhoM; 
    Double_t _m2; 
    Double_t _pt2;
    Double_t _scalar_mc; 
    Double_t _scalar_ptc;
    Double_t _A2; 
    Double_t _Apt2; 

    bool isNegativeM2 = false; 
    bool isNegativeCorrVal; 

    public : 
    jmSubtractor (jmJetVec *jet, Double_t rho, Double_t rhoM) : 
        inputJet(jet), subtractedJet(0), subtractingArea(0), _rho(rho), _rhoM(rhoM), _m2(0), _pt2(0), _scalar_mc(0), _scalar_ptc(0), _A2(0), _Apt2(0), isNegativeCorrVal(0) {
            setSubtractingArea();
            setSubtractedJet();
            calProduct();
            calScalarReconstruction(); 
        }; 
        
    virtual ~jmSubtractor() {}; 

    void setSubtractingArea();
    void setSubtractedJet();
    void calProduct ();
    void calScalarReconstruction(); 
    void resetRejectedJet();

    jmJetVec *GetSubtractedJet() const {return subtractedJet;} 
    Double_t m2() const {return _m2;}
    Double_t m();
    Double_t scalar_mc() const {return _scalar_mc;}
    Double_t scalar_ptc() const {return _scalar_ptc;}
    Double_t A2() const {return _A2;}
    Double_t Apt2() const {return _Apt2;}
};

inline void jmSubtractor::setSubtractingArea(){

    subtractingArea = new jmJetVec(); 

    Double_t ax = inputJet->Ax();
    Double_t ay = inputJet->Ay();
    Double_t az = inputJet->Az();
    Double_t ae = inputJet->Ae();

    Double_t f_ax = ax * _rho; 
    Double_t f_ay = ay * _rho; 
    Double_t f_az = az * _rho; //Double_t f_az = az * (_rho + _rhoM); 
    Double_t f_ae = ae * _rho; //Double_t f_ae = ae * (_rho + _rhoM); 

    //cout << ax << " " << ay << " " << az << " " << ae << " " <<  f_ax << " " << f_ay << " " << f_az << " " << f_ae << endl; 
    subtractingArea->setJetPxPyPzE(f_ax, f_ay, f_az, f_ae); 
}

inline void jmSubtractor::setSubtractedJet(){
    
    isNegativeCorrVal = false;

    subtractedJet = new jmJetVec(); 
    
    Double_t corr_px = inputJet->px() - subtractingArea->px(); 
    Double_t corr_py = inputJet->py() - subtractingArea->py(); 
    Double_t corr_pz = inputJet->pz() - subtractingArea->pz(); 
    Double_t corr_e = inputJet->e() - subtractingArea->e(); 
    
    if (corr_px < 0 || corr_py < 0 || corr_pz < 0 || corr_e < 0)  isNegativeCorrVal = false; 
    subtractedJet->setJetPxPyPzE(corr_px, corr_py, corr_pz, corr_e); 
}

inline void jmSubtractor::calProduct(){
    Double_t px2 = subtractedJet->px() * subtractedJet->px(); 
    Double_t py2 = subtractedJet->py() * subtractedJet->py();
    Double_t pz2 = subtractedJet->pz() * subtractedJet->pz();
    Double_t e2 = subtractedJet->e() * subtractedJet->e(); 

    _m2 = (e2 - px2 - py2 - pz2); 
    if (_m2 < 0) isNegativeM2 = true; 
    else if (_m2 > 0) isNegativeM2 = false; 
}

inline void jmSubtractor::calScalarReconstruction() { 
    
    Double_t f_px = subtractingArea->px(); 
    Double_t f_py = subtractingArea->py(); 
    Double_t f_pz = subtractingArea->pz(); 
    Double_t f_e = subtractingArea->e(); 

    
    Double_t f_rhoApt2 = f_px * f_px + f_py * f_py; 
    Double_t f_rhoMAm2  = f_e * f_e - (f_px * f_px + f_py * f_py + f_pz * f_pz); 
    _A2 = f_rhoMAm2; 
    _Apt2 = f_rhoApt2;
    if (f_rhoApt2 > 0) {
        Double_t f_rhoApt = TMath::Sqrt(f_rhoApt2);
        _scalar_ptc = inputJet->pt() - f_rhoApt; 
    }

    else if (f_rhoApt2 <= 0) { 
        //_scalar_ptc = inputJet->pt(); 
        _scalar_ptc = 0;
    }

    if (f_rhoMAm2 > 0){ 
        Double_t f_rhoMAm = TMath::Sqrt(f_rhoMAm2); 
        _scalar_mc = inputJet->m() - f_rhoMAm;
        //cout << "Scalar subtraction result for netative m2 : " << "input jet mass : " << inputJet->m() << " correction term : " << f_rhoMAm << endl; 
    }

    else if (f_rhoMAm2 < 0) { 
        //_scalar_mc = inputJet->m();
        _scalar_mc = 0;
    }
}

inline Double_t jmSubtractor::m(){ 
    if (_m2 < 0 || isNegativeCorrVal){
        return 0; 
    } 

    else {
        return TMath::Sqrt(_m2); 
    }
}

#endif // __JMSUBTRACTOR_H__
