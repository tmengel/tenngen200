#include "TGParticle.h"

TENNGEN_BEGIN_NAMESPACE

double TGParticle::M(){ return getMass(); }
double TGParticle::Px(){  return pt*cos(phi);}
double TGParticle::Py(){ return pt*sin(phi); }
double TGParticle::Pz(){ return pt*sinh(eta); }
double TGParticle::E(){ return sqrt(pt*pt*cosh(eta)*cosh(eta)+(getMass()*getMass())); }

double TGParticle::Rap(){
    double E = sqrt(pt*pt*cosh(eta)*cosh(eta)+(getMass()*getMass()));
    double PZ = pt*sinh(eta);
    //const double c = 2.99792458e8;
    return 0.5*log((E+PZ*M_C)/(E-PZ*M_C));
}
  
double TGParticle::getMass(){
    double mass;
    if(kf==211||kf==-211){ mass= 0.13957;}
    else if(kf==321||kf==-321){mass= 0.49368;}
    else if(kf==2212||kf==-2212){mass=0.93827;}
    else {mass= 0.13957;} 
    return mass;
}

ostream& operator<<(ostream& os, const TGParticle& v) {
  os << fixed << setprecision(3) <<  v.pt << " "
     << v.eta << " " <<  v.phi << " " << v.kf << " ";
  return os;
}



bool operator==(const TGParticle & a, const TGParticle & b) {
  if (a.pt != b.pt) return false;
  if (a.eta != b.eta) return false;
  if (a.phi != b.phi) return false;
  if (a.kf != b.kf) return false;
  
  if (a.userindex    != b.userindex) return false;

  return true;
}

bool operator==(const TGParticle& a, const int type) {
  if (a.kf != type) return false;
  return true;
}

void TGPartList::setListUserIndex(int u){
    
    for(int i =0; i<this->size();i++){
     (particle.at(i)).setUserIndex(u);
    }
}


TGPartList& TGPartList::operator=( const TGPartList& List) {

  // Do not copy if same.
  if (this != &List) {
    // Reset all current info in the list.
    clear();
    // Copy all the particles one by one.
    for (int i = 0; i < List.size(); ++i) add( List[i] );
  }
  return *this;

}

TGPartList& TGPartList::operator+=( const TGPartList& addList) {

  TGParticle temp;
  for (int i = 0; i < addList.size(); ++i) {
    // temp = addList[i];
    add( addList[i] );
  }

  return *this;

}

TGPartList& TGPartList::operator+=(const TGParticle& addPart){
    add(addPart);
    return *this;
}

TGPartList operator+(const TGPartList &a, const TGPartList &b){
    TGPartList sumlist(a);
    for(int i = 0; i<b.size(); i++){
        sumlist.add(b[i]);
    }
    return sumlist;
}

TENNGEN_END_NAMESPACE