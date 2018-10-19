
#include <bout/physicsmodel.hxx>
#include <bout/deriv_methods.hxx>

class WaveTest : public PhysicsModel {
public:
  int init(bool restarting) {
    a = 1.0;
    b = 2.0;
    c = 3.0;
    solver->add(a, "a");
    solver->add(b, "b");
    return 0;
  }
  
  int rhs(BoutReal time) {
    //Standard
    tester1.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);
    tester2.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);
    tester3.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);
    tester4.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);
    tester5.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);
    tester6.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);
    tester7.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);
    tester8.operator()<DIRECTION::X>(a, b, RGN_NOBNDRY);        
      
    //Upwind
    tester9.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    tester10.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    tester11.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    tester12.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    tester13.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    tester14.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    
    //Flux
    tester15.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    tester16.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);
    tester17.operator()<DIRECTION::X>(a, c, b, RGN_NOBNDRY);    
    
    ddt(a) = 0.;
    ddt(b) = 0.;
    return 0;
  }
  
private:
  Field3D a, b, c;
  //Standard
  const DDX_C2 tester1{};
  const DDX_C4 tester2{};
  const DDX_CWENO2 tester3{};
  const DDX_S2 tester4{};
  const D2DX2_C2 tester5{};
  const D2DX2_C4 tester6{};
  const DDX_CWENO3 tester7{};
  const D4DX4_C2 tester8{};
  //Upwind
  const VDDX_C2 tester9{};
  const VDDX_C4 tester10{};
  const VDDX_U1 tester11{};
  const VDDX_U2 tester12{};
  const VDDX_U3 tester13{};
  const VDDX_WENO3 tester14{};  
  //Flux
  const FDDX_U1 tester15{};
  const FDDX_C2 tester16{};
  const FDDX_C4 tester17{};  
};


BOUTMAIN(WaveTest);
