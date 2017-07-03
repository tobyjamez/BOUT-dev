#include <bout.hxx>
#include <bout/paralleltransform.hxx>
#include <derivs.hxx>

#include <chrono>

typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
typedef std::chrono::duration<double> Duration;
using namespace std::chrono;

#include <valarray>

//Convert x/y/z index tuple to single compound index
const int threeToCompound(const int &jx, const int &jy, const int &jz) {
  return (jx*mesh->LocalNy + jy)*mesh->LocalNz + jz;
}

// Y derivative using yup() and ydown() fields
const Field3D DDY_yud(const Field3D &f) {
  Field3D result;
  result.allocate();

  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++)
    result(i,j,k) = 0.5*(f.yup()(i,j+1,k) - f.ydown()(i,j-1,k));

  return result;
}

// Y derivative assuming field is aligned in Y
const Field3D DDY_aligned(const Field3D &f) {
  Field3D result;
  result.allocate();
  
  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++)
  	result(i,j,k) = 0.5*(f(i,j+1,k) - f(i,j-1,k));

  return result;
}

// Y derivative inbuilt
const Field3D DDY_inbuilt(const Field3D &f) {
  Field3D result;
  result.allocate();
  result = DDY(f);
  return result;
}

// Y derivative valarray
const Field3D DDY_valarray(const Field3D &f, const std::gslice &yplus,const std::gslice &yminus,const std::gslice &ycen) {
  Field3D result;
  result.allocate();
  result.get()[ycen]  = f.yup().get()  [yplus];
  result.get()[ycen] -= f.ydown().get()[yminus];
  result.get() *= 0.5;
  return result;
}

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  ShiftedMetric s(*mesh);

  // Read variable from mesh
  Field3D var;
  mesh->get(var, "var");
  
  // Var starts in orthogonal X-Z coordinates

  // Calculate yup and ydown
  s.calcYUpDown(var);

  // Try each approach once first in order to "warm up"

  /////////////////////////////////////////////////
  // Calculate d/dy using yup() and ydown() fields
  Field3D ddy = DDY_yud(var);
  //-----------------------------------------------

  /////////////////////////////////////////////////
  // Use the "inbuilt" DDY method
  Field3D ddy2 = DDY_inbuilt(var);
  //-----------------------------------------------

  /////////////////////////////////////////////////
  // Use a "direct" field aligned method
  // Change into field-aligned coordinates
  Field3D var_aligned = mesh->toFieldAligned(var);
  
  // var now field aligned
  Field3D ddy3 = DDY_aligned(var_aligned);
  
  // Shift back to orthogonal X-Z coordinates
  ddy3 = mesh->fromFieldAligned(ddy3);
  //-----------------------------------------------

  /////////////////////////////////////////////////
  // Use valarray with slicing method
  //First construct slices
  const std::gslice yPlus(threeToCompound(0,mesh->ystart+1,0),//Offset jy+1
			  {mesh->LocalNx,mesh->LocalNy-2*mesh->ystart,mesh->LocalNz}, 
			  {mesh->LocalNz*mesh->LocalNy,mesh->LocalNz,1});

  const std::gslice yMinus(threeToCompound(0,mesh->ystart-1,0),//Offset jy-1
			   {mesh->LocalNx,mesh->LocalNy-2*mesh->ystart,mesh->LocalNz}, //How many points in {x,y,z} = {nx,ny-2*MYG,nz}
			   {mesh->LocalNz*mesh->LocalNy,mesh->LocalNz,1}); //Stride in {x,y,z}={nz*ny,nz,1}

  const std::gslice yCen(threeToCompound(0,mesh->ystart,0),//No offset
			 {mesh->LocalNx,mesh->LocalNy-2*mesh->ystart,mesh->LocalNz},
			 {mesh->LocalNz*mesh->LocalNy,mesh->LocalNz,1});

  Field3D ddy4 = DDY_valarray(var,yPlus,yMinus,yCen);
  //-----------------------------------------------

  
  const int nrepeat=100;
  SteadyClock start;
  
  start = steady_clock::now();
  for(int i=0; i<nrepeat; i++){
    ddy = DDY_yud(var);
  };
  Duration time1 = steady_clock::now() - start;

  start = steady_clock::now();
  for(int i=0; i<nrepeat; i++){
    ddy2 = DDY_inbuilt(var);
  };
  Duration time2 = steady_clock::now() - start;

  start = steady_clock::now();
  for(int i=0; i<nrepeat; i++){
    var_aligned = mesh->toFieldAligned(var);
    ddy3 = DDY_aligned(var_aligned);
    ddy3 = mesh->fromFieldAligned(ddy3);
  };
  Duration time3 = steady_clock::now() - start;

  start = steady_clock::now();
  for(int i=0; i<nrepeat; i++){
    ddy = DDY_valarray(var,yPlus,yMinus,yCen);
  };
  Duration time4 = steady_clock::now() - start;

  output << endl;
  output << "TIMING DDY(var) per call \n==================\n";
  output << "Direct with yup/ydown: " << time1.count()/nrepeat << endl;
  output << "Inbuilt ddy: " << time2.count()/nrepeat << endl;
  output << "Direct to/from aligned: " << time3.count()/nrepeat << endl;
  output << "Valarray slices: " << time4.count()/nrepeat << endl;
  output << "-------------------------------------\n" << endl;
    
  BoutFinalise();

  return 0;
}
