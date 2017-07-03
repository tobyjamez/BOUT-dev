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

  //Now save and finish
  SAVE_ONCE4(ddy, ddy2, ddy3, ddy4);
  dump.write();
    
  BoutFinalise();

  return 0;
}

//   SteadyClock start1 = steady_clock::now();
//   for(int i=0; i<nrepeat; i++) {
//     //  ddy = DDY_yud(var);
//     ddy = DDY_inbuilt(var);
//   }
//   Duration elapsed1 = steady_clock::now() - start1;

//   SteadyClock start2;
//   Duration elapsed2;
  
//   bool usealigned=false;
//   Field3D ddy2;
//   if(usealigned){
//     // Change into field-aligned coordinates
//     Field3D var_aligned = mesh->toFieldAligned(var);
  
//     // var now field aligned
//     ddy2 = DDY_aligned(var_aligned);
  
//     // Shift back to orthogonal X-Z coordinates
//     ddy2 = mesh->fromFieldAligned(ddy2);
//   }else {
//     //    ddy2 = DDY_inbuilt(var);
//   start2 = steady_clock::now();
//   for(int i=0; i<nrepeat; i++) {

//   //ddy2 = DDY_valarray(var);
// }
//   elapsed2 = steady_clock::now() - start2;
//   }

//   SAVE_ONCE2(ddy, ddy2);



//   Field3D testing;
//   testing.allocate();
//   testing = 1.;
//   testing.get()[std::gslice(threeToCompound(0,mesh->ystart+1,0),
// 			    {mesh->LocalNx,mesh->LocalNy-2*mesh->ystart,mesh->LocalNz}, {mesh->LocalNz*mesh->LocalNy,mesh->LocalNz,1})] = 2.0;
//   testing(4,4,4) = 3.0;
//   SAVE_ONCE(testing);
    
//   dump.write();

//   output << endl;
//   output << "TIMING DDY(var) \n==================\n";
//   output << "Direct: " << elapsed1.count()/nrepeat << endl;
//   output << "Alternate: " << elapsed2.count()/nrepeat << endl;
//   output << "-------------------------------------\n" << endl;
