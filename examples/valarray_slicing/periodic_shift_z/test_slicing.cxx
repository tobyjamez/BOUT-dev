#include <bout.hxx>

#include <chrono>

typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
typedef std::chrono::duration<double> Duration;
using namespace std::chrono;

#include <valarray>
#include <slicingutils.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D indices=0.;
  Field3D indicesStore=0.;
  Field3D indicesLoad=0.;
  Field3D indexZ=0.;

  for(const auto& i: indices){
    const int jzm = i.z>0? i.z-1 : mesh->LocalNz-1;
    const int jzp = i.z<mesh->LocalNz-1? i.z+1 : 0;

    indices[i] = threeToCompound(i.x,i.y,i.z); //Single index at this point
    indexZ[i] = i.z; //Z index at this point
    
    //Single index at which to store this point -- as we're shifting z by one point left
    //we're storing at jz-1 so use jzm to calculate single index of data destination
    indicesStore[i] = threeToCompound(i.x,i.y,jzm);
    //Single index at which to load this point from -- as we're shifting z by one point left
    //we're loading from jz+1 so use jzp to calculate single index of data source
    indicesLoad[i] = threeToCompound(i.x,i.y,jzp);
  }
  SAVE_ONCE(indices);

  ///////////////////////////////////////////////
  // Method 1 : Make slices and manually merge
  //    Somewhat clunky and not very general.
  ///////////////////////////////////////////////

  //Now we make two slices, one for jz>=1 and the other for jz==0
  auto firstPart =  getRegionSlice(0,mesh->LocalNx-1,0,mesh->LocalNy-1,1,mesh->LocalNz-1);
  auto secondPart = getRegionSlice(0,mesh->LocalNx-1,0,mesh->LocalNy-1,0,0);

  //Slice the index array
  const auto first = std::valarray<BoutReal>(indices.get()[firstPart]);
  const auto second = std::valarray<BoutReal>(indices.get()[secondPart]);

  const int firstSize = first.size();
  const int secondSize = second.size();


  //This will hold our indices.
  auto merged = std::valarray<std::size_t>(indices.get().size());
  merged = 0;

  //Sanity check
  ASSERT0(merged.size()==firstSize+secondSize); 

  //Now we insert the indices from the two slices into 
  //a single index valarray. Note we have to interleave the
  //two sections appropriately in order to get things stored
  //in the right place
  int counter=0;
  int firstCount=0;
  int secondCount=0;
  do{
    for(int i=0; i<mesh->LocalNz-1; i++){
      merged[counter] = first[firstCount];
      firstCount++;
      counter++;
    }
    for(int i=0; i<1; i++){
      merged[counter] = second[secondCount];
      secondCount++;
      counter++;
    }
  }while(counter<merged.size());

  //Sanity check
  ASSERT0(merged.size()==counter); 

  Field3D ind=-1.;
  ind.get() = indexZ.get()[merged];
  SAVE_ONCE(ind);

  ///////////////////////////////////////////////
  // Method 2 : Just use "target" indices
  //    Much simpler but less clean setup 
  //    (no clear relation to slices etc.)
  //    Currently used to change indices of store
  //    but could be done to change indices of source
  //    (i.e. instead of saying where data is going
  //    say where it has come from, probably better
  //    as can be used on RHS more clearly).
  ///////////////////////////////////////////////

  //This will hold our indices.
  auto merged2 = std::valarray<std::size_t>(indices.get().size());
  merged2 = 0;

  //Sanity check
  ASSERT0(merged2.size()==firstSize+secondSize); 

  counter=0;
  for(const auto i: indicesStore){
    merged2[counter]  = indicesStore[i];
    counter++;
  }

  Field3D ind2=-1.;
  ind2.get()[merged2] = indexZ.get();
  SAVE_ONCE(ind2);

  ///////////////////////////////////////////////
  // Method 3 : Just use "source" indices
  //    Much simpler but less clean setup 
  //    (no clear relation to slices etc.)
  ///////////////////////////////////////////////

  //This will hold our indices.
  auto merged3 = std::valarray<std::size_t>(indices.get().size());
  merged3 = 0;

  //Sanity check
  ASSERT0(merged3.size()==firstSize+secondSize); 

  counter=0;
  for(const auto i: indicesLoad){
    merged3[counter]  = indicesLoad[i];
    counter++;
  }

  Field3D ind3=-1.;
  ind3.get() = indexZ.get()[merged3];
  SAVE_ONCE(ind3);

  dump.write();
  BoutFinalise();

  return 0;
}
