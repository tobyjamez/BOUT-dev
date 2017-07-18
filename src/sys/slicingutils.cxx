#include <slicingutils.hxx>
#include <bout/mesh.hxx>

//Convert x/y/z index tuple to single compound index
const int threeToCompound(const int &jx, const int &jy, const int &jz) {
  //Dependence on global mesh --> suggests maybe this should be a mesh member
  //or at least member of something that is constructed with a mesh.
  return (jx*mesh->LocalNy + jy)*mesh->LocalNz + jz;
}

//Get gslice representing specified region/range
const std::gslice getRegionSlice (const int &xStart, const int &xEnd,
				  const int &yStart, const int &yEnd,
				  const int &zStart, const int &zEnd) {

  //Get the first point in range in "single index space"
  const int startPoint = threeToCompound(xStart, yStart, zStart);

  //Calculate the effective extent of each of the three dimensions
  //Inclusive of start+end, remove +1 to make exclusive of end point.
  const unsigned long nxTot = xEnd - xStart + 1; 
  const unsigned long nyTot = yEnd - yStart + 1;
  const unsigned long nzTot = zEnd - zStart + 1;

  //Describe the logical stride in each of the three dimensions for the
  //*stored* data (i.e. not the sliced strides).
  const unsigned long zStride = 1;
  const unsigned long yStride = zStride * mesh->LocalNz;
  const unsigned long xStride = yStride * mesh->LocalNy;

  //Construct and return general slice (gslice).
  return std::gslice(startPoint, {nxTot, nyTot, nzTot}, {xStride, yStride, zStride});
}
