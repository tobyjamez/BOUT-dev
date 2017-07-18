#ifndef __SLICINGUTILS_H__
#define __SLICINGUTILS_H__

#include <valarray>

/*!
 * Convert x/y/z index tuple to single compound index
 *
 * @param[in] jx The index in x
 * @param[in] jy The index in y
 * @param[in] jz The index in z
 *
 * @param[out] The flattened compound index
 */
const int threeToCompound(const int &jx, const int &jy, const int &jz);

/*!
 * Get a std::gslice object representing a custom region
 *
 * @param[in] xStart The first index in x to include in the region (inclusive)
 * @param[in] xEnd   The last index in x to include in the region  (inclusive)
 * @param[in] yStart The first indey in y to include in the region (inclusive)
 * @param[in] yEnd   The last indey in y to include in the region  (inclusive)
 * @param[in] zStart The first indez in z to include in the region (inclusive)
 * @param[in] zEnd   The last indez in z to include in the region  (inclusive)
 * 
 * @param[out] gslice The std::gslice object representing the region.
 */
const std::gslice getRegionSlice (const int &xStart, const int &xEnd,
				  const int &yStart, const int &yEnd,
				  const int &zStart, const int &zEnd);
#endif
