
#include <bout/varray.hxx>
#include <bout_types.hxx>

template<>
std::map< int, std::vector<VArray<BoutReal>::dataPtrType > > VArray<BoutReal>::store = {};

template<>
bool VArray<BoutReal>::use_store = true;
