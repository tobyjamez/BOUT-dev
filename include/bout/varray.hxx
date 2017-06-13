#ifndef __VARRAY_H__
#define __VARRAY_H__

#include <valarray>

template<typename T>
class VArray : public std::valarray<T> {
public:
  using std::valarray<T>::valarray;

  // typedef T data_type;

  /*!
   * Release data. After this the Array is empty and any data access
   * will be invalid
   */
  void clear() {
    //std::valarray<T>::resize(0);#
    this->resize(0);
  }
  
  /*!
   * Returns true if the Array is empty
   */
  bool empty() const {
    //    return std::valarray<T>::size()==0;
    return this->size()==0;
  }
  
  /*!
   * Returns true if the data is unique to this Array.
   * 
   */
  bool unique() const {
    return true;
  }

  /*!
   * Ensures that this Array does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    return;
  }

private:
 
};

#endif // __VARRAY_H__

