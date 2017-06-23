#ifndef __VARRAY_H__
#define __VARRAY_H__

#include <map>
#include <vector>
#include <valarray>
#include <memory>
#include "output.hxx"

template<typename T>
class VArray {
public:
  typedef T data_type;
  
  VArray() : ptr(nullptr) {}
  
  VArray(int len) {
    ptr = get(len);
  }
  
  ~VArray() {
    release(ptr);
  }
  
  VArray(const VArray &other) {
    ptr = other.ptr;
  }
    
  VArray(VArray&& other) {
    ptr = other.ptr; //Should use std::move?
    other.ptr.reset(); //Release pointer and set it to nullptr, could probably call release
  }

  VArray& operator=(const VArray &other) {
    dataPtrType  old = ptr;

    // Add reference
    ptr = other.ptr;

    // Release the old data
    release(old);
    
    return *this;
  }

  VArray& operator=(VArray &&other) {
    dataPtrType  old = ptr;

    ptr = other.ptr;
    other.ptr = nullptr;

    release(old);
    
    return *this;
  }

  void clear() {
    release(ptr);
  }
  
  bool empty() const {
    return !ptr;
  }

  int size() const {
    if(!ptr)
      return 0;

    return ptr->size();
  }
  
  bool unique() const {
    return ptr.use_count() == 1;
  }

  void ensureUnique() {
    if(!ptr || unique())
      return;

    // Get a new (unique) block of data
    dataPtrType p = get(size());
  
    //Make copy of the underlying data
    p->operator=((*ptr));

    //Update the local pointer and release old
    //could probably just do ptr=p as shared_ptr should
    //handle the rest.
    dataPtrType old = ptr;
    ptr = std::move(p);
    release(old);
  }

  //////////////////////////////////////////////////////////
  // Element access

  /*!
   * Access a data element. This will fail if the VArray is empty (so ptr is null),
   * or if ind is out of bounds. For efficiency no checking is performed,
   * so the user should perform checks.
   */
  T& operator[](int ind) {
    return ptr->operator[](ind);
  }
  const T& operator[](int ind) const {
    return ptr->operator[](ind);
  }

  static void cleanup() {
    for(auto &p : store) {
      auto &v = p.second;
      for(dataPtrType a : v) {
        a.reset();
      }
      v.clear();
    }
    // Don't use the store anymore
    use_store = false;
  }

  //Type defs to help keep things brief
  typedef std::valarray<T> dataBlock;
  typedef std::shared_ptr<dataBlock>  dataPtrType;

  dataPtrType ptr; //Shared_ptr to the valarray container

private:
  
  /*!
   * This maps from array size (int) to vectors of pointers to VArrayData objects
   * For each data type T, an instance should be declared once (and once only)
   */
  static std::map< int, std::vector< dataPtrType  > > store;
  static bool use_store; ///< Should the store be used?
  
  /*!
   * Returns a pointer to an VArrayData object with no
   * references. This is either from the store, or newly allocated
   */
  dataPtrType get(int len) {
    std::vector<dataPtrType >& st = store[len];
    if(st.empty()) {
      return std::make_shared<std::valarray<T>>(len);
    }
    dataPtrType p = st.back();
    st.pop_back();
    return p;
  }

  /*!
   * Release an VArrayData object, reducing its reference count by one. 
   * If no more references, then put back into the store.
   * Note: We have to pass a reference to the pointer, otherwise we get
   * a copy of the shared_ptr, which therefore increases the use count
   * and doesn't allow us to free the pass pointer directly
   */
  void release(dataPtrType &d) {
    if(!d) return;
    
    //If this is the only reference then put back into store or free
    if(d.use_count()==1) {
      if(use_store) {
        // Put back into store
	store[d->size()].push_back(std::move(d));
      }else {
	d.reset();
      }
    } else { //If not the last reference then just set it to nullptr
      d.reset();
    }
  }
  
};

/*!
 * Create a copy of an VArray, which does not share data
 */ 
template<typename T>
VArray<T>& copy(const VArray<T> &other) {
  VArray<T> a(other);
  a.ensureUnique();
  return a;
}

#endif // __VARRAY_H__

