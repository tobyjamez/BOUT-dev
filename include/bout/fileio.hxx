
#include "bout/deprecated.hxx"
#include "bout_types.hxx"
#include "dataformat.hxx"
#include "../src/fileio/formatfactory.hxx"
#include <boutexception.hxx>
#include <output.hxx>
#include <boutcomm.hxx>

#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////
// A dataset represents a value saved to file, for example a field
// or vector. It should handle attributes as well as values.

// NOTE: Need to ensure that setting the value for a time-dependent value
// multiple times does not result in time indices getting out of sync.

class Dataset {
public:
  Dataset(DataFormat *file, std::string name, bool time_dependent)
    : file(file), name(name), time_dependent(time_dependent) {
    size = file->getSize(name);
  }
  
  /// Cast operator, allowing assignment to different variable types
  /// e.g Field3D f = dump["value"];
  template <typename T> operator T() { return get<T>(); }
  
  // Get the value as a specified type
  // Some way to compare the requested type against the type
  // in the file is needed.
  // This method can be specialised for particlular types e.g. Field3D
  template <typename T> T get() {
    
    T result;

    if (time_dependent) {
      file->read_rec(&result, name);
    } else {
      file->read(&result, name);
    }
    return result;
  }

  // Assign a value
  template <typename T> T operator=(T value) {
    if (time_dependent) {
      file->write_rec(&value, name);
    } else {
      file->write(&value, name);
    }
    
    return value;
  }
private:
  DataFormat *file;  ///< Shared pointer?
  std::string name;  ///< Name of the variable
  bool time_dependent;
  std::vector<int> size;
};

///////////////////////////////////////////////////////////////////////

class File {
public:
  static const int ReadOnly = 0;
  static const int ReadWrite = 1;
  static const int Create = 2;
  
  File(const std::string name, int flags = File::ReadOnly) {
    file = FormatFactory::getInstance()->createDataFormat(name.c_str(), false);
    if (!file)
      throw BoutException("Datafile::open: Factory failed to create a DataFormat!");

    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);

    if (flags == ReadOnly) {
      if(!file->openr(name, MYPE))
        throw BoutException("Datafile::open: Failed to open file!");
    } else {
      bool appending = !(flags & Create);
      if(!file->openw(name, MYPE, appending))
        throw BoutException("Datafile::open: Failed to open file!");
    }
  }
  
  ~File() {
    file->close();
  }
  
  Dataset operator[](const string &name) {
    return Dataset(file.get(), name, false);
  }

  template <typename T>
  Dataset createDataset(const string &name) {
    throw BoutException("Can't create this dataset type");
  }
  
private:
  std::unique_ptr<DataFormat> file; // Shared?
};

// Specialised templates
template <>
Dataset File::createDataset<int>(const string &name) {
  if (!file->addVarInt(name, false)) {
    throw BoutException("Failed to add int variable %s to Datafile", name.c_str());
  }
  return Dataset(file.get(), name, false);
}
