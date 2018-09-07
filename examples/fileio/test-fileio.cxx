#include "bout/fileio.hxx"

#include "bout.hxx"

int main(int argc, char* argv[]) {
  BoutInitialise(argc, argv);

  {
    // Open a file, specifying mode
    File dump("test.nc", File::ReadWrite | File::Create);
    
    // Create a dataset if it doesn't exist.
    // If it does exist, check if it's the same as what already exists
    
    auto dataset = dump.createDataset<int>("val");
    
    // Perhaps also allow creating a dataset by giving an example value:
    //   Field3D val;
    //   auto dataset = dump.createDataset("val", val);

    dataset = 4; // Give it a value
    
  } // File goes out of scope, should close.
  
  {
    File dump("test.nc", File::ReadWrite);

    int value = dump["val"]; // Read. This uses the cast operator
    
    output << "Value is " << value << "\n";
    
    dump["val"] = 3; // Change value (write). 

    output << "New value " << dump["val"].get<int>() << "\n";
  }

  // Time dependent variables may be more tricky, ensuring that all
  // fields have the correct time axis.

  // Perhaps time axes are a separate object:
  //
  //   auto time_axis = dump.createTimeAxis["t"];
  //   auto dataset = dump.createDataset<int>("val", time_axis);
  //
  // Now set the time index which all writes will use:
  //
  //   dump["val"] = 3;
  //   dump.timeAxis("t").increment();
  //   dump["val"] = 4;
  //
  //   val is now [3,4]
  //
  // This now has some non-local effects, since the time axis is shared
  // between several variables. On the other hand this is the point, to
  // keep variables syncronised.
  
  BoutFinalise();
  return 0;
}
