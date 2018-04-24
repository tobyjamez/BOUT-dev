#include <stdio.h>

#include "bout/array.hxx"
#include "gtest/gtest.h"

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from bout_test_main.cxx\n");
  testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  // Clean up the array store, so valgrind doesn't report false
  // positives
  Array<double>::cleanup();
  return result;
}
