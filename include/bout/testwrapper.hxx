/*!************************************************************************
 * \file physicsmodel.hxx
 *
 * @brief Base class for Physics Models
 *
 *
 *
 * Changelog:
 *
 * 2013-08 Ben Dudson <benjamin.dudson@york.ac.uk>
 *    * Initial version
 *
 **************************************************************************
 * Copyright 2013 B.D.Dudson
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __BOUT_TESTWRAPPER_H__
#define __BOUT_TESTWRAPPER_H__

#include "solver.hxx"
#include "unused.hxx"
#include <bout.hxx>
#include <msg_stack.hxx>
#include <options.hxx>

/*!
 * Macro to define a simple main() which runs
 * the given test. This should be sufficient
 * for most use cases, but a user can define their own
 * main() function if needed.
 *
 * Example
 * -------
 *
 * void test() {
 *   ..
 * };
 *
 * BOUTTEST(test);
 */
#define BOUTTEST(test)                                                                   \
  int main(int argc, char **argv) {                                                      \
    int init_err = BoutInitialise(argc, argv);                                           \
    if (init_err < 0)                                                                    \
      return 0;                                                                          \
    else if (init_err > 0)                                                               \
      return init_err;                                                                   \
    try {                                                                                \
      test();                                                                            \
    } catch (BoutException & e) {                                                        \
      output << "Error encountered\n";                                                   \
      output << e.what() << endl;                                                        \
      MPI_Abort(BoutComm::get(), 1);                                                     \
    }                                                                                    \
    BoutFinalise();                                                                      \
    return 0;                                                                            \
  }

#endif // __BOUT_TESTWRAPPER_H__
