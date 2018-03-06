serial_test_sources = [
  './field/test_field2d.cxx',
  './field/test_field3d.cxx',
  './field/test_fieldgroup.cxx',
  './field/test_vector2d.cxx',
  './field/test_vector3d.cxx',
  './sys/test_output.cxx',
  './sys/test_optionsreader.cxx',
  './sys/test_range.cxx',
  './sys/test_utils.cxx',
  './sys/test_options.cxx',
  './sys/test_array.cxx',
  './sys/test_boutexception.cxx',
  './sys/test_expressionparser.cxx',
  './sys/test_msg_stack.cxx',
  './test_extras.cxx',
  './bout_test_main.cxx',
]

test_deps = deps
test_deps += gtest_dep
test_incdir = include_directories('.')

serial_test = executable('serial_test',
                         sources: serial_test_sources,
                         include_directories: [incdir, configuration_inc, test_incdir],
                         dependencies: test_deps,
                         cpp_args: '-fopenmp', link_args : '-fopenmp',
                         link_with: libbout)

test('serial test', serial_test)
