project('BOUT++', 'cpp',
        version: '4.1.2',
        default_options : ['cpp_std=c++11'])

configuration_object = configuration_data()

fftw3 = dependency('fftw3')
netcdf = dependency('netcdf', required: false)
netcdf_cpp = dependency('netcdf-cxx4', required: false)
hdf5 = dependency('hdf5', required: false)
mpi = dependency('mpi', language: 'cpp')

with_petsc = get_option('with-petsc')
petsc_enabled = false
if with_petsc
  petsc = dependency('PETSc', required: with_petsc)
  petsc_enabled = petsc.found()
endif

if not netcdf_cpp.found()
  ncxx4_config = find_program('ncxx4-config')
  if ncxx4_config.found()
    message('Using ncxx4-config')
    netcdf_cpp_libs = run_command(ncxx4_config, '--libs').stdout().strip().split()
    netcdf_cpp_inc = run_command(ncxx4_config, '--cflags').stdout().strip().split()
    netcdf_cpp = declare_dependency(link_args: netcdf_cpp_libs,
                                    compile_args: netcdf_cpp_inc)
  endif
endif

configuration_object.set('NCDF', netcdf.found() and (not netcdf_cpp.found()))
configuration_object.set('NCDF4', netcdf_cpp.found())
configuration_object.set('HDF5', hdf5.found())
configuration_object.set('BOUT_HAS_PETSC', petsc_enabled)
configuration_object.set('BOUT_HAS_PETSC_DEV', petsc_enabled)

deps = [fftw3, netcdf, netcdf_cpp, hdf5, mpi]

if petsc_enabled
  deps += petsc
endif

incdir = include_directories('include')

configuration_inc = include_directories('.')

configure_file(input: 'boutconfig.hxx.in',
               output: 'boutconfig.hxx',
               configuration: configuration_object)

gtest_dir = './googletest/googletest'
gtest_incdir = include_directories(join_paths(gtest_dir, 'include'), is_system : true)

libgtest = static_library('gtest',
                           cpp_args : ['-w'],
                           include_directories : [include_directories(gtest_dir), gtest_incdir],
                           sources : [
                             join_paths(gtest_dir, 'src', 'gtest-all.cc'),
                             join_paths(gtest_dir, 'src', 'gtest_main.cc')
                           ])

gtest_dep = declare_dependency(dependencies : dependency('threads'),
                               include_directories : gtest_incdir,
                               link_with : libgtest)

subdir('src')
subdir('tests/unit')
