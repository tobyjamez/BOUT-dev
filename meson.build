project('BOUT++', 'cpp',
        version: '4.1.1',
        default_options : ['cpp_std=c++11'])

configuration_object = configuration_data()

fftw3 = dependency('fftw3')
netcdf = dependency('netcdf', required: false)
netcdf_cpp = dependency('netcdf-cxx4', required: false)
hdf5 = dependency('hdf5', required: false)
mpi = dependency('mpi', language: 'cpp')

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

deps = [fftw3, netcdf, netcdf_cpp, hdf5, mpi]


incdir = include_directories('include')

configuration_inc = include_directories('.')

configure_file(input: 'boutconfig.hxx.in',
               output: 'boutconfig.hxx',
               configuration: configuration_object)

src = [
  './src/invert/laplacexy/laplacexy.cxx',
  './src/invert/laplace3d/laplace3d.cxx',
  './src/invert/laplace3d/laplace3d_factory.cxx',
  './src/invert/laplace3d/impls/petsc/petsc_laplace3d.cxx',
  './src/invert/laplacexz/impls/petsc/laplacexz-petsc.cxx',
  './src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.cxx',
  './src/invert/laplacexz/laplacexz.cxx',
  './src/invert/laplace/impls/cyclic/cyclic_laplace.cxx',
  './src/invert/laplace/impls/serial_band/serial_band.cxx',
  './src/invert/laplace/impls/mumps/mumps_laplace.cxx',
  './src/invert/laplace/impls/serial_tri/serial_tri.cxx',
  './src/invert/laplace/impls/pdd/pdd.cxx',
  './src/invert/laplace/impls/shoot/shoot_laplace.cxx',
  './src/invert/laplace/impls/petsc/petsc_laplace.cxx',
  './src/invert/laplace/impls/spt/spt.cxx',
  './src/invert/laplace/impls/multigrid/multigrid_alg.cxx',
  './src/invert/laplace/impls/multigrid/multigrid_laplace.cxx',
  './src/invert/laplace/impls/multigrid/multigrid_solver.cxx',
  './src/invert/laplace/laplacefactory.cxx',
  './src/invert/laplace/invert_laplace.cxx',
  './src/invert/parderiv/impls/cyclic/cyclic.cxx',
  './src/invert/parderiv/impls/serial/serial.cxx',
  './src/invert/parderiv/invert_parderiv.cxx',
  './src/invert/parderiv/parderiv_factory.cxx',
  './src/invert/fft_fftw.cxx',
  './src/invert/lapack_routines.cxx',
  './src/mesh/interpolation.cxx',
  './src/mesh/mesh.cxx',
  './src/mesh/boundary_region.cxx',
  './src/mesh/boundary_factory.cxx',
  './src/mesh/meshfactory.cxx',
  './src/mesh/surfaceiter.cxx',
  './src/mesh/coordinates.cxx',
  './src/mesh/parallel_boundary_region.cxx',
  './src/mesh/impls/bout/boutmesh.cxx',
  './src/mesh/interpolation/bilinear.cxx',
  './src/mesh/interpolation/hermite_spline.cxx',
  './src/mesh/interpolation/interpolation_factory.cxx',
  './src/mesh/interpolation/lagrange_4pt.cxx',
  './src/mesh/data/gridfromfile.cxx',
  './src/mesh/data/gridfromoptions.cxx',
  './src/mesh/parallel/fci.cxx',
  './src/mesh/parallel/shiftedmetric.cxx',
  './src/mesh/fv_ops.cxx',
  './src/mesh/boundary_standard.cxx',
  './src/mesh/parallel_boundary_op.cxx',
  './src/mesh/index_derivs.cxx',
  './src/mesh/difops.cxx',
  './src/physics/gyro_average.cxx',
  './src/physics/physicsmodel.cxx',
  './src/physics/smoothing.cxx',
  './src/physics/sourcex.cxx',
  './src/sys/msg_stack.cxx',
  './src/sys/output.cxx',
  './src/sys/range.cxx',
  './src/sys/utils.cxx',
  './src/sys/boutcomm.cxx',
  './src/sys/petsclib.cxx',
  './src/sys/timer.cxx',
  './src/sys/boutexception.cxx',
  './src/sys/options/options_ini.cxx',
  './src/sys/options.cxx',
  './src/sys/expressionparser.cxx',
  './src/sys/comm_group.cxx',
  './src/sys/optionsreader.cxx',
  './src/sys/derivs.cxx',
  './src/sys/slepclib.cxx',
  './src/sys/stencils.cxx',
  './src/field/fieldperp.cxx',
  './src/field/initialprofiles.cxx',
  './src/field/field2d.cxx',
  './src/field/fieldgroup.cxx',
  './src/field/vecops.cxx',
  './src/field/where.cxx',
  './src/field/field.cxx',
  './src/field/field_data.cxx',
  './src/field/fieldgenerators.cxx',
  './src/field/globalfield.cxx',
  './src/field/vector3d.cxx',
  './src/field/field_factory.cxx',
  './src/field/field3d.cxx',
  './src/field/vector2d.cxx',
  './src/solver/impls/rkgeneric/impls/rkf34/rkf34.cxx',
  './src/solver/impls/rkgeneric/impls/rkf45/rkf45.cxx',
  './src/solver/impls/rkgeneric/impls/cashkarp/cashkarp.cxx',
  './src/solver/impls/rkgeneric/impls/rk4simple/rk4simple.cxx',
  './src/solver/impls/rkgeneric/rkschemefactory.cxx',
  './src/solver/impls/rkgeneric/rkgeneric.cxx',
  './src/solver/impls/rkgeneric/rkscheme.cxx',
  './src/solver/impls/karniadakis/karniadakis.cxx',
  './src/solver/impls/petsc/petsc.cxx',
  './src/solver/impls/rk4/rk4.cxx',
  './src/solver/impls/arkode/arkode.cxx',
  './src/solver/impls/imex-bdf2/imex-bdf2.cxx',
  './src/solver/impls/cvode/cvode.cxx',
  './src/solver/impls/petsc-3.1/petsc-3.1.cxx',
  './src/solver/impls/power/power.cxx',
  './src/solver/impls/petsc-3.4/petsc-3.4.cxx',
  './src/solver/impls/slepc-3.4/slepc-3.4.cxx',
  './src/solver/impls/euler/euler.cxx',
  './src/solver/impls/petsc-3.2/petsc-3.2.cxx',
  './src/solver/impls/pvode/pvode.cxx',
  './src/solver/impls/ida/ida.cxx',
  './src/solver/impls/petsc-3.3/petsc-3.3.cxx',
  './src/solver/impls/rk3-ssp/rk3-ssp.cxx',
  './src/solver/impls/petsc-3.5/petsc-3.5.cxx',
  './src/solver/impls/snes/snes.cxx',
  './src/solver/solverfactory.cxx',
  './src/solver/monitor.cxx',
  './src/solver/solver.cxx',
  './src/fileio/impls/pnetcdf/pnetcdf.cxx',
  './src/fileio/impls/netcdf/nc_format.cxx',
  './src/fileio/impls/netcdf4/ncxx4.cxx',
  './src/fileio/impls/hdf5/h5_format.cxx',
  './src/fileio/dataformat.cxx',
  './src/fileio/formatfactory.cxx',
  './src/fileio/datafile.cxx',
  './src/bout++.cxx'
]

libbout = static_library('bout++', sources: src,
                         include_directories: [incdir, configuration_inc],
                         dependencies: deps,
                         cpp_args: '-fopenmp', link_args : '-fopenmp')

serial_test_sources = [
  './tests/unit/field/test_field2d.cxx',
  './tests/unit/field/test_field3d.cxx',
  './tests/unit/field/test_fieldgroup.cxx',
  './tests/unit/field/test_vector2d.cxx',
  './tests/unit/field/test_vector3d.cxx',
  './tests/unit/sys/test_output.cxx',
  './tests/unit/sys/test_optionsreader.cxx',
  './tests/unit/sys/test_range.cxx',
  './tests/unit/sys/test_utils.cxx',
  './tests/unit/sys/test_options.cxx',
  './tests/unit/sys/test_array.cxx',
  './tests/unit/sys/test_boutexception.cxx',
  './tests/unit/sys/test_expressionparser.cxx',
  './tests/unit/sys/test_msg_stack.cxx',
  './tests/unit/test_extras.cxx',
  './tests/unit/bout_test_main.cxx',
]

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

test_deps = deps
test_deps += gtest_dep
test_incdir = include_directories('tests/unit')

serial_test = executable('serial_test',
                         sources: serial_test_sources,
                         include_directories: [incdir, configuration_inc, test_incdir],
                         dependencies: test_deps,
                         cpp_args: '-fopenmp', link_args : '-fopenmp',
                         link_with: libbout)

test('serial test', serial_test)
