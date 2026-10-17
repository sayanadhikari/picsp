set(PICSP_TEST_NAMES test_poisson test_particles test_geometry test_physics test_magnetic test_collisions)

foreach(t ${PICSP_TEST_NAMES})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE picsp_core)
  add_test(NAME ${t}_serial COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${t}>)
  add_test(NAME ${t}_mpi4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${t}>)
  set_tests_properties(${t}_serial ${t}_mpi4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
endforeach()
