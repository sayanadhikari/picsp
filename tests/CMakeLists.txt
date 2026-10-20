set(PICSP_TEST_NAMES test_poisson test_particles test_geometry test_physics test_magnetic test_collisions)

# --oversubscribe is an Open MPI flag; MPICH's Hydra launcher rejects it and oversubscribes by default
set(PICSP_OVERSUBSCRIBE "")
if(PICSP_MPIEXEC_IMPL STREQUAL "openmpi")
  set(PICSP_OVERSUBSCRIBE --oversubscribe)
endif()

foreach(t ${PICSP_TEST_NAMES})
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE picsp_core)
  # each variant runs in its own directory so the output files do not collide under ctest -j
  foreach(v serial mpi4)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/run_${t}_${v})
  endforeach()
  add_test(NAME ${t}_serial COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${t}>
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/run_${t}_serial)
  add_test(NAME ${t}_mpi4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${PICSP_OVERSUBSCRIBE} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${t}>
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/run_${t}_mpi4)
  set_tests_properties(${t}_serial PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
  set_tests_properties(${t}_mpi4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" PROCESSORS 4)
endforeach()
