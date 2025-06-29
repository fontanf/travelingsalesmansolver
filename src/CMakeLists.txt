add_subdirectory(distances)

add_library(TravelingSalesmanSolver_traveling_salesman)
target_sources(TravelingSalesmanSolver_traveling_salesman PRIVATE
    instance.cpp
    solution.cpp
    algorithm_formatter.cpp)
target_include_directories(TravelingSalesmanSolver_traveling_salesman PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_traveling_salesman PUBLIC
    TravelingSalesmanSolver_distances
    OptimizationTools::utils
    OptimizationTools::containers)
add_library(TravelingSalesmanSolver::traveling_salesman ALIAS TravelingSalesmanSolver_traveling_salesman)

add_subdirectory(algorithms)

if(TRAVELINGSALESMANSOLVER_BUILD_MAIN)
  add_executable(TravelingSalesmanSolver_main)
  target_sources(TravelingSalesmanSolver_main PRIVATE
      main.cpp)
  target_link_libraries(TravelingSalesmanSolver_main PUBLIC
      TravelingSalesmanSolver_lkh
      TravelingSalesmanSolver_concorde
      Boost::program_options)
  set_target_properties(TravelingSalesmanSolver_main PROPERTIES OUTPUT_NAME "travelingsalesmansolver")
  install(TARGETS TravelingSalesmanSolver_main)
endif()
