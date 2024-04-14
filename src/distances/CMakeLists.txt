add_library(TravelingSalesmanSolver_distances)
target_sources(TravelingSalesmanSolver_distances PRIVATE
    distances.cpp
    distances_builder.cpp)
target_include_directories(TravelingSalesmanSolver_distances PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_distances PUBLIC
    OptimizationTools::utils)
add_library(TravelingSalesmanSolver::distances ALIAS TravelingSalesmanSolver_distances)
