add_library(TravelingSalesmanSolver_candidates)
target_sources(TravelingSalesmanSolver_candidates PRIVATE
    nearest_neighbor.cpp
    alpha_nearness.cpp)
target_include_directories(TravelingSalesmanSolver_candidates PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_candidates PUBLIC
    TravelingSalesmanSolver_distances)
add_library(TravelingSalesmanSolver::candidates ALIAS TravelingSalesmanSolver_candidates)
