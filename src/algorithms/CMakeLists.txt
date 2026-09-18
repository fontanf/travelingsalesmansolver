add_library(TravelingSalesmanSolver_concorde)
target_sources(TravelingSalesmanSolver_concorde PRIVATE
    concorde.cpp)
target_include_directories(TravelingSalesmanSolver_concorde PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_concorde PUBLIC
    TravelingSalesmanSolver_traveling_salesman)
add_library(TravelingSalesmanSolver::concorde ALIAS TravelingSalesmanSolver_concorde)

add_library(TravelingSalesmanSolver_lkh)
target_sources(TravelingSalesmanSolver_lkh PRIVATE
    lkh.cpp)
target_include_directories(TravelingSalesmanSolver_lkh PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_lkh PUBLIC
    TravelingSalesmanSolver_traveling_salesman)
add_library(TravelingSalesmanSolver::lkh ALIAS TravelingSalesmanSolver_lkh)

add_library(TravelingSalesmanSolver_eax)
target_sources(TravelingSalesmanSolver_eax PRIVATE
    eax.cpp)
target_include_directories(TravelingSalesmanSolver_eax PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_eax PUBLIC
    TravelingSalesmanSolver_traveling_salesman)
add_library(TravelingSalesmanSolver::eax ALIAS TravelingSalesmanSolver_eax)
