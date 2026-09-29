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

add_library(TravelingSalesmanSolver_greedy)
target_sources(TravelingSalesmanSolver_greedy PRIVATE
    greedy.cpp)
target_include_directories(TravelingSalesmanSolver_greedy PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_greedy PUBLIC
    TravelingSalesmanSolver_traveling_salesman)
add_library(TravelingSalesmanSolver::greedy ALIAS TravelingSalesmanSolver_greedy)

add_library(TravelingSalesmanSolver_lin_kernighan)
target_sources(TravelingSalesmanSolver_lin_kernighan PRIVATE
    lin_kernighan.cpp)
target_include_directories(TravelingSalesmanSolver_lin_kernighan PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_lin_kernighan PUBLIC
    TravelingSalesmanSolver_traveling_salesman
    TravelingSalesmanSolver_greedy)
add_library(TravelingSalesmanSolver::lin_kernighan ALIAS TravelingSalesmanSolver_lin_kernighan)
