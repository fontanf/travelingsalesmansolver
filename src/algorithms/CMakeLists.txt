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

add_library(TravelingSalesmanSolver_random_permutation)
target_sources(TravelingSalesmanSolver_random_permutation PRIVATE
    random_permutation.cpp)
target_include_directories(TravelingSalesmanSolver_random_permutation PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_random_permutation PUBLIC
    TravelingSalesmanSolver_traveling_salesman)
add_library(TravelingSalesmanSolver::random_permutation ALIAS TravelingSalesmanSolver_random_permutation)

add_library(TravelingSalesmanSolver_random_walk)
target_sources(TravelingSalesmanSolver_random_walk PRIVATE
    random_walk.cpp)
target_include_directories(TravelingSalesmanSolver_random_walk PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_random_walk PUBLIC
    TravelingSalesmanSolver_traveling_salesman
    TravelingSalesmanSolver_candidates)
add_library(TravelingSalesmanSolver::random_walk ALIAS TravelingSalesmanSolver_random_walk)

add_library(TravelingSalesmanSolver_two_opt)
target_sources(TravelingSalesmanSolver_two_opt PRIVATE
    two_opt.cpp)
target_include_directories(TravelingSalesmanSolver_two_opt PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_two_opt PUBLIC
    TravelingSalesmanSolver_traveling_salesman
    TravelingSalesmanSolver_candidates
    TravelingSalesmanSolver_random_permutation)
add_library(TravelingSalesmanSolver::two_opt ALIAS TravelingSalesmanSolver_two_opt)

add_library(TravelingSalesmanSolver_greedy_edge)
target_sources(TravelingSalesmanSolver_greedy_edge PRIVATE
    greedy_edge.cpp)
target_include_directories(TravelingSalesmanSolver_greedy_edge PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_greedy_edge PUBLIC
    TravelingSalesmanSolver_traveling_salesman
    TravelingSalesmanSolver_candidates)
add_library(TravelingSalesmanSolver::greedy_edge ALIAS TravelingSalesmanSolver_greedy_edge)

add_library(TravelingSalesmanSolver_lin_kernighan)
target_sources(TravelingSalesmanSolver_lin_kernighan PRIVATE
    lin_kernighan.cpp)
target_include_directories(TravelingSalesmanSolver_lin_kernighan PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_lin_kernighan PUBLIC
    TravelingSalesmanSolver_traveling_salesman
    TravelingSalesmanSolver_greedy_edge)
add_library(TravelingSalesmanSolver::lin_kernighan ALIAS TravelingSalesmanSolver_lin_kernighan)

add_library(TravelingSalesmanSolver_eax)
target_sources(TravelingSalesmanSolver_eax PRIVATE
    eax.cpp)
target_include_directories(TravelingSalesmanSolver_eax PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(TravelingSalesmanSolver_eax PUBLIC
    TravelingSalesmanSolver_traveling_salesman
    TravelingSalesmanSolver_candidates
    TravelingSalesmanSolver_random_permutation
    TravelingSalesmanSolver_random_walk
    TravelingSalesmanSolver_two_opt
    TravelingSalesmanSolver_lin_kernighan)
add_library(TravelingSalesmanSolver::eax ALIAS TravelingSalesmanSolver_eax)
