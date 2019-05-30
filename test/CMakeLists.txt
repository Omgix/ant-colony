add_executable(test_ant test.cpp)
target_link_libraries(test_ant AntColonyBase)

add_test(test_ulysses16 test_ant ${PROJECT_SOURCE_DIR}/data/ulysses16.tsp)

set_tests_properties(test_ulysses16 PROPERTIES PASS_REGULAR_EXPRESSION "1 14 13 12 7 6 15 5 11 9 10 16 3 2 4 8")