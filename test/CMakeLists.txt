add_executable(test_ant test.cpp)
target_link_libraries(test_ant AntColonyBase)

add_test(test_p01 test_ant ${PROJECT_SOURCE_DIR}/data/p01.tsp)

set_tests_properties(test_p01 PROPERTIES PASS_REGULAR_EXPRESSION "1 13 2 15 9 5 7 3 12 14 10 8 6 4 11")