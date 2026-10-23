add_executable(testMap testMap.C)
target_link_libraries(testMap ${PROJECT_NAME})
target_compile_options(testMap PUBLIC -Wall -Wextra -Wno-unused-function)
