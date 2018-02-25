
cmake_minimum_required (VERSION 3.4)
project (NPB-Callable)


add_executable(npb-ep ./exec/main_ep.c)
target_link_libraries(npb-ep m pthread)

add_executable(npb-is ./exec/main_is.c)
target_link_libraries(npb-is m pthread)

add_executable(npb-cg ./exec/main_cg.c)
target_link_libraries(npb-cg m pthread)

add_executable(npb-mg ./exec/main_mg.c)
target_link_libraries(npb-mg m pthread)

add_executable(npb-ft ./exec/main_ft.c)
target_link_libraries(npb-ft m pthread)