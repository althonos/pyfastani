file(COPY ../vendor/boost-math DESTINATION .)
file(COPY ../vendor/cpu_features DESTINATION .)
file(COPY ../vendor/FastANI DESTINATION .)
file(COPY ./FastANI DESTINATION .)

if(NOT HAVE_BOOST)
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/boost-math ${CMAKE_CURRENT_BINARY_DIR}/boost-math/build)
endif()

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/cpu_features ${CMAKE_CURRENT_BINARY_DIR}/cpu_features/build)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/FastANI ${CMAKE_CURRENT_BINARY_DIR}/FastANI/build)
add_subdirectory(pyfastani)