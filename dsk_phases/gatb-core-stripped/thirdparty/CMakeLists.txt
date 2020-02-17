# include other smaller libraries (json, Boophf)

add_custom_target (thirdparty_copy ALL)
add_custom_command (TARGET thirdparty_copy COMMAND ${CMAKE_COMMAND} -E
    make_directory ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/json)
add_custom_command (TARGET thirdparty_copy  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/thirdparty/json/json.hpp   ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/json/)

################################################################################
#  INSTALL 
################################################################################
IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    INSTALL (DIRECTORY ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/json DESTINATION include)
ENDIF ()
