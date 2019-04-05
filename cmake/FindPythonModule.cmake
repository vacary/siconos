# Search for a Python module.
#
# Usage :
#
#  find_python_module(mpi4py REQUIRED)
#  find_python_module(sphinx)
#
#  Warning : use ${PYTHON_EXECUTABLE} as python interpreter
# 


include(FindPackageHandleStandardArgs)
function(find_python_module module)
  
  if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
    set(${module}_FIND_REQUIRED TRUE)
  endif()

  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import ${module} as name; print(name.__file__)"
    RESULT_VARIABLE ${module}_FIND_RESULT     # Return code from command above
    OUTPUT_VARIABLE ${module}_FIND_OUTPUT     # Standard output form command above
    ERROR_QUIET # Ignores quietly standard error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  
  if(NOT ${module}_FIND_RESULT) # Return code == 0 means that things have gone well
    set(${module}_file ${${module}_FIND_OUTPUT} CACHE STRING "Python ${module} module file.")
  endif()

  find_package_handle_standard_args(${module} REQUIRED_VARS ${module}_file)

endfunction(find_python_module)
