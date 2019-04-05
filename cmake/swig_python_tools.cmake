# =======================================
# Macros and functions related to swig
# and python
# =======================================

# -----------------------------------
# Build targets to generate python
# docstrings from xml doxygen output.
#
# 1) Doxygen to create xml outputs.
# 2) Doxy2swig to
#    create docstrings from xml outputs.
#
# This macro must be called for each
# component, in LibraryProjectSetup.
#
# -- Targets :
# 1) headers --> xml using doxygen
# 2) xml files --> .i using doxy2swig (one .i for each xml file)
# 3) .i files --> ${COMP}-docstrings.i (one single file)

# 1 : target ${COMP}_xml4swig
# 2 and 3 : target ${COMP}_docstrings

# Note FP : the whole process must be re-executed for any change in a header file of the component
# (unless we set a doxygen conf for each header, don't think we need to bother with that.)
# 
# Doxygen steps are driven by cmake while doxy2swig and related are hidden in build_docstrings python
# script.
# -----------------------------------
macro(doxy2swig_docstrings COMP)
  if(WITH_${COMP}_DOXY2SWIG)
    
    # -- doxygen/xml config --

    # 1 - Create COMPdoxy.config.xml in binary dir,
    #     from doxy.config (source dir), taking
    #     into account sources/headers files
    #     of current component COMP.
    #
    # 2 - Run doxygen to build xml documentation
    #
    # Results in CMAKE_BINARY_DIR/docs/build/xml-docstrings
    #
    # Warning : output path must be different from
    # output for xml used by sphinx/breathe/exhale.
    # -----------------------------------------------------
    
    # Path where xml files will be generated by doxygen as inputs for swig.
    set(DOXY2SWIG_OUTPUT ${DOXYGEN_OUTPUT}/doxy2swig-xml/${COMP})
    file(MAKE_DIRECTORY ${DOXY2SWIG_OUTPUT})

    # Set doxygen config file name
    set(DOXY_CONFIG_XML "${CMAKE_BINARY_DIR}/docs/config/${COMP}doxy2swig-xml.config")
    # Create doxygen configuration file in binary dir
    set(DOXYGEN_INPUTS ${${COMP}_DOXYGEN_INPUTS}) #  list of inputs
    set(DOXY_QUIET "YES")
    set(DOXY_WARNINGS "NO")
    set(GENERATE_HTML NO)
    set(GENERATE_XML YES)
    set(EXTRACT_ALL NO)
    set(EXTRACT_PRIVATE NO)
    set(XML_OUTPUT doxy2swig-xml/${COMP})
    message(" -- Create doxygen conf (xml for docstrings) for component ${COMP}")
    configure_file(${CMAKE_SOURCE_DIR}/docs/config/doxy2swig.config.in ${DOXY_CONFIG_XML} @ONLY)

    # -- target to generate xml files (doxygen) from current component source files --
    #
    add_custom_target(${COMP}_xml4swig
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG_XML}
      OUTPUT_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy.log ERROR_FILE ${DOXYGEN_OUTPUT}/${COMP}doxy.log
      COMMENT " -- Build xml (for swig) doc for component ${COMP} ..."
      )
    
    # -- command to build .i files from xml (doxygen outputs) for current component  --
    # Call a python function from gendoctools (build_docstrings).
    add_custom_command(OUTPUT  ${SICONOS_SWIG_ROOT_DIR}/${COMP}-docstrings.i
      DEPENDS ${COMP}_xml4swig
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
      "from gendoctools.python2rst import build_docstrings as f; f('${${COMP}_HDRS}', '${COMP}', '${DOXY_CONFIG_XML}', '${SICONOS_SWIG_ROOT_DIR}')"
      VERBATIM
      )
    add_custom_target(${COMP}_docstrings DEPENDS ${SICONOS_SWIG_ROOT_DIR}/${COMP}-docstrings.i
      COMMENT "Create swig files from xml for component ${COMP}.")
  else()
    # No doxy2swig but 
    # generate empty ${COMP}-docstrings.i file (required because of %include in swig files)
    add_custom_command(OUTPUT ${SICONOS_SWIG_ROOT_DIR}/${COMP}-docstrings.i
      COMMAND touch
      ARGS ${SICONOS_SWIG_ROOT_DIR}/${COMP}-docstrings.i
      )
    add_custom_target(${COMP}_docstrings DEPENDS ${SICONOS_SWIG_ROOT_DIR}/${COMP}-docstrings.i
      COMMENT "This target has no real effect since swig/docstrings is not activated (PYTHON_WRAPPER and/or DOXY2SWIG is OFF).")

  endif()
  

endmacro()

# ----------------------------------------------------------------------
# Build a swig module from .i files
#
# Usage :
# add_siconos_swig_submodule(full_name)
#
# full_name must be a path to file.i, path relative to siconos
# python package root.
# For example, to build a swig module 'bullet' in siconos.mechanics.contact_detection,
# call this macro with full_name = mechanics/contact_detection/bullet
# and to build a swig module for kernel in siconos,
# call this macro with full_name = ./kernel
#
# ----------------------------------------------------------------------
macro(add_siconos_swig_sub_module fullname)

  get_filename_component(_name ${fullname} NAME)
  get_filename_component(_path ${fullname} PATH)
  message(" -- Build module ${_name} in directory ${_path} for parent ${COMPONENT}")
  # Add component dependencies to the current submodule deps.
  if(DEFINED SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)
    set(SWIG_MODULE_${_name}_EXTRA_DEPS ${SWIG_MODULE_${COMPONENT}_EXTRA_DEPS})
  endif()

  if(WITH_DOXY2SWIG)
    list(APPEND SWIG_MODULE_${_name}_EXTRA_DEPS ${COMPONENTS}_docstrings)
  endif()
  
  # add as dependencies all the i files
  file(GLOB ${_name}_I_FILES ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/*.i)
  foreach(_f IN LISTS ${_name}_I_FILES)
    list(APPEND SWIG_MODULE_${_name}_EXTRA_DEPS ${_f})
  endforeach()

  # set properties for current '.i' file
  set(swig_file ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/${_name}.i)

  # set output dir
  set(CMAKE_SWIG_OUTDIR "${SICONOS_SWIG_ROOT_DIR}/${_path}")

  # compile flags
  foreach(_dir IN LISTS ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)
    set(${COMPONENT}_SWIG_DEFS "-I${_dir};${${COMPONENT}_SWIG_DEFS}")
  endforeach()

  # extra per-module flags if any
  set(${COMPONENT}_SWIG_DEFS_${_name} "${${COMPONENT}_SWIG_DEFS};${${COMPONENT}_SWIG_DEFS_${_name}}")

  IF(WITH_CXX AND (BUILD_AS_CPP OR NOT ${COMPONENT} MATCHES "numerics"))
    set_source_files_properties(${swig_file}
      PROPERTIES SWIG_FLAGS "${${COMPONENT}_SWIG_DEFS_${_name}}" CPLUSPLUS ON)
  ELSE(WITH_CXX AND (BUILD_AS_CPP OR NOT ${COMPONENT} MATCHES "numerics"))
    # C compilation, pass SWIG_FLAGS.
    IF(${COMPONENT} MATCHES "numerics")
      set_source_files_properties(${swig_file}
        PROPERTIES SWIG_FLAGS "${${COMPONENT}_SWIG_DEFS_${_name}}")
    ENDIF()
  ENDIF()

  # --- build swig module ---
  if(CMAKE_VERSION VERSION_LESS 3.8.0)
    swig_add_module(${_name} python ${swig_file})
  else()
    set(ADDITIONAL_SWIG_DEFINES ${ADDITIONAL_SWIG_DEFINES} -DBOOST_NOEXCEPT)
    swig_add_library(${_name} LANGUAGE python SOURCES ${swig_file})
  endif()
  # WARNING ${swig_generated_file_fullname} is overriden 
  set(${_name}_generated_file_fullname ${swig_generated_file_fullname})
  set_source_files_properties( ${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing")
  # Set path for the library generated by swig for the current module --> siconos python package path
  set_property(TARGET ${SWIG_MODULE_${_name}_REAL_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/${_path})
  message(" -- ${_name} generated (swig) file will be ${${_name}_generated_file_fullname}")

  # Set the SONAME for the SWIG module to the Siconos SONAME
  set_target_properties(${SWIG_MODULE_${name}_REAL_NAME} PROPERTIES
    NO_SONAME OFF
    VERSION "${SICONOS_SOVERSION}"
    SOVERSION "${SICONOS_SOVERSION_MAJOR}")

  IF(MSVC AND ${COMPONENT} MATCHES "kernel")
    set_source_files_properties(${${_name}_generated_file_fullname} PROPERTIES COMPILE_FLAGS "/bigobj")
  ENDIF()

  # Add a post-build step that prepends utf-8 coding indicator to .py files
  add_custom_command(TARGET ${SWIG_MODULE_${_name}_REAL_NAME}
    POST_BUILD COMMAND sh -c "(echo '# -*- coding: utf-8 -*-'; cat ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.py) > ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.tmp; mv ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.tmp ${SICONOS_SWIG_ROOT_DIR}/${_path}/${_name}.py" VERBATIM)

  # Check dependencies and then link ...
  add_dependencies(${SWIG_MODULE_${_name}_REAL_NAME} ${COMPONENT})

  if(UNIX AND NOT APPLE)
    # do not link against the Python library on unix, it is useless
    swig_link_libraries(${_name} ${${COMPONENT}_LINK_LIBRARIES} ${COMPONENT})
  else()
    swig_link_libraries(${_name} ${Python3_LIBRARIES} ${${COMPONENT}_LINK_LIBRARIES} ${COMPONENT})
  endif()

  # set dependency of sphinx apidoc to this target
  if(WITH_DOCUMENTATION AND WITH_${COMPONENT}_DOXY2SWIG)
    include(doc_tools)
    docstrings2rst(${_path} ${_name})
  endif()
  
  # Copy __init__.py file if needed
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py.in)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py.in
      ${SICONOS_SWIG_ROOT_DIR}/${_path}/__init__.py @ONLY)
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_path}/__init__.py
      ${SICONOS_SWIG_ROOT_DIR}/${_path}/__init__.py COPYONLY)
  endif()

  # --- install python files and target ---
  # install path ...
  set(DEST "${SICONOS_PYTHON_INSTALL_DIR}/${SICONOS_PYTHON_PACKAGE}/${_path}")

  #install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${_name}.py DESTINATION ${DEST})
  install(TARGETS ${SWIG_MODULE_${_name}_REAL_NAME} LIBRARY DESTINATION ${DEST})
  
endmacro()

macro(swig_module_setup modules_list)
  if(WITH_${COMPONENT}_PYTHON_WRAPPER)
    # we should use lowercase name for python module (pep8 ...)
    message(" -- Prepare python bindings for component ${COMPONENT} ...")
    # build python bindings
    include_directories(${SICONOS_SWIG_INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)
    if(HAVE_SICONOS_IO)
      include_directories(${CMAKE_SOURCE_DIR}/io/src)
      if(HAVE_SICONOS_MECHANICS)
	include_directories(${CMAKE_SOURCE_DIR}/io/src/mechanics)
      endif()
    endif()
    foreach(module IN LISTS ${COMPONENT}_PYTHON_MODULES)
      add_siconos_swig_sub_module(${module})
    endforeach()
  endif()
endmacro()

include(tools4tests)
