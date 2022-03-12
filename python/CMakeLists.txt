set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6 3.7 3.8 3.9)
find_package(PythonInterp 3)

execute_process (
  COMMAND git describe --abbrev=4 --dirty --always --tags
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# file(WRITE  ${CMAKE_CURRENT_SOURCE_DIR}/symd/version.py "__version__ = ${GIT_COMMIT_HASH}")

if(NOT PYTHONINTERP_FOUND)
    message("Could not find Python interpreter. Not installing python package")
else()

    # detect virtualenv and set Pip args accordingly
    if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
        set(_pip_args)
    else()
        set(_pip_args "--user")
    endif()



    # install PyPI Python package using pip
    execute_process(
            COMMAND ${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} ${_pip_args})
endif()