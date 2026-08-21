if(NOT DEFINED FEMASTER_EXECUTABLE OR NOT EXISTS "${FEMASTER_EXECUTABLE}")
    message(FATAL_ERROR "FEMASTER_EXECUTABLE must name an existing FEMaster binary")
endif()
if(NOT DEFINED FEMASTER_SOURCE_DIR)
    message(FATAL_ERROR "FEMASTER_SOURCE_DIR is required")
endif()
if(NOT DEFINED FEMASTER_LOADSTEP_DIR)
    set(FEMASTER_LOADSTEP_DIR "${CMAKE_CURRENT_BINARY_DIR}/pretension_loadsteps")
endif()

file(REMOVE_RECURSE "${FEMASTER_LOADSTEP_DIR}")
file(MAKE_DIRECTORY "${FEMASTER_LOADSTEP_DIR}")

function(run_pretension_loadstep_case case_name source_deck)
    set(input_file "${FEMASTER_LOADSTEP_DIR}/${case_name}.inp")
    configure_file("${source_deck}" "${input_file}" COPYONLY)
    execute_process(
        COMMAND "${FEMASTER_EXECUTABLE}" "${input_file}"
        WORKING_DIRECTORY "${FEMASTER_SOURCE_DIR}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE stdout
        ERROR_VARIABLE stderr
    )
    set(log "${stdout}\n${stderr}")
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "${case_name}: FEMaster exited with ${result}\n${log}")
    endif()
    string(FIND "${log}" "[FAIL]" failure_position)
    if(NOT failure_position EQUAL -1)
        message(FATAL_ERROR "${case_name}: reported a failed post-check\n${log}")
    endif()
    string(REGEX MATCHALL "\\[PASS\\] lagrange constraints" constraint_passes "${log}")
    list(LENGTH constraint_passes pass_count)
    if(pass_count LESS 2)
        message(FATAL_ERROR "${case_name}: expected repeated constraint passes\n${log}")
    endif()
    set(${case_name}_log "${log}" PARENT_SCOPE)
    message(STATUS "Pretension loadsteps PASS: ${case_name}")
endfunction()

run_pretension_loadstep_case(
    displacement_steps
    "${FEMASTER_SOURCE_DIR}/examples/44_pretension_loadsteps/44_pretension_loadsteps.inp")
foreach(required
        "solved mean gap = 0.02"
        "solved mean gap = 0.05"
        "locked gap = 0.05")
    string(FIND "${displacement_steps_log}" "${required}" match_position)
    if(match_position EQUAL -1)
        message(FATAL_ERROR "displacement_steps: missing '${required}'")
    endif()
endforeach()

run_pretension_loadstep_case(
    force_lock
    "${FEMASTER_SOURCE_DIR}/examples/45_pretension_force_lock/45_pretension_force_lock.inp")
string(REGEX MATCH "solved mean gap = ([^\r\n]+)" solved_match "${force_lock_log}")
if(NOT solved_match)
    message(FATAL_ERROR "force_lock: solved gap was not reported")
endif()
set(solved_gap "${CMAKE_MATCH_1}")
string(FIND "${force_lock_log}" "locked gap = ${solved_gap}" lock_position)
if(lock_position EQUAL -1)
    message(FATAL_ERROR "force_lock: solved gap ${solved_gap} was not transferred to LOCK")
endif()

file(REMOVE_RECURSE "${FEMASTER_LOADSTEP_DIR}")
message(STATUS "Pretension loadstep regression completed")
