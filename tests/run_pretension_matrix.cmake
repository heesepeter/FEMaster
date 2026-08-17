if(NOT DEFINED FEMASTER_EXECUTABLE OR NOT EXISTS "${FEMASTER_EXECUTABLE}")
    message(FATAL_ERROR "FEMASTER_EXECUTABLE must name an existing FEMaster binary")
endif()
if(NOT DEFINED FEMASTER_SOURCE_DIR)
    message(FATAL_ERROR "FEMASTER_SOURCE_DIR is required")
endif()
if(NOT DEFINED FEMASTER_MATRIX_DIR)
    set(FEMASTER_MATRIX_DIR "${CMAKE_CURRENT_BINARY_DIR}/pretension_matrix")
endif()

set(template_file
    "${FEMASTER_SOURCE_DIR}/examples/41_mixed_hex_block/41_mixed_hex_block.inp")
if(NOT EXISTS "${template_file}")
    message(FATAL_ERROR "Missing pretension matrix template: ${template_file}")
endif()

file(REMOVE_RECURSE "${FEMASTER_MATRIX_DIR}")
file(MAKE_DIRECTORY "${FEMASTER_MATRIX_DIR}")
file(READ "${template_file}" template_deck)

set(axis_names oblique_xyz oblique_x oblique_y exact_vertex near_axis)
set(axis_oblique_xyz "0.3, 0.2, 1.0")
set(axis_oblique_x   "1.0, 0.37, 0.23")
set(axis_oblique_y   "0.23, 1.0, 0.47")
set(axis_exact_vertex "1.0, 0.4, 0.2")
set(axis_near_axis   "0.013, 0.021, 1.0")

set(connectivity_1 "1, 2, 5, 4, 10, 11, 14, 13, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39")
set(connectivity_2 "2, 3, 6, 5, 11, 12, 15, 14")
set(connectivity_3 "4, 5, 8, 7, 13, 14, 17, 16")
set(connectivity_4 "5, 6, 9, 8, 14, 15, 18, 17")
set(connectivity_5 "10, 11, 14, 13, 19, 20, 23, 22")
set(connectivity_6 "11, 12, 15, 14, 20, 21, 24, 23")
set(connectivity_7 "13, 14, 17, 16, 22, 23, 26, 25")
set(connectivity_8 "14, 15, 18, 17, 23, 24, 27, 26, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51")

set(order_names original permuted)
set(ids_original "1;2;3;4;5;6;7;8")
set(ids_permuted "8;3;6;1;7;2;5;4")

set(case_count 0)
foreach(axis_name IN LISTS axis_names)
    set(axis "${axis_${axis_name}}")
    foreach(order_name IN LISTS order_names)
        set(deck "${template_deck}")
        string(REPLACE "0.3, 0.2, 1.0" "${axis}" deck "${deck}")
        set(ids "${ids_${order_name}}")

        foreach(old_id RANGE 1 8)
            math(EXPR list_index "${old_id} - 1")
            list(GET ids ${list_index} new_id)
            string(REPLACE
                "${old_id}, ${connectivity_${old_id}}"
                "${new_id}, ${connectivity_${old_id}}"
                deck "${deck}")
        endforeach()
        list(GET ids 0 surface_element_id)
        foreach(surface_id RANGE 1 6)
            string(REPLACE
                "${surface_id}, 1, ${surface_id}"
                "${surface_id}, ${surface_element_id}, ${surface_id}"
                deck "${deck}")
        endforeach()

        set(case_name "${axis_name}_${order_name}")
        set(input_file "${FEMASTER_MATRIX_DIR}/${case_name}.inp")
        file(WRITE "${input_file}" "${deck}")
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
        foreach(required
                "validated 12 shared C3D8 closure face(s)"
                "quality-gate failures = 0"
                "[PASS] lagrange equilibrium"
                "[PASS] lagrange constraints")
            string(FIND "${log}" "${required}" match_position)
            if(match_position EQUAL -1)
                message(FATAL_ERROR "${case_name}: missing '${required}'\n${log}")
            endif()
        endforeach()
        string(FIND "${log}" "[FAIL]" failure_position)
        if(NOT failure_position EQUAL -1)
            message(FATAL_ERROR "${case_name}: reported a failed post-check\n${log}")
        endif()

        file(REMOVE
            "${input_file}"
            "${FEMASTER_MATRIX_DIR}/${case_name}.res"
            "${FEMASTER_MATRIX_DIR}/${case_name}.frd")
        math(EXPR case_count "${case_count} + 1")
        message(STATUS "Pretension matrix PASS: ${case_name}")
    endforeach()
endforeach()

file(REMOVE_RECURSE "${FEMASTER_MATRIX_DIR}")
message(STATUS "Pretension matrix completed: ${case_count} cases")
