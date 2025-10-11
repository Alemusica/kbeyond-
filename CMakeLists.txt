cmake_minimum_required(VERSION 3.18)
project(kbeyond_tilde LANGUAGES CXX)

# --- Resolve MAX_SDK_ROOT ---
if (NOT DEFINED MAX_SDK_ROOT)
  if (DEFINED ENV{MAX_SDK_ROOT})
    set(MAX_SDK_ROOT $ENV{MAX_SDK_ROOT})
  else()
    message(FATAL_ERROR "Setta MAX_SDK_ROOT o passa -DMAX_SDK_ROOT=/path/al/max-sdk")
  endif()
endif()

# --- Supporta entrambi i layout (max-sdk e max-sdk-base) ---
set(_c74_candidates
  "${MAX_SDK_ROOT}/max-sdk-base/c74support"
  "${MAX_SDK_ROOT}/source/max-sdk-base/c74support"
  "${MAX_SDK_ROOT}/source/c74support"
  "${MAX_SDK_ROOT}/c74support"
)

set(C74SUPPORT "")
foreach(candidate IN LISTS _c74_candidates)
  if (EXISTS "${candidate}/max-includes/ext.h")
    set(C74SUPPORT "${candidate}")
    break()
  endif()
endforeach()

# Alcune snapshot del Max SDK (in particolare i checkout git) possono usare
# layout leggermente diversi. Come fallback cerchiamo qualsiasi percorso che
# contenga `c74support/max-includes/ext.h` al suo interno.
if (NOT C74SUPPORT)
  file(GLOB_RECURSE _ext_h_paths LIST_DIRECTORIES FALSE
    "${MAX_SDK_ROOT}/*.h")
  foreach(_ext_h IN LISTS _ext_h_paths)
    if (_ext_h MATCHES "/c74support/max-includes/ext\\.h$")
      list(APPEND _c74support_hits "${_ext_h}")
    endif()
  endforeach()
  list(LENGTH _c74support_hits _ext_count)
  if (_ext_count GREATER 0)
    list(SORT _c74support_hits)
    list(GET _c74support_hits 0 _ext_h)
    get_filename_component(_max_includes "${_ext_h}" DIRECTORY)
    get_filename_component(C74SUPPORT "${_max_includes}" DIRECTORY)
  endif()
endif()

if (NOT C74SUPPORT)
  message(FATAL_ERROR "ext.h non trovato sotto ${MAX_SDK_ROOT}. Verifica che il Max SDK sia stato installato correttamente.")
endif()

if (NOT EXISTS "${C74SUPPORT}/msp-includes/z_dsp.h")
  message(FATAL_ERROR "z_dsp.h non trovato sotto ${C74SUPPORT}/msp-includes")
endif()

add_library(kbeyond~ MODULE
    source/kbeyond_tilde.cpp
)

target_compile_features(kbeyond~ PRIVATE cxx_std_17)

target_include_directories(kbeyond~ PRIVATE
    ${C74SUPPORT}/max-includes
    ${C74SUPPORT}/msp-includes
    ${CMAKE_CURRENT_SOURCE_DIR}/source
)

set_target_properties(kbeyond~ PROPERTIES
    PREFIX ""
    OUTPUT_NAME "kbeyond~"
)

if(APPLE)
    target_compile_definitions(kbeyond~ PRIVATE "C74_EXPORT=__attribute__((visibility(\"default\")))")
    set_target_properties(kbeyond~ PROPERTIES
        BUNDLE TRUE
        BUNDLE_EXTENSION "mxo"
        SUFFIX ".mxo"
    )
    target_link_directories(kbeyond~ PRIVATE
        ${C74SUPPORT}/msp-includes
    )
    target_link_options(kbeyond~ PRIVATE
        "-F${C74SUPPORT}/max-includes"
    )
    target_link_libraries(kbeyond~ PRIVATE
        "-framework Carbon"
        "-framework Cocoa"
        "-framework MaxAPI"
        "-framework MaxAudioAPI"
    )
elseif(WIN32)
    target_compile_definitions(kbeyond~ PRIVATE "C74_EXPORT=__declspec(dllexport)" WIN_VERSION=1)
    set_target_properties(kbeyond~ PROPERTIES SUFFIX ".mxe64")
    target_link_directories(kbeyond~ PRIVATE
        ${C74SUPPORT}/max-includes
        ${C74SUPPORT}/msp-includes
    )
    target_link_libraries(kbeyond~ PRIVATE
        MaxAPI
        MaxAudioAPI
    )
else()
    target_compile_definitions(kbeyond~ PRIVATE "C74_EXPORT=")
endif()

target_compile_options(kbeyond~ PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/EHsc /Zc:__cplusplus /permissive->
    $<$<CXX_COMPILER_ID:GNU,Clang>:-fvisibility=hidden>
)

option(KBEYOND_BUILD_TESTS "Build kbeyond~ unit tests" OFF)
if(KBEYOND_BUILD_TESTS)
    enable_testing()
    add_executable(kbeyond_tests
        tests/early_reflections_test.cpp
        source/kbeyond_tilde.cpp)
    target_compile_definitions(kbeyond_tests PRIVATE KBEYOND_UNIT_TEST)
    target_include_directories(kbeyond_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source)
    target_compile_features(kbeyond_tests PRIVATE cxx_std_17)
    add_test(NAME kbeyond_early_mono_impulse COMMAND kbeyond_tests)
endif()

