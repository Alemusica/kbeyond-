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

# Fallback: cerca ext.h ovunque sotto MAX_SDK_ROOT
if (NOT C74SUPPORT)
  file(GLOB_RECURSE _ext_h_paths LIST_DIRECTORIES FALSE "${MAX_SDK_ROOT}/*.h")
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

# ====== UNICA DEFINIZIONE TARGET ======
add_library(kbeyond_tilde MODULE
  source/kbeyond_tilde.cpp
  source/dsp/prime_modes.cpp
  source/dsp/buffers.cpp
  source/dsp/filters.cpp
  source/dsp/early.cpp
  source/dsp/fdn.cpp
  source/dsp/decay.cpp
  source/dsp/detector.cpp
  source/dsp/mod.cpp
  source/dsp/mixing.cpp
)

target_compile_features(kbeyond_tilde PRIVATE cxx_std_17)

target_include_directories(kbeyond_tilde PRIVATE
  ${C74SUPPORT}/max-includes
  ${C74SUPPORT}/msp-includes
  ${CMAKE_CURRENT_SOURCE_DIR}/source
)

set_target_properties(kbeyond_tilde PROPERTIES
  PREFIX ""
  OUTPUT_NAME "kbeyond~"
)

if(APPLE)
  target_compile_definitions(kbeyond_tilde PRIVATE "C74_EXPORT=__attribute__((visibility(\"default\")))")
  set_target_properties(kbeyond_tilde PROPERTIES
    BUNDLE TRUE
    BUNDLE_EXTENSION "mxo"
    SUFFIX ".mxo"
  )

  # NIENTE framework Max: risolvi i simboli a runtime (robusto su macOS)
  target_link_options(kbeyond_tilde PRIVATE
    "-Wl,-undefined,dynamic_lookup"
    "-Wl,-dead_strip"
    "-Wl,-headerpad_max_install_names"
  )

  target_link_libraries(kbeyond_tilde PRIVATE
    "-framework Cocoa"
  )

  # Firma post-build (opzionale)
  set(CODESIGN_IDENTITY "" CACHE STRING "codesign identity for macOS")
  if (CODESIGN_IDENTITY)
    add_custom_command(TARGET kbeyond_tilde POST_BUILD
      COMMAND /usr/bin/codesign -s "${CODESIGN_IDENTITY}" --deep --force --timestamp=none
              "$<TARGET_BUNDLE_DIR:kbeyond_tilde>"
      COMMENT "Code signing kbeyond~.mxo with identity: ${CODESIGN_IDENTITY}"
      VERBATIM
    )
  endif()

elseif(WIN32)
  target_compile_definitions(kbeyond_tilde PRIVATE "C74_EXPORT=__declspec(dllexport)" WIN_VERSION=1)
  set_target_properties(kbeyond_tilde PROPERTIES SUFFIX ".mxe64")
  target_link_directories(kbeyond_tilde PRIVATE
    ${C74SUPPORT}/max-includes
    ${C74SUPPORT}/msp-includes
  )
  target_link_libraries(kbeyond_tilde PRIVATE
    MaxAPI
    MaxAudioAPI
  )
else()
  target_compile_definitions(kbeyond_tilde PRIVATE "C74_EXPORT=")
endif()

target_compile_options(kbeyond_tilde PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/EHsc /Zc:__cplusplus /permissive->
  $<$<CXX_COMPILER_ID:GNU,Clang>:-fvisibility=hidden>
)

option(KBEYOND_BUILD_TESTS "Build kbeyond~ unit tests" OFF)
if(KBEYOND_BUILD_TESTS)
    enable_testing()
    add_executable(kbeyond_tests
        tests/early_reflections_test.cpp
        tests/motion_detector_test.cpp
        source/kbeyond_tilde.cpp
        source/dsp/prime_modes.cpp
        source/dsp/buffers.cpp
        source/dsp/filters.cpp
        source/dsp/early.cpp
        source/dsp/fdn.cpp
        source/dsp/decay.cpp
        source/dsp/detector.cpp
        source/dsp/mod.cpp
        source/dsp/mixing.cpp)
    target_compile_definitions(kbeyond_tests PRIVATE KBEYOND_UNIT_TEST)
    target_include_directories(kbeyond_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source)
    target_compile_features(kbeyond_tests PRIVATE cxx_std_17)
    add_test(NAME kbeyond_early_mono_impulse COMMAND kbeyond_tests)
endif()

