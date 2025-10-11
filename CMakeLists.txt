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
if (EXISTS "${MAX_SDK_ROOT}/source/c74support/max-includes/ext.h")
  set(C74SUPPORT "${MAX_SDK_ROOT}/source/c74support")
elseif (EXISTS "${MAX_SDK_ROOT}/c74support/max-includes/ext.h")
  set(C74SUPPORT "${MAX_SDK_ROOT}/c74support")
else()
  message(FATAL_ERROR "ext.h non trovato sotto ${MAX_SDK_ROOT}")
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

