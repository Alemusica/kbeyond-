cmake_minimum_required(VERSION 3.18)
project(kbeyond_tilde LANGUAGES CXX)

if(NOT DEFINED ENV{MAX_SDK_ROOT})
    message(FATAL_ERROR "MAX_SDK_ROOT environment variable must point to the root of the Max SDK")
endif()

set(MAXSDK_ROOT $ENV{MAX_SDK_ROOT})

set(MAXSDK_MAX_INCLUDES "${MAXSDK_ROOT}/source/c74support/max-includes")
set(MAXSDK_MSP_INCLUDES "${MAXSDK_ROOT}/source/c74support/msp-includes")

if(NOT EXISTS "${MAXSDK_MAX_INCLUDES}/ext.h")
    message(FATAL_ERROR "Could not find ext.h in ${MAXSDK_MAX_INCLUDES}. Check that MAX_SDK_ROOT points to a valid Max SDK installation.")
endif()

if(NOT EXISTS "${MAXSDK_MSP_INCLUDES}/z_dsp.h")
    message(FATAL_ERROR "Could not find z_dsp.h in ${MAXSDK_MSP_INCLUDES}. Check that MAX_SDK_ROOT points to a valid Max SDK installation.")
endif()

add_library(kbeyond~ MODULE
    source/kbeyond_tilde.cpp
)

target_compile_features(kbeyond~ PRIVATE cxx_std_17)

target_include_directories(kbeyond~ PRIVATE
    ${MAXSDK_MAX_INCLUDES}
    ${MAXSDK_MSP_INCLUDES}
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
        ${MAXSDK_MSP_INCLUDES}
    )
    target_link_options(kbeyond~ PRIVATE
        "-F${MAXSDK_ROOT}/source/c74support/max-includes"
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
        ${MAXSDK_ROOT}/source/c74support/max-includes
        ${MAXSDK_ROOT}/source/c74support/msp-includes
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

