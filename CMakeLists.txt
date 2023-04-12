cmake_minimum_required(VERSION 3.18)

project(bitpunch C)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/bitpunch)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-Wall -Wextra -O2")

add_definitions(-DBPU_CONF_FULL -DERROR_L)
# uncomment to enable random generation of error vector in encryption
# add_definitions(-DGENERATE_ERROR_VECTOR)


set(SOURCES src/bitpunch/bitpunch.h
        src/bitpunch/code/codectx.c
        src/bitpunch/code/codectx.h
        src/bitpunch/code/goppa/goppa.c
        src/bitpunch/code/goppa/goppa.h
        src/bitpunch/code/goppa/goppatypes.c
        src/bitpunch/code/goppa/goppatypes.h
        src/bitpunch/code/qcmdpc/qcmdpc.c
        src/bitpunch/code/qcmdpc/qcmdpc.h
        src/bitpunch/code/qcmdpc/qcmdpctypes.c
        src/bitpunch/code/qcmdpc/qcmdpctypes.h
        src/bitpunch/config.h
        src/bitpunch/crypto/cca2/mecspointcheval.c
        src/bitpunch/crypto/cca2/mecspointcheval.h
        src/bitpunch/crypto/hash/polarssl/polarssl/sha512.h
        src/bitpunch/crypto/hash/polarssl/sha512.c
        src/bitpunch/crypto/hash/sha512.c
        src/bitpunch/crypto/hash/sha512.h
        src/bitpunch/crypto/mecs.c
        src/bitpunch/crypto/mecs.h
        src/bitpunch/crypto/mecsbasic/mecsbasic.c
        src/bitpunch/crypto/mecsbasic/mecsbasic.h
        src/bitpunch/crypto/mecsctx.c
        src/bitpunch/crypto/mecsctx.h
        src/bitpunch/crypto/padding/padding.c
        src/bitpunch/crypto/padding/padding.h
        src/bitpunch/debugio.c
        src/bitpunch/debugio.h
        src/bitpunch/errorcodes.h
        src/bitpunch/math/bigint.c
        src/bitpunch/math/bigint.h
        src/bitpunch/math/gf2.c
        src/bitpunch/math/gf2.h
        src/bitpunch/math/gf2x.c
        src/bitpunch/math/gf2x.h
        src/bitpunch/math/int.c
        src/bitpunch/math/int.h
        src/bitpunch/math/mathctx.c
        src/bitpunch/math/perm.c
        src/bitpunch/math/perm.h
        src/bitpunch/math/uni.c
        src/bitpunch/math/uni.h
        src/bitpunch/prng/prng.c
        src/bitpunch/prng/prng.h
        src/bitpunch/tools.c
        src/bitpunch/tools.h
        src/bitpunch/version.h
        src/main.c)

if(CMAKE_BUILD_TYPE MATCHES Release)
    message("release mode")
else()
    message("debug mode")
    add_definitions(-D_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
endif()

# add_link_options(-lm)

add_executable(bitpunch-main ${SOURCES})

target_link_libraries(bitpunch-main m)