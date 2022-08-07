#!/bin/bash
Pexsi_DIR=/software/pexsi_v2.0.0
PARMETIS_DIR=/software/parmetis-4.0.3
SuperLU_DIR=/software/superlu_dist-7.2.0
INSTALL_DIR=/software/pexsi_v2.0.0


BDIR=build
rm -rf ${BDIR}

cmake -H${Pexsi_DIR} -B${BDIR} \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_BUILD_TYPE=Release \
    -DPEXSI_DEBUG_LEVEL=1 \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DSuperLU_DIST_PREFIX=${SuperLU_DIR} \
    -DParMETIS_PREFIX=${PARMETIS_DIR}

make -j 12 -C ${BDIR} install

