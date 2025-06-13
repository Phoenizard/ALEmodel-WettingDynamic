
set(GMSH_ROOT /usr/local)
if ( DEFINED GMSH_ROOT)
    find_path( GMSH_INCLUDE_DIR Gmsh.h PATHS ${GMSH_ROOT}/include/gmsh
            NO_DEFAULT_PATH )
endif()

find_path( GMSH_INCLUDE_DIR Gmsh.h )
set(GMSH_INCLUDE_DIR /usr/local/include)
set(GMSH_LIB /usr/local/lib/libgmsh.so)
include_directories( ${GMSH_INCLUDE_DIR} )

find_library( GMSH_LIB NAMES Gmsh PATHS ${GMSH_ROOT}/lib )

add_amdis_executable(
  NAME "ALEModelAmdis2"
  SOURCES ale_model.cc)
target_link_libraries( ALEModelAmdis2 PUBLIC ${GMSH_LIB} )


#target_link_libraries(ALEModelAmdis2 PUBLIC Problems AMDiS /home/dune/gmsh/gmsh/build/libgmsh.a)