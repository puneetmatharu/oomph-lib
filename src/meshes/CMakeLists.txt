# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered meshes subdirectory")

# Define the headers
set(HEADERS
    one_d_mesh.h
    simple_rectangular_quadmesh.h
    rectangular_quadmesh.h
    quarter_circle_sector_mesh.h
    fish_mesh.h
    one_d_lagrangian_mesh.h
    two_layer_spine_mesh.h
    single_layer_spine_mesh.h
    bretherton_spine_mesh.h
    quarter_tube_mesh.h
    eighth_sphere_mesh.h
    channel_spine_mesh.h
    collapsible_channel_mesh.h
    triangle_mesh.h
    geompack_mesh.h
    tetgen_mesh.h
    simple_rectangular_tri_mesh.h
    simple_cubic_tet_mesh.h
    simple_cubic_mesh.h
    single_layer_cubic_spine_mesh.h
    channel_with_leaflet_mesh.h
    cylinder_with_flag_mesh.h
    fsi_driven_cavity_mesh.h
    hermite_element_quad_mesh.h
    circular_shell_mesh.h
    tube_mesh.h
    backward_step_mesh.h
    full_circle_mesh.h
    horizontal_single_layer_spine_mesh.h
    xda_tet_mesh.h
    thin_layer_brick_on_tet_mesh.h
    brick_from_tet_mesh.h
    quarter_pipe_mesh.h
    annular_mesh.h
    rectangle_with_hole_mesh.h
    quad_from_triangle_mesh.h
    refineable_tetgen_mesh.h
    gmsh_tet_mesh.h
    rectangle_with_moving_cylinder_mesh.h
    extruded_cube_mesh_from_quad_mesh_with_macro_elements.h)

# Define the template implementation headers (these are not source files as
# we're not trying to build them)
set(TEMPLATE_IMPL_HEADERS
    one_d_mesh.template.cc
    simple_rectangular_quadmesh.template.cc
    rectangular_quadmesh.template.cc
    quarter_circle_sector_mesh.template.cc
    fish_mesh.template.cc
    one_d_lagrangian_mesh.template.cc
    two_layer_spine_mesh.template.cc
    single_layer_spine_mesh.template.cc
    bretherton_spine_mesh.template.cc
    quarter_tube_mesh.template.cc
    eighth_sphere_mesh.template.cc
    channel_spine_mesh.template.cc
    collapsible_channel_mesh.template.cc
    triangle_mesh.template.cc
    geompack_mesh.template.cc
    tetgen_mesh.template.cc
    simple_rectangular_tri_mesh.template.cc
    simple_cubic_tet_mesh.template.cc
    simple_cubic_mesh.template.cc
    single_layer_cubic_spine_mesh.template.cc
    channel_with_leaflet_mesh.template.cc
    cylinder_with_flag_mesh.template.cc
    fsi_driven_cavity_mesh.template.cc
    hermite_element_quad_mesh.template.cc
    circular_shell_mesh.template.cc
    tube_mesh.template.cc
    backward_step_mesh.template.cc
    full_circle_mesh.template.cc
    horizontal_single_layer_spine_mesh.template.cc
    xda_tet_mesh.template.cc
    thin_layer_brick_on_tet_mesh.template.cc
    brick_from_tet_mesh.template.cc
    quarter_pipe_mesh.template.cc
    annular_mesh.template.cc
    rectangle_with_hole_mesh.template.cc
    quad_from_triangle_mesh.template.cc
    refineable_tetgen_mesh.template.cc
    gmsh_tet_mesh.template.cc
    rectangle_with_moving_cylinder_mesh.template.cc
    extruded_cube_mesh_from_quad_mesh_with_macro_elements.template.cc)

# Define the non-mesh headers
set(NON_MESH_HEADERS
    collapsible_channel_domain.h
    fish_domain.h
    channel_with_leaflet_domain.h
    quarter_tube_domain.h
    eighth_sphere_domain.h
    quarter_circle_sector_domain.h
    cylinder_with_flag_domain.h
    topologically_rectangular_domain.h
    tube_domain.h
    full_circle_domain.h
    quarter_pipe_domain.h
    annular_domain.h
    rectangle_with_hole_domain.h)

# Define the sources
set(SOURCES
    one_d_mesh.cc
    simple_rectangular_quadmesh.cc
    rectangular_quadmesh.cc
    quarter_circle_sector_mesh.cc
    fish_mesh.cc
    one_d_lagrangian_mesh.cc
    two_layer_spine_mesh.cc
    single_layer_spine_mesh.cc
    bretherton_spine_mesh.cc
    quarter_tube_mesh.cc
    eighth_sphere_mesh.cc
    channel_spine_mesh.cc
    collapsible_channel_mesh.cc
    triangle_mesh.cc
    geompack_mesh.cc
    tetgen_mesh.cc
    simple_rectangular_tri_mesh.cc
    simple_cubic_tet_mesh.cc
    simple_cubic_mesh.cc
    single_layer_cubic_spine_mesh.cc
    channel_with_leaflet_mesh.cc
    cylinder_with_flag_mesh.cc
    fsi_driven_cavity_mesh.cc
    hermite_element_quad_mesh.cc
    circular_shell_mesh.cc
    tube_mesh.cc
    backward_step_mesh.cc
    full_circle_mesh.cc
    horizontal_single_layer_spine_mesh.cc
    xda_tet_mesh.cc
    thin_layer_brick_on_tet_mesh.cc
    brick_from_tet_mesh.cc
    quarter_pipe_mesh.cc
    annular_mesh.cc
    rectangle_with_hole_mesh.cc
    quad_from_triangle_mesh.cc
    refineable_tetgen_mesh.cc
    gmsh_tet_mesh.cc
    rectangle_with_moving_cylinder_mesh.cc
    extruded_cube_mesh_from_quad_mesh_with_macro_elements.cc
    collapsible_channel_domain.cc
    fish_domain.cc
    channel_with_leaflet_domain.cc
    quarter_tube_domain.cc
    eighth_sphere_domain.cc
    quarter_circle_sector_domain.cc
    cylinder_with_flag_domain.cc
    topologically_rectangular_domain.cc
    tube_domain.cc
    full_circle_domain.cc
    quarter_pipe_domain.cc
    annular_domain.cc
    rectangle_with_hole_domain.cc)

# Set the name of the library
set(LIBNAME meshes)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

# The meshes library relies on the generic library.
target_link_libraries(meshes INTERFACE oomph::generic)

message(VERBOSE "Leaving meshes subdirectory")
# ------------------------------------------------------------------------------
