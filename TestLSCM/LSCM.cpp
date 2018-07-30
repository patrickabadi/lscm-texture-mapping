
#include "LSCM.h"
#include "pngio.h"
#include "ImageDraw.h"
#include "ScopeTimer.h"
#include "Helpers.h"

using namespace TEX_MAPPER;

IndexedMesh::IndexedMesh () 
  : in_facet ( false )
{
  facet_ptr.push_back ( 0 );
}

/*IndexedMesh::IndexedMesh ( const GL::Mesh& mesh )
  : IndexedMesh ()
{
  clear ();

  auto vcount = mesh.vertexCount ();
  auto tcount = mesh.triangleCount ();

  DebugOut ( "IndexedMesh::IndexedMesh vertices %d, triangles %d", vcount, tcount );

  for (int i = 0; i < vcount; i++)
  {
    auto v = mesh.vertex ( i );
    add_vertex ( vec3 ( v.x, v.y, v.z ), vec2 ( 0.0, 0.0 ) );
  }

  for (int i = 0; i < tcount; i++)
  {
    auto triangle = mesh.triangle ( i );

    begin_facet ();
    add_vertex_to_facet ( triangle.i1 );
    add_vertex_to_facet ( triangle.i2 );
    add_vertex_to_facet ( triangle.i3 );
    end_facet ();
  }
}*/

void IndexedMesh::load ( const std::string& file_name )
{
  /*std::ifstream input ( file_name.c_str () );
  clear ();
  while (input)
  {
    std::string line;
    std::getline ( input, line );
    std::stringstream line_input ( line );
    std::string keyword;
    line_input >> keyword;
    if (keyword == "v")
    {
      vec3 p;
      line_input >> p;
      add_vertex ( p, vec2 ( 0.0, 0.0 ) );
    }
    else if (keyword == "vt")
    {
      // Ignore tex vertices
    }
    else if (keyword == "f")
    {
      begin_facet ();
      while (line_input)
      {
        std::string s;
        line_input >> s;
        if (s.length () > 0)
        {
          std::stringstream v_input ( s.c_str () );
          NLuint index;
          v_input >> index;
          add_vertex_to_facet ( index - 1 );
          char c;
          v_input >> c;
          if (c == '/')
          {
            v_input >> index;
            // Ignore tex vertex index
          }
        }
      }
      end_facet ();
    }
  }
  std::cout << "Loaded " << vertex.size () << " vertices and "
    << nb_facets () << " facets" << std::endl;*/

  ScopeTimer t ( "IndexedMesh::Load" );

  std::ifstream in ( file_name, std::ios::in );
  if (!in)
  {
    DebugOut ( "Failed to load mesh with filename: " + file_name );
    return;
  }

  clear ();

  float x, y, z;
  float r, g, b;
  int i1, i2, i3, in1, in2, in3, iuv1, iuv2, iuv3;
  std::string line;
  int ret;

  while (std::getline ( in, line ))
  {
    if (line.empty () || line.size() < 3)
      continue;

    auto rest = line.substr ( 2 );

    int triangleMode = -1;

    if (line[0] == 'v')
    {
      // normal
      if (line[1] == 'n')
      {
        //sscanf_s ( rest.c_str (), "%f %f %f", &x, &y, &z );
        //normals.push_back ( Coordinate ( x, y, z ) );

      }
      else
      {
        ret = sscanf_s ( rest.c_str (), "%f %f %f %f %f %f", &x, &y, &z, &r, &g, &b );

        add_vertex ( vec3(x,y,z), vec2 ( 0.0, 0.0 ) );

        if (ret = 6)
        {
          colors.push_back ( color ( r, g, b ) );
        }

      }

    }
    else if (line[0] == 'f')
    {
      // triangleMode
      // -1: unknown so try stuff
      // 1: vertex and uv included
      // 2: vertex and normal include
      // 3: vertex and uv and normal included

      switch (triangleMode)
      {
      case -1:
        ret = sscanf_s ( rest.c_str (), "%d//%d %d//%d %d//%d", &i1, &in1, &i2, &in2, &i3, &in3 );
        if (ret == 6)
        {
          triangleMode = 2;
        }
        else
        {
          ret = sscanf_s ( rest.c_str (), "%d/%d %d/%d %d/%d", &i1, &iuv1, &i2, &iuv2, &i3, &iuv3 );
          if (ret == 6)
          {
            triangleMode = 1;
          }
          else
          {
            sscanf_s ( rest.c_str (), "%d/%d/%d %d/%d/%d %d/%d/%d", &i1, &iuv1, &in1, &i2, &iuv2, &in2, &i3, &iuv3, &in3 );
            triangleMode = 3;
          }
        }
        break;
      case 1:
        sscanf_s ( rest.c_str (), "%d/%d %d/%d %d/%d", &i1, &iuv1, &i2, &iuv2, &i3, &iuv3 );
        break;
      case 2:
        sscanf_s ( rest.c_str (), "%d//%d %d//%d %d//%d", &i1, &in1, &i2, &in2, &i3, &in3 );
        break;
      case 3:
        sscanf_s ( rest.c_str (), "%d/%d/%d %d/%d/%d %d/%d/%d", &i1, &iuv1, &in1, &i2, &iuv2, &in2, &i3, &iuv3, &in3 );
        break;
      default:
        break;
      }

      begin_facet ();
      add_vertex_to_facet ( i1 - 1 );
      add_vertex_to_facet ( i2 - 1 );
      add_vertex_to_facet ( i3 - 1 );
      end_facet ();

    }    
  }

  DebugOut ( "Loaded %d vertices and %d facets", vertex.size (), nb_facets () );
}

void IndexedMesh::save ( const std::string& file_name )
{
  std::ofstream out ( file_name.c_str () );
  for (NLuint v = 0; v < nb_vertices (); ++v)
  {
    out << "v " << vertex[v].point << std::endl;
  }
  for (NLuint v = 0; v < nb_vertices (); ++v)
  {
    out << "vt " << vertex[v].tex_coord << std::endl;
  }
  for (NLuint f = 0; f < nb_facets (); ++f)
  {
    NLuint nv = facet_nb_vertices ( f );
    out << "f ";
    for (NLuint lv = 0; lv < nv; ++lv)
    {
      NLuint v = facet_vertex ( f, lv );
      out << (v + 1) << "/" << (v + 1) << " ";
    }
    out << std::endl;
  }
  for (NLuint v = 0; v < nb_vertices (); ++v)
  {
    if (vertex[v].locked)
    {
      out << "# anchor " << v + 1 << std::endl;
    }
  }
}

void IndexedMesh::SaveTexture ( const std::string filename, int textureSize )
{
  auto png = common::pngio ( textureSize, textureSize, common::png_color_type::RGB_A );
  auto draw = ImageDraw ( png );

  ScopeTimer t ( "SaveTexture" );

  color defaultColor ( 255, 0, 0 );
  color* c1;
  color* c2;
  color* c3;

  for (NLuint f = 0; f < nb_facets (); ++f)
  {
    NLuint nv = facet_nb_vertices ( f );
    if (nv != 3)
    {
      DebugOut ( "triangle count changed" );
    }
      
    auto v1 = vertex[facet_vertex ( f, 0 )];
    auto v2 = vertex[facet_vertex ( f, 1 )];
    auto v3 = vertex[facet_vertex ( f, 2 )];

    if (colors.size () > 0)
    {
      c1 = &colors[facet_vertex ( f, 0 )];
      c2 = &colors[facet_vertex ( f, 1 )];
      c3 = &colors[facet_vertex ( f, 2 )];
    }
    else
    {
      c1 = c2 = c3 = &defaultColor;
    }

    draw.InterpolatedTriangle
    (
      vec2 ( v1.tex_coord.x * (double)textureSize, v1.tex_coord.y * textureSize ), *c1,
      vec2 ( v2.tex_coord.x * (double)textureSize, v2.tex_coord.y * textureSize ), *c2,
      vec2 ( v3.tex_coord.x * (double)textureSize, v3.tex_coord.y * textureSize ), *c3
    );

  }

  png.Save ( filename.c_str () );

}

/*void IndexedMesh::SaveUV ( GL::Mesh& mesh )
{
  mesh.ClearUV ();

  int vcount = (int)nb_vertices ();
  mesh.ResizeUV ( vcount );  

  DebugOut ( "IndexedMesh::SaveUV vertices %d", vcount);

  double minX, minY;
  double maxX, maxY;

  minX = std::numeric_limits<double>::max ();
  minY = std::numeric_limits<double>::max ();
  maxX = std::numeric_limits<double>::min ();
  maxY = std::numeric_limits<double>::min ();

  for (int v = 0; v < vcount; ++v)
  {
#ifdef _DEBUG
    if (v < 100)
      DebugOut ( "SaveUV %d: %f, %f", v, vertex[v].tex_coord.x, vertex[v].tex_coord.y );
#endif

    mesh.SetUV(v, Coordinate ( vertex[v].tex_coord.x, vertex[v].tex_coord.y, 0 ));

    minX = std::min ( minX, vertex[v].tex_coord.x );
    minY = std::min ( minY, vertex[v].tex_coord.y );
    maxX = std::max ( maxX, vertex[v].tex_coord.x );
    maxY = std::max ( maxY, vertex[v].tex_coord.y );

  }

  DebugOut("min (%f,%f) max (%f,%f)", minX, minY, maxX, maxY);
}*/

void LSCM::apply ()
{
  const int nb_eigens = 10;
  nlNewContext ();
  if (_spectral)
  {
    if (nlInitExtension ( "ARPACK" ))
    {
      std::cout << "ARPACK extension initialized"
        << std::endl;
    }
    else
    {
      std::cout << "Could not initialize ARPACK extension"
        << std::endl;
      exit ( -1 );
    }
    nlEigenSolverParameteri ( NL_EIGEN_SOLVER, NL_ARPACK_EXT );
    nlEigenSolverParameteri ( NL_NB_EIGENS, nb_eigens );
    nlEnable ( NL_VERBOSE );
  }
  NLuint nb_vertices = NLuint ( _mesh->vertex.size () );
  if (!_spectral)
  {
    project ();
  }
  nlSolverParameteri ( NL_NB_VARIABLES, NLint ( 2 * nb_vertices ) );
  nlSolverParameteri ( NL_LEAST_SQUARES, NL_TRUE );
  nlSolverParameteri ( NL_MAX_ITERATIONS, NLint ( 5 * nb_vertices ) );
  if (_spectral)
  {
    nlSolverParameterd ( NL_THRESHOLD, 0.0 );
  }
  else
  {
    nlSolverParameterd ( NL_THRESHOLD, 1e-6 );
  }
  nlBegin ( NL_SYSTEM );
  mesh_to_solver ();
  nlBegin ( NL_MATRIX );
  setup_lscm ();
  nlEnd ( NL_MATRIX );
  nlEnd ( NL_SYSTEM );
  std::cout << "Solving ..." << std::endl;

  if (_spectral)
  {
    nlEigenSolve ();
    for (NLuint i = 0; i < nb_eigens; ++i)
    {
      std::cerr << "[" << i << "] "
        << nlGetEigenValue ( i ) << std::endl;
    }

    // Find first "non-zero" eigenvalue
    double small_eigen = ::fabs ( nlGetEigenValue ( 0 ) );
    _eigen = 1;
    for (NLuint i = 1; i < nb_eigens; ++i)
    {
      if (::fabs ( nlGetEigenValue ( i ) ) / small_eigen > 1e3)
      {
        _eigen = i;
        break;
      }
    }
  }
  else
  {
    nlSolve ();
  }

  solver_to_mesh ();
  normalize_uv ();

  if (!_spectral)
  {
    double time;
    NLint iterations;
    nlGetDoublev ( NL_ELAPSED_TIME, &time );
    nlGetIntegerv ( NL_USED_ITERATIONS, &iterations );
    std::cout << "Solver time: " << time << std::endl;
    std::cout << "Used iterations: " << iterations << std::endl;
  }

  nlDeleteContext ( nlGetCurrent () );
}

void LSCM::setup_lscm ()
{
  for (NLuint f = 0; f < _mesh->nb_facets (); ++f)
  {
    setup_lscm ( f );
  }
}

/**
* \brief Creates the LSCM equations in OpenNL, related
*  with a given facet.
* \param[in] f the index of the facet.
* \details no-need to triangulate the facet,
*   we do that "virtually", by creating triangles
*   radiating around vertex 0 of the facet.
*   (however, this may be invalid for concave facets)
*/
void LSCM::setup_lscm ( NLuint f )
{
  NLuint nv = _mesh->facet_nb_vertices ( f );
  for (NLuint i = 1; i < nv - 1; ++i)
  {
    setup_conformal_map_relations (
      _mesh->facet_vertex ( f, 0 ),
      _mesh->facet_vertex ( f, i ),
      _mesh->facet_vertex ( f, i + 1 )
    );
  }
}

void LSCM::project_triangle (
  const vec3& p0,
  const vec3& p1,
  const vec3& p2,
  vec2& z0,
  vec2& z1,
  vec2& z2
)
{
  vec3 X = p1 - p0;
  X.normalize ();
  vec3 Z = cross ( X, (p2 - p0) );
  Z.normalize ();
  vec3 Y = cross ( Z, X );
  const vec3& O = p0;

  double x0 = 0;
  double y0 = 0;
  double x1 = (p1 - O).length ();
  double y1 = 0;
  double x2 = dot ( (p2 - O), X );
  double y2 = dot ( (p2 - O), Y );

  z0 = vec2 ( x0, y0 );
  z1 = vec2 ( x1, y1 );
  z2 = vec2 ( x2, y2 );
}

void LSCM::setup_conformal_map_relations (
  NLuint v0, NLuint v1, NLuint v2
)
{

  const vec3& p0 = _mesh->vertex[v0].point;
  const vec3& p1 = _mesh->vertex[v1].point;
  const vec3& p2 = _mesh->vertex[v2].point;

  vec2 z0, z1, z2;
  project_triangle ( p0, p1, p2, z0, z1, z2 );
  vec2 z01 = z1 - z0;
  vec2 z02 = z2 - z0;
  double a = z01.x;
  double b = z01.y;
  double c = z02.x;
  double d = z02.y;
  assert ( b == 0.0 );

  // Note  : 2*id + 0 --> u
  //         2*id + 1 --> v
  NLuint u0_id = 2 * v0;
  NLuint v0_id = 2 * v0 + 1;
  NLuint u1_id = 2 * v1;
  NLuint v1_id = 2 * v1 + 1;
  NLuint u2_id = 2 * v2;
  NLuint v2_id = 2 * v2 + 1;

  // Note : b = 0

  // Real part
  nlBegin ( NL_ROW );
  nlCoefficient ( u0_id, -a + c );
  nlCoefficient ( v0_id, b - d );
  nlCoefficient ( u1_id, -c );
  nlCoefficient ( v1_id, d );
  nlCoefficient ( u2_id, a );
  nlEnd ( NL_ROW );

  // Imaginary part
  nlBegin ( NL_ROW );
  nlCoefficient ( u0_id, -b + d );
  nlCoefficient ( v0_id, -a + c );
  nlCoefficient ( u1_id, -d );
  nlCoefficient ( v1_id, -c );
  nlCoefficient ( v2_id, a );
  nlEnd ( NL_ROW );
}

void LSCM::solver_to_mesh ()
{
  for (NLuint i = 0; i < _mesh->vertex.size (); ++i)
  {
    Vertex& it = _mesh->vertex[i];
    double u = _spectral ? nlMultiGetVariable ( 2 * i, _eigen )
      : nlGetVariable ( 2 * i );
    double v = _spectral ? nlMultiGetVariable ( 2 * i + 1, _eigen )
      : nlGetVariable ( 2 * i + 1 );
    it.tex_coord = vec2 ( u, v );
  }
}

void LSCM::normalize_uv ()
{
  double u_min = 1e30, v_min = 1e30, u_max = -1e30, v_max = -1e30;
  for (NLuint i = 0; i < _mesh->vertex.size (); ++i)
  {
    u_min = std::min ( u_min, _mesh->vertex[i].tex_coord.x );
    v_min = std::min ( v_min, _mesh->vertex[i].tex_coord.y );
    u_max = std::max ( u_max, _mesh->vertex[i].tex_coord.x );
    v_max = std::max ( v_max, _mesh->vertex[i].tex_coord.y );
  }
  double l = std::max ( u_max - u_min, v_max - v_min );
  for (NLuint i = 0; i < _mesh->vertex.size (); ++i)
  {
    _mesh->vertex[i].tex_coord.x -= u_min;
    _mesh->vertex[i].tex_coord.x /= l;
    _mesh->vertex[i].tex_coord.y -= v_min;
    _mesh->vertex[i].tex_coord.y /= l;
  }
}

void LSCM::mesh_to_solver ()
{
  for (NLuint i = 0; i < _mesh->vertex.size (); ++i)
  {
    Vertex& it = _mesh->vertex[i];
    double u = it.tex_coord.x;
    double v = it.tex_coord.y;
    nlSetVariable ( 2 * i, u );
    nlSetVariable ( 2 * i + 1, v );
    if (!_spectral && it.locked)
    {
      nlLockVariable ( 2 * i );
      nlLockVariable ( 2 * i + 1 );
    }
  }
}

void LSCM::project ()
{
  // Get bbox
  unsigned int i;

  double xmin = 1e30;
  double ymin = 1e30;
  double zmin = 1e30;
  double xmax = -1e30;
  double ymax = -1e30;
  double zmax = -1e30;

  for (i = 0; i < _mesh->vertex.size (); i++)
  {
    const Vertex& v = _mesh->vertex[i];
    xmin = std::min ( v.point.x, xmin );
    ymin = std::min ( v.point.y, ymin );
    zmin = std::min ( v.point.z, zmin );

    xmax = std::max ( v.point.x, xmax );
    ymax = std::max ( v.point.y, ymax );
    zmax = std::max ( v.point.z, zmax );
  }

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  vec3 V1, V2;

  // Find shortest bbox axis
  if (dx <= dy && dx <= dz)
  {
    if (dy > dz)
    {
      V1 = vec3 ( 0, 1, 0 );
      V2 = vec3 ( 0, 0, 1 );
    }
    else
    {
      V2 = vec3 ( 0, 1, 0 );
      V1 = vec3 ( 0, 0, 1 );
    }
  }
  else if (dy <= dx && dy <= dz)
  {
    if (dx > dz)
    {
      V1 = vec3 ( 1, 0, 0 );
      V2 = vec3 ( 0, 0, 1 );
    }
    else
    {
      V2 = vec3 ( 1, 0, 0 );
      V1 = vec3 ( 0, 0, 1 );
    }
  }
  else if (dz <= dx && dz <= dy)
  {
    if (dx > dy)
    {
      V1 = vec3 ( 1, 0, 0 );
      V2 = vec3 ( 0, 1, 0 );
    }
    else
    {
      V2 = vec3 ( 1, 0, 0 );
      V1 = vec3 ( 0, 1, 0 );
    }
  }

  // Project onto shortest bbox axis,
  // and lock extrema vertices

  Vertex* vxmin = NULL;
  double  umin = 1e30;
  Vertex* vxmax = NULL;
  double  umax = -1e30;

  for (i = 0; i < _mesh->vertex.size (); i++)
  {
    Vertex& V = _mesh->vertex[i];
    double u = dot ( V.point, V1 );
    double v = dot ( V.point, V2 );
    V.tex_coord = vec2 ( u, v );
    if (u < umin)
    {
      vxmin = &V;
      umin = u;
    }
    if (u > umax)
    {
      vxmax = &V;
      umax = u;
    }
  }

  vxmin->locked = true;
  vxmax->locked = true;
}