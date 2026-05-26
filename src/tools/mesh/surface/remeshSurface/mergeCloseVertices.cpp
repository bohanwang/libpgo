#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <unordered_map>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;

// --- Types ---
using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;
using VD = Mesh::Vertex_index;
using FD = Mesh::Face_index;

// A simple 3D integer key for hash grid cells
struct I3
{
  int x, y, z;
  bool operator==(const I3 &o) const { return x == o.x && y == o.y && z == o.z; }
};
struct I3Hash
{
  std::size_t operator()(const I3 &k) const noexcept
  {
    // mix from boost::hash_combine style
    std::size_t h = 1469598103934665603ull;
    auto mix = [&](std::size_t v) {
      h ^= v + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    };
    mix(std::hash<int>{}(k.x));
    mix(std::hash<int>{}(k.y));
    mix(std::hash<int>{}(k.z));
    return h;
  }
};

// Merge (snap) vertices closer than eps; returns number of merged vertices.
std::size_t mergeCloseVertices(Mesh &mesh, double eps)
{
  if (num_vertices(mesh) == 0)
    return 0;

  const double cell = eps;                   // grid cell size
  const double eps2 = eps * eps;             // compare squared distances
  auto vpm = get(CGAL::vertex_point, mesh);  // property map

  // 1) Cluster vertices using a hash-grid and pick a representative per cluster
  std::unordered_map<I3, std::vector<VD>, I3Hash> grid;
  std::vector<VD> rep_of(mesh.number_of_vertices(), VD());  // map old vertex -> representative
  std::vector<bool> is_rep(mesh.number_of_vertices(), false);

  auto to_cell = [&](const Point &p) -> I3 {
    return I3{ (int)std::floor(p.x() / cell), (int)std::floor(p.y() / cell), (int)std::floor(p.z() / cell) };
  };

  std::size_t idx_counter = 0;  // for indexing VD to rep_of
  std::unordered_map<VD, std::size_t> v_to_idx;
  v_to_idx.reserve(mesh.number_of_vertices());
  for (VD v : vertices(mesh))
    v_to_idx[v] = idx_counter++;

  for (VD v : vertices(mesh)) {
    const Point &p = get(vpm, v);
    I3 c = to_cell(p);

    // Search current and neighboring 26 cells for an existing representative within eps
    VD chosen = VD();  // null
    for (int dx = -1; dx <= 1 && chosen == VD(); ++dx)
      for (int dy = -1; dy <= 1 && chosen == VD(); ++dy)
        for (int dz = -1; dz <= 1 && chosen == VD(); ++dz) {
          I3 n{ c.x + dx, c.y + dy, c.z + dz };
          auto it = grid.find(n);
          if (it == grid.end())
            continue;
          for (VD u : it->second) {
            const Point &q = get(vpm, u);
            const double d2 = CGAL::squared_distance(p, q);
            if (d2 <= eps2) {
              chosen = u;
              break;
            }
          }
        }

    if (chosen == VD()) {
      // Become a representative for this cell
      grid[c].push_back(v);
      is_rep[v_to_idx[v]] = true;
      rep_of[v_to_idx[v]] = v;
    }
    else {
      // Map to existing representative
      rep_of[v_to_idx[v]] = chosen;
    }
  }

  // Count how many will be merged (non-representatives)
  std::size_t merged = 0;
  for (VD v : vertices(mesh))
    if (rep_of[v_to_idx[v]] != v)
      ++merged;
  if (merged == 0)
    return 0;

  // 2) Build a new mesh with remapped vertices; drop degenerate faces
  Mesh out;
  // Map representative old VD -> new VD
  std::unordered_map<VD, VD> rep_old_to_new;
  rep_old_to_new.reserve(num_vertices(mesh));

  // Create new vertices for representatives
  for (VD v : vertices(mesh)) {
    if (rep_of[v_to_idx[v]] == v) {  // is representative
      VD nv = out.add_vertex(get(vpm, v));
      rep_old_to_new[v] = nv;
    }
  }

  auto rep_new = [&](VD old_v) -> VD {
    VD r = rep_of[v_to_idx[old_v]];
    return rep_old_to_new[r];
  };

  // Recreate faces with remapped vertices; skip degenerate faces
  std::size_t kept_faces = 0;
  for (FD f : faces(mesh)) {
    std::array<VD, 3> tri;
    int k = 0;
    for (auto h : CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
      tri[k++] = target(h, mesh);
    }
    if (k != 3)
      continue;  // ensure triangle

    VD a = rep_new(tri[0]);
    VD b = rep_new(tri[1]);
    VD c = rep_new(tri[2]);

    // Drop faces that became degenerate after merging (duplicate vertices or zero area)
    if (a == b || b == c || c == a)
      continue;

    // Optional: filter near-zero area triangles for extra robustness
    const Point &pa = out.point(a);
    const Point &pb = out.point(b);
    const Point &pc = out.point(c);
    Kernel::Vector_3 n = CGAL::cross_product(pb - pa, pc - pa);
    const double area2 = CGAL::sqrt(n.squared_length());  // 2*area magnitude
    if (area2 == 0.0)
      continue;

    if (out.add_face(a, b, c) == Mesh::null_face()) {
      // If orientation/duplication blocks insertion, try flipped order once
      if (out.add_face(a, c, b) == Mesh::null_face()) {
        // give up on this face
        continue;
      }
    }
    ++kept_faces;
  }

  // 3) Post-repair: remove isolated vertices, fix orientation, remove degeneracies
  PMP::remove_isolated_vertices(out);
  PMP::remove_degenerate_faces(out);
  PMP::remove_degenerate_edges(out);
  PMP::stitch_borders(out);
  PMP::orient(out);

  mesh = std::move(out);
  return merged;
}

// --- Minimal CLI: read mesh, merge, write mesh ---
// Supports OFF/OBJ/PLY/… (whatever your CGAL build supports via read_polygon_mesh/write_polygon_mesh)
int main(int argc, char **argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <input_mesh> <output_mesh>\n";
    std::cerr << "Example: " << argv[0] << " in.obj out.obj\n";
    return 1;
  }
  const std::string in = argv[1];
  const std::string out = argv[2];
  double eps = 0;

  Mesh mesh;
  if (!CGAL::IO::read_polygon_mesh(in, mesh) || CGAL::is_empty(mesh)) {
    std::cerr << "Failed to read mesh from '" << in << "'.\n";
    return 2;
  }

  if (argc == 4) {
    eps = std::stod(argv[3]);
  }
  else {
    // default: 0.5 of the min edge length
    double min_len = 1e100;
    for (FD f : faces(mesh)) {
      std::array<Point, 10> tri;
      int k = 0;
      for (auto h : CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
        tri[k] = mesh.point(target(h, mesh));
        k = (k + 1) % 10;
      }

      for (int i = 0; i < k; ++i) {
        int j = (i + 1) % k;
        double d = std::sqrt(CGAL::squared_distance(tri[i], tri[j]));
        if (d < min_len)
          min_len = d;
      }
    }
    eps = 0.5 * min_len;
    std::cout << "Using eps = " << eps << " (0.5 * min edge length)\n";
  }

  std::size_t merged = mergeCloseVertices(mesh, eps);
  if (!CGAL::IO::write_polygon_mesh(out, mesh,
        CGAL::parameters::stream_precision(17))) {
    std::cerr << "Failed to write mesh to '" << out << "'.\n";
    return 3;
  }

  std::cout << "Merged vertices within eps=" << eps << ": " << merged << "\n";
  std::cout << "Final mesh: " << num_vertices(mesh) << " V, "
            << num_faces(mesh) << " F\n";
  return 0;
}