# Simulation Tick: Pre-Tick Setup

This walkthrough documents the data flow and function call chain in `src/api/runSimCore.cpp` before the solver starts advancing a simulation tick.

It currently stops after the initial plastic state vector is created. It does not yet cover:

- time integration
- collision and contact handling
- nonlinear solve iterations
- state update at the end of the tick

## Summary

The setup pipeline covered here does the following:

1. Load a volumetric mesh from disk.
2. Convert it into the solver-side `SimulationMesh`.
3. Load and scale the surface mesh.
4. Embed surface vertices into the volumetric mesh.
5. Build the sparse interpolation matrix `W`.
6. Create and initialize deformation models.
7. Create the FEM assembler.
8. Initialize the packed plastic state vector.

---

## 1. Load the Volumetric Mesh

### Role

This stage reads the tetrahedral or cubic volumetric mesh from `.veg` or `.vegb` and stores it in the `VolumetricMesh` family.

### Files and Functions

- File: `src/api/runSimCore.cpp`
- File: `src/core/scene/volumetricMesh/volumetricMesh.h`
- File: `src/core/scene/volumetricMesh/volumetricMesh.cpp`
- File: `src/core/scene/volumetricMesh/tetMesh.cpp`
- File: `src/core/scene/volumetricMesh/cubicMesh.cpp`
- Key constructors:
  - `TetMesh(...)`
  - `CubicMesh(...)`
  - `VolumetricMesh::loadFromAscii(...)`
  - `VolumetricMesh::loadFromBinaryGeneric(...)`

### Main Data Structures

`VolumetricMesh` stores the geometry and material partitioning:

```cpp
int                numVertices = 0;
std::vector<Vec3d> vertices;

int              numElementVertices = 0;
int              numElements        = 0;
std::vector<int> elements;

int                    numMaterials = 0;
int                    numSets      = 0;
int                    numRegions   = 0;
std::vector<Material*> materials;
std::vector<Set>       sets;
std::vector<Region>    regions;
std::vector<int>       elementMaterial;
```

Source: `src/core/scene/volumetricMesh/volumetricMesh.h`

### Geometry Layout

The geometry is stored as:

- `vertices`: global vertex table
- `elements`: flattened connectivity array

In index form, the connectivity of element $e$ is
$$
\mathcal{E}_e = (v_{e,0}, v_{e,1}, \dots, v_{e,k-1}),
$$
where $k = \texttt{numElementVertices}$ and each $v_{e,j}$ is a global vertex index.

The accessors make the layout explicit:

```cpp
inline const Vec3d& getVertex(int element, int vertex) const {
    return vertices[elements[element * numElementVertices + vertex]];
}

inline int getVertexIndex(int element, int vertex) const {
    return elements[element * numElementVertices + vertex];
}
```

Source: `src/core/scene/volumetricMesh/volumetricMesh.h`

### Material Mapping

Materials are assigned indirectly:

- `Set`: a named set of element ids
- `Region`: links one set to one material
- `elementMaterial[el]`: expanded per-element material index

The expansion happens here:

```cpp
void VolumetricMesh::propagateRegionsToElements() {
    for (int regionIndex = 0; regionIndex < numRegions; regionIndex++) {
        const Region& region        = regions[regionIndex];
        int           materialIndex = region.getMaterialIndex();

        const std::set<int>& setElements = sets[region.getSetIndex()].getElements();
        for (const auto& elt : setElements)
            elementMaterial[elt] = materialIndex;
    }
}
```

Source: `src/core/scene/volumetricMesh/volumetricMesh.cpp`

### `CubicMesh` Local Vertex Order

For cubic/hexahedral elements, local vertex order is fixed:

```text
000, 100, 110, 010, 001, 101, 111, 011
```

Equivalently, these are the eight corners of the unit cube:
$$
(0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1).
$$

It is encoded in both documentation and construction logic:

```cpp
int vtxI[8] = {0, 1, 1, 0, 0, 1, 1, 0};
int vtxJ[8] = {0, 0, 1, 1, 0, 0, 1, 1};
int vtxK[8] = {0, 0, 0, 0, 1, 1, 1, 1};
```

Source: `src/core/scene/volumetricMesh/cubicMesh.cpp`

### ASCII Loading Logic

When a `.veg` file is loaded, `VolumetricMesh::loadFromAscii(...)` uses a two-pass parser:

- first pass:
  - reads vertices and elements
  - identifies `TET` vs `CUBIC`
  - counts `*MATERIAL`, `*SET`, `*REGION`
- second pass:
  - parses material/set/region bodies
  - fills `elementMaterial`

The mesh type detection happens inside the `*ELEMENTS` section:

```cpp
if (strncmp(lineBuffer, "TET", 3) == 0) {
    if (elementType_)
        *elementType_ = TET;
} else if (strncmp(lineBuffer, "CUBIC", 5) == 0) {
    if (elementType_)
        *elementType_ = CUBIC;
} else {
    printf("Error: unknown mesh type %s in file %s\n", lineBuffer, filename);
    throw 3;
}
```

Source: `src/core/scene/volumetricMesh/volumetricMesh.cpp`

---

## 2. Convert to `SimulationMesh`

### Role

The scene-side `VolumetricMesh` representation is converted into the solver-side `SimulationMesh`.

### Files and Functions

- File: `src/api/runSimCore.cpp`
- File: `src/core/energy/solidDeformationModel/simulationMesh.h`
- File: `src/core/energy/solidDeformationModel/simulationMesh.cpp`
- Functions:
  - `SimulationMesh::createFromTetMesh(...)`
  - `SimulationMesh::createFromCubicMesh(...)`
  - `SimulationMesh::createTet(...)`
  - `SimulationMesh::createCubic(...)`
  - `SimulationMesh::createTyped<N>(...)`

### `SimulationMesh` Layout

Solver-side mesh storage:

```cpp
std::vector<EigenSupport::V3d>              vertices_;
std::vector<std::vector<int>>               elements_;
std::vector<ElementMaterialBinding>         elementMaterialBindings_;
std::vector<std::unique_ptr<SimulationMeshMaterial>> materials_;
SimulationMeshType                          meshType_;
```

Source: `src/core/energy/solidDeformationModel/simulationMesh.h`

### Cubic Conversion

For cubic meshes, conversion copies:

- all global vertices
- each hexahedral element's 8 global vertex ids
- per-element material assignments

At this stage only topology and constitutive metadata are transferred. State-dependent quantities such as displacement, deformation gradient, stress, and plastic state are still unset.

Key code:

```cpp
for (int ei = 0; ei < cubicMesh.getNumElements(); ei++) {
    elements.push_back(
        {cubicMesh.getVertexIndex(ei, 0), cubicMesh.getVertexIndex(ei, 1), cubicMesh.getVertexIndex(ei, 2),
         cubicMesh.getVertexIndex(ei, 3), cubicMesh.getVertexIndex(ei, 4), cubicMesh.getVertexIndex(ei, 5),
         cubicMesh.getVertexIndex(ei, 6), cubicMesh.getVertexIndex(ei, 7)});

    const VolumetricMeshes::VolumetricMesh::ENuMaterial* mat =
        downcastENuMaterial(cubicMesh.getElementMaterial(ei));
    materialOwners.push_back(std::make_unique<SimulationMeshENuMaterial>(mat->getE(), mat->getNu()));
    materials.push_back(materialOwners.back().get());
    elementMaterialBindings.push_back({ei, -1});
}

return createCubic(vertices, elements, elementMaterialBindings, materials);
```

Source: `src/core/energy/solidDeformationModel/simulationMesh.cpp`

### Conversion Call Chain

The `SimulationMesh` creation path is layered:

1. `SimulationMesh::createFromTetMesh(...)` or `SimulationMesh::createFromCubicMesh(...)`
   - scene-side mesh to solver-side raw arrays
2. `SimulationMesh::createTet(...)` or `SimulationMesh::createCubic(...)`
   - mesh-type-specific wrapper
3. `SimulationMesh::createTyped<N>(...)`
   - common implementation that allocates and fills the final `SimulationMesh`

So the two outer conversion functions are adapters, while `createTyped<N>(...)` is the actual constructor logic.

### `createFromTetMesh(...)` / `createFromCubicMesh(...)`

These two functions do the same kind of work for different element arities:

- copy all global vertices into a `std::vector<Vertex>`
- copy every element's connectivity into:
  - `TetElement = std::array<int, 4>` for tets
  - `CubicElement = std::array<int, 8>` for cubes
- convert source materials into `SimulationMeshMaterial` objects
- build `elementMaterialBindings`
- forward everything to `createTet(...)` or `createCubic(...)`

For tetrahedra, the connectivity packing is:

```cpp
elements.push_back(
    {tetMesh.getVertexIndex(ei, 0), tetMesh.getVertexIndex(ei, 1), tetMesh.getVertexIndex(ei, 2),
     tetMesh.getVertexIndex(ei, 3)});
```

For cubic meshes, the connectivity packing is:

```cpp
elements.push_back(
    {cubicMesh.getVertexIndex(ei, 0), cubicMesh.getVertexIndex(ei, 1), cubicMesh.getVertexIndex(ei, 2),
     cubicMesh.getVertexIndex(ei, 3), cubicMesh.getVertexIndex(ei, 4), cubicMesh.getVertexIndex(ei, 5),
     cubicMesh.getVertexIndex(ei, 6), cubicMesh.getVertexIndex(ei, 7)});
```

The material conversion path is the same in both:

```cpp
const VolumetricMeshes::VolumetricMesh::ENuMaterial* mat =
    downcastENuMaterial(sourceMesh.getElementMaterial(ei));
materialOwners.push_back(std::make_unique<SimulationMeshENuMaterial>(mat->getE(), mat->getNu()));
materials.push_back(materialOwners.back().get());
elementMaterialBindings.push_back({ei, -1});
```

The important thing here is that these functions are not yet constructing the final `SimulationMesh` object in place. They are building the four inputs expected by the common constructor:

- `vertices`
- `elements`
- `elementMaterialBindings`
- `materials`

### `createTet(...)` / `createCubic(...)`

These are very thin type wrappers:

```cpp
std::unique_ptr<SimulationMesh> SimulationMesh::createTet(
    std::span<const Vertex> vertices, std::span<const TetElement> elements,
    std::span<const ElementMaterialBinding> elementMaterialBindings,
    std::span<const SimulationMeshMaterial* const> materials) {
    return createTyped<4>(SimulationMeshType::TET, vertices, elements, elementMaterialBindings, materials);
}

std::unique_ptr<SimulationMesh> SimulationMesh::createCubic(
    std::span<const Vertex> vertices, std::span<const CubicElement> elements,
    std::span<const ElementMaterialBinding> elementMaterialBindings,
    std::span<const SimulationMeshMaterial* const> materials) {
    return createTyped<8>(SimulationMeshType::CUBIC, vertices, elements, elementMaterialBindings, materials);
}
```

Source: `src/core/energy/solidDeformationModel/simulationMesh.cpp`

Their only job is to:

- set the runtime mesh type (`TET` or `CUBIC`)
- lock in the compile-time local arity (`4` or `8`)

They do not contain any additional geometry or material logic.

### `createTyped<N>(...)`

`createTyped<N>(...)` is the real shared implementation for all supported mesh types.

Its job is:

1. validate the mesh-type / element-arity pair
2. validate that every element has a material binding
3. allocate a new `SimulationMesh`
4. copy vertices into `vertices_`
5. copy element connectivity into `elements_`
6. clone all material objects into owned storage
7. copy and validate `elementMaterialBindings_`
8. set the runtime `meshType_`

The function starts by validating consistency:

```cpp
if (expectedElementVertices(meshType) != static_cast<int>(N)) {
    throw std::invalid_argument("mesh type and typed element arity do not match");
}
if (elements.size() != elementMaterialBindings.size()) {
    throw std::invalid_argument("element/material binding size mismatch");
}
```

Then it allocates the mesh and copies the vertex table:

```cpp
auto mesh = std::unique_ptr<SimulationMesh>(new SimulationMesh());
mesh->vertices_.assign(vertices.begin(), vertices.end());
```

Then it converts the typed fixed-size element arrays into the internal storage format:

```cpp
mesh->elements_.assign(elements.size(), std::vector<int>(N, 0));
for (size_t ei = 0; ei < elements.size(); ei++) {
    memcpy(mesh->elements_[ei].data(), elements[ei].data(), sizeof(int) * N);
}
```

So externally the connectivity is passed as:

- `std::array<int, 4>` for tet
- `std::array<int, 8>` for cubic

but internally it is normalized to:

```cpp
std::vector<std::vector<int>> elements_;
```

Then it clones material objects into owned storage:

```cpp
mesh->materials_.clear();
mesh->materials_.reserve(materials.size());
for (size_t mi = 0; mi < materials.size(); mi++) {
    mesh->materials_.push_back(materials[mi]->clone());
}
```

This is why the temporary `materialOwners` vector used by `createFromTetMesh(...)` and `createFromCubicMesh(...)` is safe: the final `SimulationMesh` keeps its own clones.

Finally it copies the material bindings and validates that each `primary` or `secondary` value points to a valid material slot:

```cpp
mesh->elementMaterialBindings_.assign(elementMaterialBindings.begin(), elementMaterialBindings.end());
for (const ElementMaterialBinding& binding : mesh->elementMaterialBindings_) {
    validateMaterialIndex(binding.primary, static_cast<int>(mesh->materials_.size()), "primary");
    if (binding.secondary >= 0) {
        validateMaterialIndex(binding.secondary, static_cast<int>(mesh->materials_.size()), "secondary");
    }
}
```

and completes construction with:

```cpp
mesh->meshType_ = meshType;
return mesh;
```

### Data-Ownership Summary

The ownership transition across these functions is:

- source `TetMesh` / `CubicMesh`
  - owns scene-side geometry and materials
- temporary conversion vectors in `createFromTetMesh(...)` / `createFromCubicMesh(...)`
  - stage copied connectivity and temporary material wrappers
- final `SimulationMesh`
  - owns:
    - copied vertex array
    - copied connectivity
    - cloned `SimulationMeshMaterial` objects
    - copied element-material bindings

So after `createTyped<N>(...)` returns, the solver-side mesh is fully self-contained.

### Meaning of `elementMaterialBindings`

`ElementMaterialBinding` stores indices into the `materials` array:

```cpp
struct ElementMaterialBinding {
    int primary = -1;
    int secondary = -1;
};
```

Source: `src/core/energy/solidDeformationModel/simulationMesh.h`

The binding is validated and copied during `createTyped<N>(...)`:

```cpp
mesh->materials_.reserve(materials.size());
for (size_t mi = 0; mi < materials.size(); mi++) {
    mesh->materials_.push_back(materials[mi]->clone());
}

mesh->elementMaterialBindings_.assign(elementMaterialBindings.begin(), elementMaterialBindings.end());
for (const ElementMaterialBinding& binding : mesh->elementMaterialBindings_) {
    validateMaterialIndex(binding.primary, static_cast<int>(mesh->materials_.size()), "primary");
    if (binding.secondary >= 0) {
        validateMaterialIndex(binding.secondary, static_cast<int>(mesh->materials_.size()), "secondary");
    }
}
```

Source: `src/core/energy/solidDeformationModel/simulationMesh.cpp`

### Important Note

At the current code revision, `createFromTetMesh(...)` and `createFromCubicMesh(...)` still create one `SimulationMeshMaterial` object per element and set `primary = ei`.

That means:

- physically identical materials may be duplicated
- `primary` is a material array index, not a source mesh element id

---

## 3. Load and Scale the Surface Mesh

### Role

The surface mesh is used for rendering, output, and surface-side interpolation.

### Files and Functions

- File: `src/api/runSimCore.cpp`
- File: `src/core/scene/mesh/triMeshGeo.h`
- File: `src/core/scene/mesh/triMeshGeoLoader.cpp`
- Functions:
  - `TriMeshGeo::load(...)`
  - `TriMeshGeo::pos(...)`
  - `TriMeshGeo::numVertices()`

### Surface Mesh Data Structure

`TriMeshGeo` stores:

```cpp
std::vector<Vec3d> positions_;
std::vector<Vec3i> triangles_;
```

Source: `src/core/scene/mesh/triMeshGeo.h`

So:

- `positions_` is the global surface vertex table
- `triangles_` stores triples of vertex indices

The accessors are:

```cpp
int numVertices() const { return static_cast<int>(positions_.size()); }
int numTriangles() const { return static_cast<int>(triangles_.size()); }

const Vec3d& pos(int vtxID) const { return positions_[vtxID]; }
Vec3d&       pos(int vtxID) { return positions_[vtxID]; }
const Vec3i& tri(int triID) const { return triangles_[triID]; }
```

### OBJ Loading

The loader uses `tinyobjloader`:

```cpp
for (int vi = 0; vi < (int)attrib.vertices.size() / 3; vi++) {
    positions_.emplace_back(asVec3d(attrib.vertices.data() + vi * 3));
}

for (size_t s = 0; s < shapes.size(); s++) {
    size_t index_offset = 0;
    for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
        PGO_ALOG(fv == 3ull);

        Vec3i tri;
        for (size_t v = 0; v < fv; v++) {
            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
            tri[v]               = int(idx.vertex_index);
        }
        index_offset += fv;

        triangles_.emplace_back(tri);
    }
}
```

Source: `src/core/scene/mesh/triMeshGeoLoader.cpp`

### Scaling and Flattening

In `runSimCore.cpp`:

```cpp
Mesh::TriMeshGeo surfaceMesh;
if (surfaceMesh.load(mesh.surfaceMeshFilename) != true)
    return 1;

for (int vi = 0; vi < surfaceMesh.numVertices(); vi++) {
    surfaceMesh.pos(vi) *= mesh.scale;
}

int     surfn  = surfaceMesh.numVertices();
int     surfn3 = surfn * 3;
ES::VXd surfaceRestPositions(surfn3);
for (int vi = 0; vi < surfn; vi++) {
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
}
```

This transforms the vertex table from:

```text
[Vec3d, Vec3d, Vec3d, ...]
```

into:

```text
[x0 y0 z0 x1 y1 z1 ...]
```

If the surface has $m$ vertices with rest positions $s_i \in \mathbb{R}^3$, the flattened vector is
$$
\mathbf{s}_{\text{rest}} =
\begin{bmatrix}
s_0 \\
s_1 \\
\vdots \\
s_{m-1}
\end{bmatrix}
\in \mathbb{R}^{3m}.
$$

---

## 4. Embed Surface Vertices into the Volumetric Mesh

### Role

This stage computes, for each surface vertex, how it is represented as a linear combination of volumetric mesh vertices.

### Files and Functions

- File: `src/api/runSimCore.cpp`
- File: `src/core/scene/interpolationCoordinates/barycentricCoordinates.h`
- File: `src/core/scene/interpolationCoordinates/barycentricCoordinates.cpp`
- File: `src/core/scene/mesh/boundingBox.h`
- File: `src/core/scene/mesh/boundingVolumeTree.h`
- File: `src/core/scene/mesh/boundingVolumeTree.cpp`
- Functions:
  - `BarycentricCoordinates::BarycentricCoordinates(...)`
  - `BarycentricCoordinates::initializeInterpolationWeights(...)`
  - `BoundingBoxBVTree::buildByInertiaPartition(...)`
  - `BoundingBoxBVTree::getClosestBoundingBoxes(...)`
  - `VolumetricMesh::containsVertex(...)`
  - `VolumetricMesh::computeBarycentricWeights(...)`

### Constructor Entry Point

```cpp
InterpolationCoordinates::BarycentricCoordinates bc(
    surfaceMesh.numVertices(),
    surfaceRestPositions.data(),
    volumetricMesh.get());
```

This calls:

```cpp
BarycentricCoordinates::initializeInterpolationWeights(...)
```

### Embedding Storage

The embedding is stored in three arrays:

```cpp
int  numLocations       = 0;
int  numElementVertices = 0;
int  numCageVertices    = 0;

std::vector<int>    indices;
std::vector<double> weights;
std::vector<int>    elements;
```

Source: `src/core/scene/interpolationCoordinates/barycentricCoordinates.h`

The layout is:

- `elements[i]`
  - host element id for surface point `i`
- `indices[i * k + j]`
  - `j`-th volumetric vertex used by surface point `i`
- `weights[i * k + j]`
  - corresponding interpolation weight

where `k = numElementVertices`.

For surface point $p_i$, the embedding relation is
$$
p_i \approx \sum_{j=0}^{k-1} w_{i,j}\, x_{a_{i,j}},
$$
where $a_{i,j}$ is the `j`-th volumetric vertex index in `indices` and $x_{a_{i,j}} \in \mathbb{R}^3$ is its position.

### Build Per-Element Bounding Boxes

Key code:

```cpp
Mesh::BoundingBoxBVTree        bvTree;
std::vector<Mesh::BoundingBox> elementBBs(volumetricMesh->getNumElements());
for (int ei = 0; ei < volumetricMesh->getNumElements(); ei++) {
    BasicAlgorithms::ArrayRef<int> indexRef(volumetricMesh->getNumElementVertices(),
                                            volumetricMesh->getVertexIndices(ei));
    elementBBs[ei] = Mesh::BoundingBox(volumetricMesh->getVertices(), indexRef);
}
bvTree.buildByInertiaPartition(elementBBs);
```

Source: `src/core/scene/interpolationCoordinates/barycentricCoordinates.cpp`

This does:

1. For each volumetric element, gather its local vertex ids.
2. Use those global vertex positions to compute one AABB.
3. Build a binary BVH over all element AABBs.

### `BoundingBox` and `BoundingBoxBVTree`

`BoundingBox` stores:

```cpp
Vec3d bmin_, bmax_;
Vec3d center_, halfSides_;
```

Source: `src/core/scene/mesh/boundingBox.h`

`BoundingBoxBVTree` is built on top of `BoundingVolumeTreeBase`. Each node stores:

```cpp
struct Node {
    int              depth;
    std::vector<int> indices;
    BoundingBox      bb;
    std::vector<int> childrenIDs;
};
```

Source: `src/core/scene/mesh/boundingVolumeTreeBase.h`

The `buildByInertiaPartition(...)` routine:

- starts with one root node covering all element boxes
- computes per-box volume and inertia
- recursively partitions boxes into two groups along a dominant inertia direction
- builds a binary hierarchy for fast nearest-box queries

Core code:

```cpp
std::vector<int> rootIndices(numBBs);
std::iota(rootIndices.begin(), rootIndices.end(), 0);
BoundingBox rootBB(numBBs, boundingBoxes.data());
nodes.emplace_back(0, std::move(rootIndices));
nodes[0].bb = rootBB;

std::vector<double> bbVolumes(numBBs);
std::vector<Mat3d>  bbInertia(numBBs);
for (int i = 0; i < numBBs; i++) {
    bbVolumes[i] = boundingBoxes[i].volume();
    bbInertia[i] = getBoxInertiaTensorAroundCOM(boundingBoxes[i].sides(), bbVolumes[i]);
}
```

Source: `src/core/scene/mesh/boundingVolumeTree.cpp`

### Parallel Per-Vertex Embedding

After the BVH is built, every surface point is processed independently:

```cpp
tbb::parallel_for(0, numLocations, [&](int i) {
    ES::V3d                       pos = ES::Mp<const ES::V3d>(locations + 3 * i);
    thread_local std::vector<int> closestBBIDs;

    closestBBIDs.clear();
    bvTree.getClosestBoundingBoxes(elementBBs, pos, closestBBIDs);
    PGO_ALOG(closestBBIDs.size() > 0);

    bool posInsideElement = true;
    int  targetElementID  = -1;
    for (int eleID : closestBBIDs) {
        if (volumetricMesh->containsVertex(eleID, pos)) {
            targetElementID = eleID;
            break;
        }
    }

    if (targetElementID < 0) {
        posInsideElement = false;
        double closestDistance2 = DBL_MAX;
        for (int eleID : closestBBIDs) {
            Vec3d  center = volumetricMesh->getElementCenter(eleID);
            double dist2  = (pos - center).squaredNorm();
            if (dist2 < closestDistance2) {
                closestDistance2 = dist2;
                targetElementID  = eleID;
            }
        }
        numExternalVertices++;
    }

    elements[i] = targetElementID;

    memcpy(indices.data() + i * numElementVertices, volumetricMesh->getVertexIndices(targetElementID),
           sizeof(int) * numElementVertices);

    volumetricMesh->computeBarycentricWeights(targetElementID, pos, weights.data() + i * numElementVertices);
});
```

Source: `src/core/scene/interpolationCoordinates/barycentricCoordinates.cpp`

This means:

1. query the BVH for nearby element bounding boxes
2. find a candidate element that actually contains the point
3. if no element contains it, fall back to the nearest candidate element center
4. store the host element id
5. copy that element's global vertex ids into `indices`
6. compute interpolation weights into `weights`

### Weight Semantics

The actual formula depends on mesh type:

- `TetMesh`: tetrahedral barycentric coordinates, 4 weights
- `CubicMesh`: trilinear hexahedral weights, 8 weights

For cubic meshes, the weights are:

```cpp
weights[0] = (1 - alpha) * (1 - beta) * (1 - gamma);  // f000
weights[1] = (alpha) * (1 - beta) * (1 - gamma);      // f100
weights[2] = (alpha) * (beta) * (1 - gamma);          // f110
weights[3] = (1 - alpha) * (beta) * (1 - gamma);      // f010

weights[4] = (1 - alpha) * (1 - beta) * (gamma);  // f001
weights[5] = (alpha) * (1 - beta) * (gamma);      // f101
weights[6] = (alpha) * (beta) * (gamma);          // f111
weights[7] = (1 - alpha) * (beta) * (gamma);      // f011
```

Source: `src/core/scene/volumetricMesh/cubicMesh.cpp`

This is standard trilinear interpolation on the unit cube. If the local normalized coordinates are $(\alpha,\beta,\gamma)$, then
$$
p(\alpha,\beta,\gamma)
=
\sum_{c \in \{0,1\}^3} N_c(\alpha,\beta,\gamma)\, x_c,
$$
with shape functions
$$
N_{abc}(\alpha,\beta,\gamma)
=
\alpha^a (1-\alpha)^{1-a}
\beta^b (1-\beta)^{1-b}
\gamma^c (1-\gamma)^{1-c}.
$$

For tetrahedral meshes, the same idea uses barycentric coordinates:
$$
p = \lambda_0 x_0 + \lambda_1 x_1 + \lambda_2 x_2 + \lambda_3 x_3,
\qquad
\lambda_0 + \lambda_1 + \lambda_2 + \lambda_3 = 1.
$$

So one surface point always binds to one host element, but interpolates from multiple volumetric vertices of that element.

---

## 5. Build the Sparse Interpolation Matrix $W$

### Role

The embedding arrays are turned into a sparse matrix that maps volumetric vertex displacements to surface vertex displacements.

### Files and Functions

- File: `src/core/scene/interpolationCoordinates/barycentricCoordinates.cpp`
- File: `src/core/external/eigenSupport/EigenSupport.h`
- File: `src/core/external/eigenSupport/EigenSupport.cpp`
- Functions:
  - `BarycentricCoordinates::generateInterpolationMatrix()`
  - `EigenSupport::createWeightMatrix(...)`

### Matrix Generation Entry Point

```cpp
ES::SpMatD BarycentricCoordinates::generateInterpolationMatrix() const {
    PGO_ALOG(numCageVertices > 0);

    return ES::createWeightMatrix(numLocations * 3, numCageVertices * 3, numLocations, numElementVertices, nullptr,
                                  getEmbeddingVertexIndices().data(), getEmbeddingWeights().data(), 3);
}
```

Source: `src/core/scene/interpolationCoordinates/barycentricCoordinates.cpp`

The call arguments come directly from the embedding arrays built in the previous stage:

- `numLocations * 3`: one `x/y/z` row triplet per surface vertex
- `numCageVertices * 3`: one `x/y/z` column triplet per volumetric vertex
- `numLocations`: number of embedded surface vertices
- `numElementVertices`: number of interpolation weights per surface vertex
  - `4` for tetrahedral elements
  - `8` for cubic elements
- `nullptr`: use the default row vertex ids `0,1,2,...,numLocations-1`
- `getEmbeddingVertexIndices().data()`: the flattened `indices` array
- `getEmbeddingWeights().data()`: the flattened `weights` array
- `3`: inflate scalar weights to 3D vector form

### Matrix Dimensions and Meaning

`W` has shape:

```text
(3 * #surfaceVertices) x (3 * #volumetricVertices)
```

and satisfies:

```text
u_surface = W * u_volume
```

In symbols,
$$
\mathbf{u}_{\text{surf}} = W\,\mathbf{u}_{\text{vol}},
\qquad
W \in \mathbb{R}^{3m \times 3n},
$$
where $m = \#\text{surface vertices}$ and $n = \#\text{volumetric vertices}$.

At the scalar level, each surface vertex $i$ is embedded as
$$
p_i \approx \sum_{j=0}^{k-1} w_{i,j}\,x_{a_{i,j}},
$$
where:

- $k = \texttt{numElementVertices}$
- $a_{i,j}$ is the $j$-th volumetric vertex index stored in `indices[i * k + j]`
- $w_{i,j}$ is the matching interpolation weight stored in `weights[i * k + j]`

`W` is the matrix form of this rule, lifted from scalar coordinates to 3D vector coordinates.

### Sparse Assembly Logic

The actual sparse entries are created here:

```cpp
for (int vi = 0; vi < numVertices; vi++) {
    for (int vj = 0; vj < numWeightsPerVertex; vj++) {
        entries.emplace_back(vertexIndices[vi] * 3,
                             interpolationVertexIndices[vi * numWeightsPerVertex + vj] * 3,
                             interpolationWeights[vi * numWeightsPerVertex + vj]);
        entries.emplace_back(vertexIndices[vi] * 3 + 1,
                             interpolationVertexIndices[vi * numWeightsPerVertex + vj] * 3 + 1,
                             interpolationWeights[vi * numWeightsPerVertex + vj]);
        entries.emplace_back(vertexIndices[vi] * 3 + 2,
                             interpolationVertexIndices[vi * numWeightsPerVertex + vj] * 3 + 2,
                             interpolationWeights[vi * numWeightsPerVertex + vj]);
    }
}
```

Source: `src/core/external/eigenSupport/EigenSupport.cpp`

The helper first resolves the row-side vertex ids:

```cpp
std::vector<int> vertexIndices;
if (vtxIdx == nullptr) {
    vertexIndices.resize(numVertices);
    std::iota(vertexIndices.begin(), vertexIndices.end(), 0);
} else {
    vertexIndices.assign(vtxIdx, vtxIdx + numVertices);
}
```

Source: `src/core/external/eigenSupport/EigenSupport.cpp`

Because `generateInterpolationMatrix()` passes `nullptr`, the row vertex id for surface vertex `vi` is just `vi`.

For one scalar interpolation weight `w`, the code writes:

- one entry for `x -> x`
- one entry for `y -> y`
- one entry for `z -> z`

So for surface vertex `vi` and its `vj`-th interpolation neighbor,
$$
\begin{aligned}
W(3\,vi + 0,\; 3\,a_{vi,vj} + 0) &= w_{vi,vj}, \\
W(3\,vi + 1,\; 3\,a_{vi,vj} + 1) &= w_{vi,vj}, \\
W(3\,vi + 2,\; 3\,a_{vi,vj} + 2) &= w_{vi,vj}.
\end{aligned}
$$

No cross-coordinate entries are created. There is no `x -> y`, `x -> z`, or `y -> z` coupling in `W`.

### Relationship Between `indices` and Volumetric Vertices

The array `indices` does not store coordinates. It stores global volumetric vertex ids.

If the volumetric mesh vertex table is
$$
X = [x_0, x_1, \dots, x_{n-1}],
$$
where each $x_r \in \mathbb{R}^3$ is one volumetric vertex position, then
$$
\texttt{indices}[i \cdot k + j] = a_{i,j}
$$
means:

- surface vertex `i`
- uses the `j`-th local interpolation corner of its host element
- and that corner refers to global volumetric vertex `a_{i,j}`

So `indices` is a lookup table from:
$$
(\text{surface vertex id},\; \text{local corner id})
\;\longrightarrow\;
\text{global volumetric vertex id}.
$$

For example, if a host cubic element has local connectivity
$$
(10, 11, 16, 15, 35, 36, 41, 40),
$$
then those numbers are not coordinates. They are row numbers into the global volumetric vertex table:
$$
x_{10}, x_{11}, x_{16}, x_{15}, x_{35}, x_{36}, x_{41}, x_{40}.
$$

The interpolation rule always uses coordinates retrieved from the global vertex table:
$$
p_i \approx \sum_{j=0}^{k-1} w_{i,j}\,x_{\texttt{indices}[i \cdot k + j]}.
$$

### Worked Example

Suppose one surface vertex `i = 0` is embedded in a tetrahedral element with:

- volumetric vertex ids `[10, 13, 20, 25]`
- scalar weights `[0.1, 0.2, 0.3, 0.4]`

Then `indices` and `weights` store:

```text
indices  = [10, 13, 20, 25]
weights  = [0.1, 0.2, 0.3, 0.4]
```

and the first three rows of `W` get these nonzeros:

```text
W(0, 30) = 0.1   W(0, 39) = 0.2   W(0, 60) = 0.3   W(0, 75) = 0.4
W(1, 31) = 0.1   W(1, 40) = 0.2   W(1, 61) = 0.3   W(1, 76) = 0.4
W(2, 32) = 0.1   W(2, 41) = 0.2   W(2, 62) = 0.3   W(2, 77) = 0.4
```

This is because the global volumetric displacement vector is packed as
$$
\mathbf{u}_{\text{vol}}
=
[u_{0,x},u_{0,y},u_{0,z},u_{1,x},u_{1,y},u_{1,z},\dots]^T.
$$

So global volumetric vertex id `r` always maps to the three consecutive slots
$$
(3r,\; 3r+1,\; 3r+2),
$$
which store
$$
(u_{r,x},\; u_{r,y},\; u_{r,z}).
$$

In particular:
$$
\begin{aligned}
10 &\mapsto (30,31,32), \\
13 &\mapsto (39,40,41), \\
20 &\mapsto (60,61,62), \\
25 &\mapsto (75,76,77).
\end{aligned}
$$

So when `indices` contains `10`, that means this interpolation term reads from volumetric vertex `10`, whose packed `x/y/z` entries live at positions `30/31/32` in the global vector.

That is why:

- surface row `0` (`x`) couples to columns `30,39,60,75`
- surface row `1` (`y`) couples to columns `31,40,61,76`
- surface row `2` (`z`) couples to columns `32,41,62,77`

The same mapping can be summarized as:

| Volumetric vertex id | Packed columns in $\mathbf{u}_{\text{vol}}$ | Physical meaning |
| --- | --- | --- |
| `10` | `30, 31, 32` | $u_{10,x}, u_{10,y}, u_{10,z}$ |
| `13` | `39, 40, 41` | $u_{13,x}, u_{13,y}, u_{13,z}$ |
| `20` | `60, 61, 62` | $u_{20,x}, u_{20,y}, u_{20,z}$ |
| `25` | `75, 76, 77` | $u_{25,x}, u_{25,y}, u_{25,z}$ |

Equivalently, the corresponding volumetric vertex positions are:
$$
x_{10},\; x_{13},\; x_{20},\; x_{25} \in \mathbb{R}^3,
$$
and the scalar interpolation rule for the surface rest position is:
$$
p_0 \approx 0.1\,x_{10} + 0.2\,x_{13} + 0.3\,x_{20} + 0.4\,x_{25}.
$$

Multiplying this matrix by the packed volumetric displacement vector
$$
\mathbf{u}_{\text{vol}}
=
[u_{0,x},u_{0,y},u_{0,z},u_{1,x},u_{1,y},u_{1,z},\dots]^T
$$
produces:
$$
\begin{aligned}
u_{\text{surf},0,x} &= 0.1\,u_{10,x} + 0.2\,u_{13,x} + 0.3\,u_{20,x} + 0.4\,u_{25,x}, \\
u_{\text{surf},0,y} &= 0.1\,u_{10,y} + 0.2\,u_{13,y} + 0.3\,u_{20,y} + 0.4\,u_{25,y}, \\
u_{\text{surf},0,z} &= 0.1\,u_{10,z} + 0.2\,u_{13,z} + 0.3\,u_{20,z} + 0.4\,u_{25,z}.
\end{aligned}
$$

For cubic elements the pattern is identical, except each surface vertex contributes `8` weights instead of `4`, so each `x/y/z` row receives `8` nonzeros.

### Final Sparse Matrix Construction

After all triplets are collected, Eigen builds the sparse matrix in one shot:

```cpp
A.resize(numRows, numCols);
A.setFromTriplets(entries.begin(), entries.end());
```

and the wrapper returns it as `W`:

```cpp
SpMatD W;
pgo::EigenSupport::createWeightMatrix(numRows, numCols, numVertices, numWeightsPerVertex, vtxIdx,
                                      interpolationVertexIndices, interpolationWeights, inflate, W);
return W;
```

Source: `src/core/external/eigenSupport/EigenSupport.cpp`
---

## 6. Initialize the Deformation Model Manager

### Role

This stage creates the per-element constitutive models used by the solid deformation system.

### Files and Functions

- File: `src/api/runSimCore.cpp`
- File: `src/core/energy/solidDeformationModel/deformationModelManager.h`
- File: `src/core/energy/solidDeformationModel/deformationModelManager.cpp`
- Functions:
  - `DeformationModelManager::setMesh(...)`
  - `DeformationModelManager::init(...)`

### Construction and Mesh Registration

```cpp
std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm =
    std::make_shared<SolidDeformationModel::DeformationModelManager>();

dmm->setMesh(simMesh.get(), nullptr, nullptr);
```

Source: `src/api/runSimCore.cpp`

`setMesh(...)` registers the `SimulationMesh` and optional fiber direction arrays:

```cpp
void setMesh(const SimulationMesh* simulationMesh,
             const double* elementFiberDirections = nullptr,
             const double* vertexFiberDirections = nullptr);
```

Source: `src/core/energy/solidDeformationModel/deformationModelManager.h`

Here both fiber direction inputs are `nullptr`, so no anisotropic fiber frames are provided externally.

### Model Initialization

```cpp
dmm->init(
    pgo::SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
    solver.elasticMaterial,
    1);
```

This configures:

- plastic model type: `VOLUMETRIC_DOF6`
- elastic material type: `solver.elasticMaterial`
- SPD-related stabilization: enabled by `1`

Internally, `DeformationModelManagerImpl` stores:

```cpp
const SimulationMesh* simulationMesh;
std::vector<DeformationModel*> elementFEMs;
std::vector<ElasticModel*>     elementMaterials;
```

plus the concrete elastic and plastic model pools.

Source: `src/core/energy/solidDeformationModel/deformationModelManager.cpp`

After `init(...)`, each element has a concrete `DeformationModel` that combines:

- element geometry
- elastic constitutive behavior
- plastic constitutive behavior

### Object Relationship Chain

The runtime object stack built so far can be summarized as:

```text
SimulationMesh
  -> ElasticModel
  -> PlasticModel
  -> DeformationModel
  -> Assembler
```

The meaning of each layer is:

- `SimulationMesh`
  - solver-side mesh geometry and connectivity
  - provides rest vertex positions, element type, local vertex indices, and per-element material bindings
- `ElasticModel`
  - constitutive law for the recoverable elastic response
  - examples include linear, Stable Neo-Hookean, invariant-based StVK, volume penalty, and their combinations
- `PlasticModel`
  - constitutive law for the per-element internal plastic state
  - examples include `PlasticModel3DConstant`, `PlasticModel3D3DOF`, and `PlasticModel3D6DOF`
- `DeformationModel`
  - the full element-level FEM object
  - combines rest geometry, one elastic model, and one plastic model
  - computes element energy, forces, and tangent matrices
- `Assembler`
  - loops over all element-level `DeformationModel`s
  - assembles local contributions into global energy, gradient, and Hessian

In particular, `elementMaterials[ele]` and `elementFEMs[ele]` are not independent. The deformation model is constructed *from* the elastic material object and the plastic model object for that same element.

For example, for cubic elements:

```cpp
data->elementFEMs[ele] =
    new CubicMeshDeformationModel(restPosition.data(), data->elementMaterials[ele], pm);
```

and for tetrahedral elements:

```cpp
data->elementFEMs[ele] =
    new TetMeshDeformationModel(restPosition.data(), restPosition.data() + 3, restPosition.data() + 6,
                                restPosition.data() + 9, data->elementMaterials[ele], pm);
```

Source: `src/core/energy/solidDeformationModel/deformationModelManager.cpp`

So the dependency direction is:

```text
SimulationMesh  -> provides local geometry
ElasticModel    -> provides elastic constitutive response
PlasticModel    -> provides internal-variable parameterization
DeformationModel -> wraps geometry + elastic model + plastic model
Assembler       -> consumes all per-element deformation models
```

At the element level, the assembled local energy can be viewed abstractly as
$$
E_e = E_e(\mathbf{x}_e, \mathbf{X}_e, \mathbf{a}_e; \mathcal{M}^{\text{elastic}}_e, \mathcal{M}^{\text{plastic}}_e),
$$
where:

- $\mathbf{x}_e$ is the current element state
- $\mathbf{X}_e$ is the rest element geometry
- $\mathbf{a}_e$ is the plastic internal state
- $\mathcal{M}^{\text{elastic}}_e$ is the element's elastic model object
- $\mathcal{M}^{\text{plastic}}_e$ is the element's plastic model object

The assembler then lifts these element energies to the global level by summing over all elements.

---

## 7. Initialize the FEM Assembler

### Role

This stage creates the object that assembles element-level energy terms into global vectors and sparse matrices.

### Files and Functions

- File: `src/api/runSimCore.cpp`
- File: `src/core/energy/solidDeformationModel/deformationModelAssembler.h`
- File: `src/core/energy/solidDeformationModel/deformationModelAssembler.cpp`
- Functions:
  - `DeformationModelAssembler::DeformationModelAssembler(...)`
  - `computeEnergy(...)`
  - `computeGradient(...)`
  - `computeHessian(...)`

### Construction

```cpp
std::vector<double> elementWeights(simMesh->getNumElements(), 1.0);
std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> assembler =
    std::make_shared<SolidDeformationModel::DeformationModelAssembler>(dmm, elementWeights.data());
```

Source: `src/api/runSimCore.cpp`

The constructor signature is:

```cpp
DeformationModelAssembler(std::shared_ptr<const DeformationModelManager> dm,
                          const double* elementFlags = nullptr);
```

Source: `src/core/energy/solidDeformationModel/deformationModelAssembler.h`

`elementWeights` is copied into the assembler as `elementFlags`. At this call site, every element is given weight `1.0`, so all elements contribute normally.

At the mathematical level, the assembled global energy has the form
$$
E(\mathbf{x}, \mathbf{a}, \mathbf{b})
=
\sum_{e=0}^{N_e-1} \omega_e\, E_e(\mathbf{x}_e, \mathbf{a}_e, \mathbf{b}_e),
$$
where $\omega_e = \texttt{elementWeights[e]}$, $\mathbf{x}_e$ is the local element state, $\mathbf{a}_e$ is the plastic parameter block, and $\mathbf{b}_e$ is the elastic parameter block.

### Constructor Responsibilities

The constructor extracts global dimensions:

```cpp
nele    = deformationModelManager->getMesh()->getNumElements();
nvtx    = deformationModelManager->getMesh()->getNumVertices();
neleVtx = deformationModelManager->getMesh()->getNumElementVertices();
n3      = 3 * nvtx;
```

Source: `src/core/energy/solidDeformationModel/deformationModelAssembler.cpp`

It also copies rest positions:

```cpp
restPositions.resize(n3);
for (int vi = 0; vi < deformationModelManager->getMesh()->getNumVertices(); vi++) {
    ES::V3d p;
    deformationModelManager->getMesh()->getVertex(vi, p.data());
    restPositions.segment<3>(vi * 3) = p;
}
```

and caches per-element FEM model pointers:

```cpp
for (int i = 0; i < nele; i++) {
    femModels.push_back(deformationModelManager->getDeformationModel(i));
    data->elementCacheData.push_back(femModels.back()->allocateCacheData());
}
```

### Sparse Template Precomputation

The assembler also precomputes a sparsity pattern for the global stiffness matrix:

```cpp
for (int ele = 0; ele < nele; ele++) {
    const int* vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);

    for (int vi = 0; vi < neleVtx; vi++) {
        for (int vj = 0; vj < neleVtx; vj++) {
            for (int dofi = 0; dofi < 3; dofi++) {
                for (int dofj = 0; dofj < 3; dofj++) {
                    if (vertexIndices[vi] >= 0 && vertexIndices[vj] >= 0)
                        entries.emplace_back(vertexIndices[vi] * 3 + dofi, vertexIndices[vj] * 3 + dofj, 1.0);
                }
            }
        }
    }
}

KTemplate.resize(n3, n3);
KTemplate.setFromTriplets(entries.begin(), entries.end());
```

Source: `src/core/energy/solidDeformationModel/deformationModelAssembler.cpp`

So before any actual state-dependent solve happens, the assembler already knows the global sparsity structure it will fill later.

This template corresponds to the sparsity structure of the global tangent matrix
$$
K(\mathbf{x}) = \frac{\partial^2 E}{\partial \mathbf{x}^2}.
$$

---

## 8. Initialize the Plastic State

### Role

This stage creates the packed per-element plastic state vector and initializes every element to "no plastic deformation".

### Files and Functions

- File: `src/api/runSimCore.cpp`
- File: `src/core/energy/solidDeformationModel/deformationModelManager.h`
- File: `src/core/energy/solidDeformationModel/deformationModelEnergy.h`
- File: `src/core/energy/solidDeformationModel/deformationModelEnergy.cpp`
- File: `src/core/energy/solidDeformationModel/plasticModel3DDeformationGradient.md`
- Key runtime calls:
  - `DeformationModelManager::getDeformationModel(...)`
  - `DeformationModel::getPlasticModel()`
  - `PlasticModel3DDeformationGradient::toParam(...)`
  - `DeformationModelEnergy::DeformationModelEnergy(...)`
  - `DeformationModelEnergy::setPlasticParams(...)`

### Initialization Code

```cpp
int n    = simMesh->getNumVertices();
int n3   = n * 3;
int nele = simMesh->getNumElements();

ES::VXd plasticity(nele * 6);
ES::M3d I = ES::M3d::Identity();
for (int ei = 0; ei < nele; ei++) {
    const SolidDeformationModel::PlasticModel3DDeformationGradient* pm =
        dynamic_cast<const SolidDeformationModel::PlasticModel3DDeformationGradient*>(
            dmm->getDeformationModel(ei)->getPlasticModel());
    if (!pm) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic model is not of type PlasticModel3DDeformationGradient.");
        return 1;
    }
    pm->toParam(I.data(), plasticity.data() + ei * dmm->getNumPlasticParameters());
}
```

Source: `src/api/runSimCore.cpp`

### Why `6` Parameters

The plasticity vector uses:

```cpp
ES::VXd plasticity(nele * 6);
```

because the chosen plastic model type is:

```cpp
DeformationModelPlasticMaterial::VOLUMETRIC_DOF6
```

The `6` is the number of plastic internal parameters per element. It is not related to the number of element vertices:

- cubic element geometry: 8 vertices
- plastic model state: 6 parameters

These are different notions:

- geometric DOFs come from mesh topology
- plastic parameters come from constitutive model parameterization

### Identity Plastic State

`I` is the identity deformation gradient:

```cpp
ES::M3d I = ES::M3d::Identity();
```

This represents the initial "no plastic deformation" state.

In multiplicative plasticity notation, the total deformation gradient is split as
$$
F = F_e F_p,
$$
and the initialization here sets
$$
F_p = I,
$$
which means the initial state has no accumulated plastic deformation.

This is an explicit initialization choice baked into the setup code: the simulation starts by assuming every element has zero plastic history. In that sense, the code is hard-coding the initial plastic state to
$$
F_p = I.
$$

What is *not* hard-coded is the raw parameter vector itself. The code does not directly write six numbers such as `[0,0,0,0,0,0]` into `plasticity`. Instead, it hard-codes the physical initial condition "no plastic deformation" and then asks the plastic model to encode that state into its own parameterization:

```cpp
pm->toParam(I.data(), plasticity.data() + ei * dmm->getNumPlasticParameters());
```

So the logic is:

1. choose the initial physical plastic state `F_p = I`
2. convert that matrix state into the model's packed parameter form
3. store the encoded result in the global `plasticity` vector

If the code later needs to support nontrivial initial plastic history, this is the place that would need to change. For example, instead of always using `I`, the setup could load an elementwise initial `F_p` field from disk, from a previous frame, or from a preprocessing stage.

For each element:

1. fetch the deformation model from `dmm`
2. fetch its plastic model
3. verify it is a `PlasticModel3DDeformationGradient`
4. convert `F_p = I` into the model's packed parameter form with `toParam(...)`

So the final `plasticity` vector is:

```text
[element0 plastic params | element1 plastic params | ...]
```

with each block initialized from the identity plastic deformation gradient.

If each element has $q = \texttt{dmm->getNumPlasticParameters()}$ plastic parameters, then the packed global plastic state is
$$
\mathbf{a} =
\begin{bmatrix}
\mathbf{a}_0 \\
\mathbf{a}_1 \\
\vdots \\
\mathbf{a}_{N_e-1}
\end{bmatrix}
\in \mathbb{R}^{q N_e}.
$$

### Build the Packed Rest Position Vector

After the plastic state is initialized, the code builds the global rest-position vector:

```cpp
ES::VXd restPosition(n3);
for (int vi = 0; vi < n; vi++) {
    double p[3];
    simMesh->getVertex(vi, p);
    restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
}
```

Source: `src/api/runSimCore.cpp`

This converts the mesh vertex table into one packed vector:
$$
\mathbf{x}_{\text{rest}}
=
[x_{0},y_{0},z_{0},x_{1},y_{1},z_{1},\dots,x_{n-1},y_{n-1},z_{n-1}]^T
\in \mathbb{R}^{3n}.
$$

The source of truth here is still `simMesh`. The code simply flattens the solver-side mesh coordinates into the vector format expected by the energy and optimizer code.

### Create `DeformationModelEnergy`

The next step wraps the assembler, the rest configuration, and the plastic state into one global energy object:

```cpp
std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy =
    std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, &restPosition, 0);
elasticEnergy->setPlasticParams(plasticity);
```

Source: `src/api/runSimCore.cpp`

At this point:

- `assembler` knows how to assemble per-element energies, gradients, and Hessians
- `restPosition` supplies the reference configuration
- `plasticity` supplies the current per-element internal variables

So `elasticEnergy` becomes the object that can evaluate the global solid energy for a current state vector.

At a high level, it represents a function of the form
$$
E(\mathbf{x}, \mathbf{a})
=
\sum_{e=0}^{N_e-1} \omega_e\,E_e(\mathbf{x}_e,\mathbf{X}_e,\mathbf{a}_e),
$$
where:

- $\mathbf{x}_e$ is the current element state extracted from the global unknown
- $\mathbf{X}_e$ is the element's rest configuration extracted from `restPosition`
- $\mathbf{a}_e$ is the element's packed plastic parameter block
- $\omega_e$ is the element assembly weight

### What `setPlasticParams(...)` Does

The call

```cpp
elasticEnergy->setPlasticParams(plasticity);
```

copies or registers the packed plastic state with the energy object so that later evaluations use the correct per-element internal variables.

So this final part of the setup does three things in sequence:

1. initialize the per-element plastic internal state to `F_p = I`
2. flatten the mesh vertices into the packed rest-position vector
3. create the global deformation energy object and attach the plastic state

After this point, the solver has all the static ingredients needed to start evaluating:

- total elastic/plastic energy
- internal force / energy gradient
- tangent stiffness / Hessian

---

## Key Call Chain

- `TetMesh(...)` or `CubicMesh(...)`
- `SimulationMesh::createFromTetMesh(...)` or `SimulationMesh::createFromCubicMesh(...)`
- `TriMeshGeo::load(...)`
- `BarycentricCoordinates::initializeInterpolationWeights(...)`
- `BoundingBoxBVTree::buildByInertiaPartition(...)`
- `BoundingBoxBVTree::getClosestBoundingBoxes(...)`
- `VolumetricMesh::computeBarycentricWeights(...)`
- `BarycentricCoordinates::generateInterpolationMatrix()`
- `DeformationModelManager::setMesh(...)`
- `DeformationModelManager::init(...)`
- `DeformationModelAssembler::DeformationModelAssembler(...)`
- `PlasticModel3DDeformationGradient::toParam(...)`

## Validation Checklist

- Confirm `W.rows() == 3 * surfaceMesh.numVertices()`.
- Confirm `W.cols() == 3 * volumetricMesh->getNumVertices()`.
- For each embedded surface point, confirm `indices` has exactly `numElementVertices` entries.
- For cubic meshes, confirm a surface point binds to one host cube but interpolates from 8 cube vertices.
- Confirm `plasticity.size() == nele * dmm->getNumPlasticParameters()`.
- Confirm `dynamic_cast<const PlasticModel3DDeformationGradient*>` succeeds under this configuration.

## Current Scope Boundary

This walkthrough currently stops before:

- evaluating the total energy at a particular state `x`
- assembling the global force vector
- assembling the global Hessian for a state
- time integration and nonlinear solve
- collision and contact response
