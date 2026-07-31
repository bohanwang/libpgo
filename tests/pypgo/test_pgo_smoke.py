import json
import math
from pathlib import Path
import tomllib

import numpy as np
import pytest

import pypgo


EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples"
TORUS_VEG = EXAMPLES_DIR / "assets" / "common" / "torus-tet.veg"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SIM_ASSETS_DIR = Path(__file__).with_name("assets")

TET_VERTICES = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ],
    dtype=np.float32,
)
TET_ELEMENTS = np.array([[0, 1, 2, 3]], dtype=np.int32)

TRI_VERTICES = TET_VERTICES[:3].copy()
TRI_ELEMENTS = np.array([[0, 1, 2]], dtype=np.int32)


def _load_torus_mesh():
    tetmesh = pypgo.create_tetmeshgeo_from_file(str(TORUS_VEG))
    try:
        vertices = pypgo.tetmeshgeo_get_vertices(tetmesh)
        tets = pypgo.tetmeshgeo_get_tets(tetmesh)
    finally:
        pypgo.destroy_tetmeshgeo(tetmesh)

    return np.reshape(vertices, (-1, 3)), np.reshape(tets, (-1, 4))


def _build_repeated_mesh(vertices: np.ndarray, tets: np.ndarray, repeats: int = 10):
    repeated_vertices = np.concatenate([vertices for _ in range(repeats)], axis=0)
    repeated_tets = np.concatenate(
        [tets + repeat_index * len(vertices) for repeat_index in range(repeats)],
        axis=0,
    )

    tetmesh = pypgo.create_tetmeshgeo(
        repeated_vertices.flatten().astype(np.float32),
        repeated_tets.flatten().astype(np.int32),
    )
    return tetmesh, repeated_vertices, repeated_tets


def _create_unit_tetmeshgeo():
    return pypgo.create_tetmeshgeo(TET_VERTICES.ravel(), TET_ELEMENTS.ravel())


def _create_unit_trimeshgeo():
    return pypgo.create_trimeshgeo(TRI_VERTICES.ravel(), TRI_ELEMENTS.ravel())


def _simulation_config(
    mesh_type: str,
    contact_model: str,
    sim_type: str,
    output_path: Path,
) -> dict:
    config = {
        f"{mesh_type}-mesh": str(SIM_ASSETS_DIR / f"box-{mesh_type}.veg"),
        "surface-mesh": str(SIM_ASSETS_DIR / "box-surface.obj"),
        "fixed-vertices": [
            {
                "filename": str(SIM_ASSETS_DIR / "box-fixed.txt"),
                "movement": [0.0, 0.0, 0.0],
                "coeff": 1.0e5,
            }
        ],
        "g": [0.0, -0.1, 0.0],
        "init-vel": [0.0, 0.0, 0.0],
        "init-disp": [0.0, 0.0, 0.0],
        "scale": 1.0,
        "timestep": 0.001,
        "num-timestep": 2 if sim_type == "dynamic" else 1,
        "damping-params": [0.0, 0.0],
        "sim-type": sim_type,
        "solver-eps": 1.0e-4,
        "solver-max-iter": 5,
        "elastic-material": "stable-neo",
        "dump-interval": 1,
        "output": str(output_path),
        "contact-model": contact_model,
    }
    if contact_model == "ipc":
        config.update(
            {
                "ipc-dhat": 0.002,
                "ipc-kappa": 3000.0,
                "use-floor": True,
                "floor-axis": "y",
                "floor-height": 0.0,
                "floor-kappa": 3000.0,
            }
        )
    else:
        config.update(
            {
                "contact-stiffness": 1000.0,
                "contact-sample": 1,
                "contact-friction-coeff": 0.0,
                "contact-vel-eps": 1.0e-5,
            }
        )
    return config


def _assert_finite_obj(path: Path):
    vertices = []
    for line in path.read_text().splitlines():
        if line.startswith("v "):
            vertices.append([float(value) for value in line.split()[1:4]])
    assert vertices
    assert all(math.isfinite(value) for vertex in vertices for value in vertex)


def test_package_version_matches_project_metadata():
    with (PROJECT_ROOT / "pyproject.toml").open("rb") as pyproject_file:
        project_version = tomllib.load(pyproject_file)["project"]["version"]

    assert pypgo.__version__ == project_version


@pytest.mark.parametrize("mesh_type", ["tet", "cubic"])
@pytest.mark.parametrize("contact_model", ["sampled", "ipc"])
@pytest.mark.parametrize("sim_type", ["dynamic", "static"])
def test_python_simulation_entry_point_smoke(
    mesh_type: str,
    contact_model: str,
    sim_type: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capfd: pytest.CaptureFixture[str],
):
    case_name = f"{mesh_type}-{contact_model}-{sim_type}"
    output_path = (
        tmp_path / f"{case_name}.obj"
        if contact_model == "sampled" and sim_type == "static"
        else tmp_path / case_name
    )
    config_path = tmp_path / f"{case_name}.json"
    config_path.write_text(
        json.dumps(
            _simulation_config(mesh_type, contact_model, sim_type, output_path),
            indent=2,
        )
        + "\n"
    )
    # The sampled runner writes a diagnostic fv.obj in its current directory.
    # Keep that legacy side effect inside pytest's temporary directory.
    monkeypatch.chdir(tmp_path)

    assert pypgo.run_sim_from_config(str(config_path)) == 0
    captured_output = capfd.readouterr()
    if contact_model == "sampled" and sim_type == "static":
        assert "Static sampled contact is not enforced" in (
            captured_output.out + captured_output.err
        )

    if contact_model == "sampled" and sim_type == "static":
        surface_path = output_path
    else:
        surface_path = output_path / "ret0000.obj"
        state_path = output_path / "deform0000.u"
        assert state_path.is_file()
        assert state_path.stat().st_size > 0
    assert surface_path.is_file()
    _assert_finite_obj(surface_path)


def test_tetmeshgeo_memory_round_trip_and_barycentric_query():
    tetmesh = _create_unit_tetmeshgeo()
    try:
        assert pypgo.tetmeshgeo_get_num_vertices(tetmesh) == 4
        assert pypgo.tetmeshgeo_get_num_tets(tetmesh) == 1
        np.testing.assert_allclose(
            pypgo.tetmeshgeo_get_vertices(tetmesh).reshape(-1, 3),
            TET_VERTICES,
        )
        np.testing.assert_array_equal(
            pypgo.tetmeshgeo_get_tets(tetmesh).reshape(-1, 4),
            TET_ELEMENTS,
        )

        weights, element_indices = pypgo.tetmesh_barycentric_weights(
            tetmesh,
            np.array([0.25, 0.25, 0.25], dtype=np.float32),
        )
        np.testing.assert_allclose(weights.reshape(-1, 4), [[0.25] * 4])
        np.testing.assert_array_equal(element_indices, [0])
    finally:
        pypgo.destroy_tetmeshgeo(tetmesh)


def test_tetmesh_update_and_file_round_trip(tmp_path: Path):
    tetmesh = pypgo.create_tetmesh(
        TET_VERTICES.ravel(),
        TET_ELEMENTS.ravel(),
        1.0e5,
        0.45,
        1000.0,
    )
    translated_vertices = TET_VERTICES + np.array([1.0, 2.0, 3.0], dtype=np.float32)
    updated_tetmesh = pypgo.update_tetmesh_vertices(tetmesh, translated_vertices)
    output_path = tmp_path / "unit-tet.veg"

    try:
        assert translated_vertices.shape == (4, 3)
        np.testing.assert_allclose(
            pypgo.get_tetmesh_vertex_positions(updated_tetmesh),
            translated_vertices,
        )
        np.testing.assert_array_equal(
            pypgo.get_tetmesh_element_indices(updated_tetmesh),
            TET_ELEMENTS,
        )
        pypgo.save_tetmesh_to_file(updated_tetmesh, str(output_path))

        loaded_tetmesh = pypgo.create_tetmesh_from_file(str(output_path))
        try:
            np.testing.assert_allclose(
                pypgo.get_tetmesh_vertex_positions(loaded_tetmesh),
                translated_vertices,
            )
            np.testing.assert_array_equal(
                pypgo.get_tetmesh_element_indices(loaded_tetmesh),
                TET_ELEMENTS,
            )
        finally:
            pypgo.destroy_tetmesh(loaded_tetmesh)
    finally:
        pypgo.destroy_tetmesh(updated_tetmesh)
        pypgo.destroy_tetmesh(tetmesh)


def test_tetmesh_update_rejects_wrong_vertex_count():
    tetmesh = pypgo.create_tetmesh(
        TET_VERTICES.ravel(),
        TET_ELEMENTS.ravel(),
        1.0e5,
        0.45,
        1000.0,
    )
    try:
        with pytest.raises(ValueError, match="row count"):
            pypgo.update_tetmesh_vertices(tetmesh, TET_VERTICES[:3])
    finally:
        pypgo.destroy_tetmesh(tetmesh)


def test_tetmesh_creation_forcecasts_integer_elements():
    tetmesh = pypgo.create_tetmesh(
        TET_VERTICES.ravel(),
        TET_ELEMENTS.astype(np.int64).ravel(),
        1.0e5,
        0.45,
        1000.0,
    )
    try:
        np.testing.assert_array_equal(
            pypgo.get_tetmesh_element_indices(tetmesh),
            TET_ELEMENTS,
        )
    finally:
        pypgo.destroy_tetmesh(tetmesh)


def test_mesh_creation_rejects_malformed_flat_arrays():
    with pytest.raises(ValueError, match="multiple of 3"):
        pypgo.create_tetmeshgeo(
            TET_VERTICES.ravel()[:-1],
            TET_ELEMENTS.ravel(),
        )

    with pytest.raises(ValueError, match="multiple of 4"):
        pypgo.create_tetmesh(
            TET_VERTICES.ravel(),
            TET_ELEMENTS.ravel()[:-1],
            1.0e5,
            0.45,
            1000.0,
        )

    with pytest.raises(ValueError, match="multiple of 3"):
        pypgo.create_trimeshgeo(
            TRI_VERTICES.ravel(),
            TRI_ELEMENTS.ravel()[:-1],
        )


def test_trimeshgeo_memory_round_trip_and_closest_distance():
    trimesh = _create_unit_trimeshgeo()
    try:
        assert pypgo.trimeshgeo_get_num_vertices(trimesh) == 3
        assert pypgo.trimeshgeo_get_num_triangles(trimesh) == 1
        np.testing.assert_allclose(
            pypgo.trimeshgeo_get_vertices(trimesh).reshape(-1, 3),
            TRI_VERTICES,
        )
        np.testing.assert_array_equal(
            pypgo.trimeshgeo_get_triangles(trimesh).reshape(-1, 3),
            TRI_ELEMENTS,
        )

        squared_distances = pypgo.trimesh_closest_distances(
            trimesh,
            np.array([0.25, 0.25, 2.0, 0.25, 0.25, 0.0], dtype=np.float32),
        )
        np.testing.assert_allclose(squared_distances, [4.0, 0.0])
    finally:
        pypgo.destroy_trimeshgeo(trimesh)


def test_tet_gradient_sparse_matrix_is_finite():
    tetmesh = _create_unit_tetmeshgeo()
    try:
        gradient = pypgo.create_tet_gradient_matrix(tetmesh)
        try:
            rows = pypgo.sparse_matrix_get_row_indices(gradient)
            cols = pypgo.sparse_matrix_get_col_indices(gradient)
            values = pypgo.sparse_matrix_get_values(gradient)

            assert pypgo.sparse_matrix_get_num_entries(gradient) > 0
            assert rows.shape == cols.shape == values.shape
            assert np.isfinite(values).all()
            assert np.max(np.abs(values)) > 0
        finally:
            pypgo.destroy_sparse_matrix(gradient)
    finally:
        pypgo.destroy_tetmeshgeo(tetmesh)


@pytest.mark.parametrize("path", [TORUS_VEG])
def test_can_load_tetmesh_from_example_file(path: Path):
    assert path.is_file(), f"missing test fixture: {path}"

    vertices, tets = _load_torus_mesh()

    assert vertices.ndim == 2
    assert vertices.shape[1] == 3
    assert vertices.shape[0] > 0
    assert np.isfinite(vertices).all()

    assert tets.ndim == 2
    assert tets.shape[1] == 4
    assert tets.shape[0] > 0
    assert np.issubdtype(tets.dtype, np.integer)


def test_python_api_smoke_workflow():
    vertices, tets = _load_torus_mesh()
    tetmesh, repeated_vertices, repeated_tets = _build_repeated_mesh(vertices, tets)

    try:
        laplacian = pypgo.create_element_laplacian_matrix(tetmesh, 0, 9, 0)
        try:
            row_indices = pypgo.sparse_matrix_get_row_indices(laplacian)
            col_indices = pypgo.sparse_matrix_get_col_indices(laplacian)
            values = pypgo.sparse_matrix_get_values(laplacian)

            assert pypgo.sparse_matrix_get_num_entries(laplacian) > 0
            assert row_indices.shape == col_indices.shape == values.shape
            assert row_indices.ndim == 1
            assert np.isfinite(values).all()
            assert np.max(values) > 0
            quadratic_form = pypgo.conjugate_mv(
                laplacian,
                np.ones(repeated_tets.shape[0] * 9, dtype=np.float32),
            )
            assert np.isfinite(quadratic_form)
        finally:
            pypgo.destroy_sparse_matrix(laplacian)

        biharmonic = pypgo.create_tet_biharmonic_gradient_matrix(tetmesh, 1, 0)
        try:
            row_indices = pypgo.sparse_matrix_get_row_indices(biharmonic)
            col_indices = pypgo.sparse_matrix_get_col_indices(biharmonic)
            values = pypgo.sparse_matrix_get_values(biharmonic)

            assert pypgo.sparse_matrix_get_num_entries(biharmonic) > 0
            assert row_indices.shape == col_indices.shape == values.shape
            assert row_indices.ndim == 1
            assert np.isfinite(values).all()
            assert np.max(values) > 0
        finally:
            pypgo.destroy_sparse_matrix(biharmonic)

        gradient_mats = pypgo.create_tet_gradient_per_element_matrix(tetmesh)
        assert gradient_mats.ndim == 3
        assert gradient_mats.shape[0] == repeated_tets.shape[0]
        # The C API fills a 9x12 element gradient matrix in Eigen's column-major
        # storage, and the Python binding exposes that buffer as a C-contiguous
        # (12, 9) array for each tet. Transposing recovers the mathematical 9x12
        # element gradient used throughout libpgo.
        assert gradient_mats.shape[1:] == (12, 9)
        assert np.isfinite(gradient_mats).all()

        first_gradient = np.transpose(gradient_mats[0, :, :])
        assert first_gradient.shape == (9, 12)
        assert np.isfinite(first_gradient).all()
        assert np.max(np.abs(first_gradient)) > 0
        assert repeated_vertices.shape[0] == vertices.shape[0] * 10
    finally:
        pypgo.destroy_tetmeshgeo(tetmesh)
