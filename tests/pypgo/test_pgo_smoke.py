from pathlib import Path

import numpy as np
import pytest

import pypgo


EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples" / "legacy" / "tet"
TORUS_VEG = EXAMPLES_DIR / "torus.veg"


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
