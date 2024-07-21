import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal

from symd.groups import projectors2d, projectors3d, project_cell

def get_projectors(dimension):
    return projectors2d if dimension == 2 else projectors3d


@pytest.mark.parametrize("name,proj,dimension", [
    (name, proj, dimension) for dimension in [2, 3] for name, proj in get_projectors(dimension).items()
])
def test_projector_shape(name, proj, dimension):
    print(f"Testing {name} {proj} {dimension} projector")
    expected_shape = (dimension**2, dimension**2)
    assert proj.shape == expected_shape, f"Projector {name} has incorrect shape: {proj.shape}, expected: {expected_shape}"

@pytest.mark.parametrize("name,proj", [
    (name, proj) for dim in [2, 3] for name, proj in get_projectors(dim).items()
])

def test_projector_properties(name, proj):
    dimension = int(np.sqrt(proj.shape[0]))
    cells = [np.random.rand(dimension, dimension) for _ in range(100)]
    for cell in cells:
        projected_cell = project_cell(cell, proj)
        check_properties(projected_cell, name)
        check_idempotence(cell, proj)

def check_idempotence(cell, projector):
    return # Skip this test for now
    projected_once = project_cell(cell, projector)
    projected_twice = project_cell(projected_once, projector)
    assert_array_almost_equal(projected_once, projected_twice, err_msg=f"Projector failed idempotence check")

def check_properties(cell, lattice_type):
    property_checks = {
        'Square': check_square_properties,
        'Rectangular': check_rectangular_properties,
        'Hexagonal': check_hexagonal_properties,
        'Oblique': check_oblique_properties,
        'Cubic': check_cubic_properties,
        'Tetragonal': check_tetragonal_properties,
        'Orthorhombic': check_orthorhombic_properties,
        'Trigonal': check_trigonal_properties,
        'Monoclinic': check_monoclinic_properties,
        'Triclinic': check_triclinic_properties
    }
    check_function = property_checks.get(lattice_type, lambda x: None)
    check_function(cell)

def check_square_properties(cell):
    assert np.allclose(cell[0, 0], cell[1, 1]), "Square: a != b"
    assert np.allclose(cell[0, 1], 0) and np.allclose(cell[1, 0], 0), "Square: non-orthogonal"

def check_rectangular_properties(cell):
    assert np.allclose(cell[0, 1], 0) and np.allclose(cell[1, 0], 0), "Rectangular: non-orthogonal"

def check_hexagonal_properties(cell):
    dimension = int(np.sqrt(cell.shape[0]))
    assert np.allclose(np.linalg.norm(cell[:, 0]), np.linalg.norm(cell[:, 1])), "Hexagonal: a != b"
    angle = np.arccos(np.dot(cell[:, 0], cell[:, 1]) / (np.linalg.norm(cell[:, 0]) * np.linalg.norm(cell[:, 1])))
    expected_angle = np.pi / 3 if dimension == 2 else np.pi * 2 / 3
    assert np.allclose(angle, expected_angle, atol=1e-5), f"Hexagonal: incorrect angle (got {np.degrees(angle)})"

def check_oblique_properties(cell):
    # Oblique lattice has no specific constraints
    pass

def check_cubic_properties(cell):
    assert np.allclose(cell[0, 0], cell[1, 1]) and np.allclose(cell[1, 1], cell[2, 2]), "Cubic: a != b != c"
    assert_array_almost_equal(cell - np.diag(np.diag(cell)), np.zeros((3, 3)), err_msg="Cubic: non-orthogonal")

def check_tetragonal_properties(cell):
    assert np.allclose(cell[0, 0], cell[1, 1]), "Tetragonal: a != b"
    assert_array_almost_equal(cell - np.diag(np.diag(cell)), np.zeros((3, 3)), err_msg="Tetragonal: non-orthogonal")

def check_orthorhombic_properties(cell):
    assert_array_almost_equal(cell - np.diag(np.diag(cell)), np.zeros((3, 3)), err_msg="Orthorhombic: non-orthogonal")

def check_trigonal_properties(cell):
    angles = [np.arccos(np.dot(y, z) / (np.linalg.norm(y) * np.linalg.norm(z))) 
              for y, z in [(cell[:,1], cell[:,2]), (cell[:,0], cell[:,2]), (cell[:,0], cell[:,1])]]
    
    alpha, beta, gamma = np.degrees(angles)
    
    lengths = [np.linalg.norm(v) for v in cell.T]
    a, b, c = lengths

    assert np.isclose(a, b), f"a ({a:.4f}) is not equal to b ({b:.4f})"
    assert np.isclose(alpha, beta), f"Alpha ({alpha:.2f}°) is not equal to Beta ({beta:.2f}°)"
    assert np.isclose(alpha, gamma), f"Alpha ({alpha:.2f}°) is not equal to Gamma ({gamma:.2f}°)"
    assert not np.isclose(alpha, 90), f"All angles are 90°: Alpha = {alpha:.2f}°"

def check_monoclinic_properties(cell):
    angles = [np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))) 
              for v1, v2 in [(cell[:,1], cell[:,2]), (cell[:,0], cell[:,2]), (cell[:,0], cell[:,1])]]
    
    alpha, beta, gamma = np.degrees(angles)
    print(alpha, beta, gamma)
    
    assert np.isclose(alpha, 90), f"Alpha angle is not 90°: {alpha:.2f}°"
    assert np.isclose(gamma, 90), f"Gamma angle is not 90°: {gamma:.2f}°"
    assert not np.isclose(beta, 90), f"Beta angle is 90°: {beta:.2f}°"
    
    return True
def check_triclinic_properties(cell):
    # Triclinic lattice has no specific constraints
    pass

# Run the tests
if __name__ == "__main__":
    pytest.main([__file__])