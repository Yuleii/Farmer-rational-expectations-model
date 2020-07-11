import numpy as np
import pytest
from numpy.testing import assert_almost_equal
from function import solve_farmer_model


@pytest.fixture
def setup_solve_farmer_model():
    out = {}
    out["alpha"] = 0.3
    out["beta"] = 0.95
    out["delta"] = 0.1
    out["rho"] = 0.9
    out["T"] = 10
    out["inital_productivity_state"] = [0, -0.005, -0.009,
                                        -0.013, -0.022, -0.021, -0.019, -0.011, -0.012, -0.003]

    return out


@pytest.fixture
def expected_solve_farmer_model():
    out = {}
    out["s_t"] = np.array(
        [[0.], [-0.005], [-0.009], [-0.013], [-0.022], [-0.021], [-0.019],
            [-0.011], [-0.012], [-0.003]]
    )

    out["c_t"] = np.array(
        [[0.], [-0.00264342], [-0.00552647], [-0.00889971], [-0.01532715],
            [-0.01758031], [-0.01870046], [-0.01598844], [-0.01655944], [-0.01199041]]
    )
    out["k_t"] = np.array(
        [[0.], [0.], [-0.0014633], [-0.0038602], [-0.00703945], [-0.01233762],
            [-0.01648484], [-0.01937491], [-0.01945551], [-0.01981571], [-0.01748361]]
    )
    out["y_t"] = np.array(
        [[0.], [-0.005], [-0.00943899], [-0.01415806], [-0.02411184], [-0.02470129],
            [-0.02394545], [-0.01681247], [-0.01783665], [-0.00894471]]
    )

    return out


def test_solve_farmer_model_s_t(setup_solve_farmer_model, expected_solve_farmer_model):
    calc_s_t, calc_c_t, calc_k_t, calc_y_t = solve_farmer_model(**setup_solve_farmer_model)
    assert_almost_equal(calc_s_t, expected_solve_farmer_model["s_t"], decimal=6)


def test_solve_farmer_model_c_t(setup_solve_farmer_model, expected_solve_farmer_model):
    calc_s_t, calc_c_t, calc_k_t, calc_y_t = solve_farmer_model(**setup_solve_farmer_model)
    assert_almost_equal(calc_c_t, expected_solve_farmer_model["c_t"], decimal=6)


def test_solve_farmer_model_k_t(setup_solve_farmer_model, expected_solve_farmer_model):
    calc_s_t, calc_c_t, calc_k_t, calc_y_t = solve_farmer_model(**setup_solve_farmer_model)
    assert_almost_equal(calc_k_t, expected_solve_farmer_model["k_t"], decimal=6)


def test_solve_farmer_model_y_t(setup_solve_farmer_model, expected_solve_farmer_model):
    calc_s_t, calc_c_t, calc_k_t, calc_y_t = solve_farmer_model(**setup_solve_farmer_model)
    assert_almost_equal(calc_y_t, expected_solve_farmer_model["y_t"], decimal=6)
