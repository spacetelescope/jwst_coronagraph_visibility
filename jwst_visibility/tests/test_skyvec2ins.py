from os.path import abspath, dirname, join, isdir
import numpy as np

from .. import skyvec2ins

NPOINTS = 360
NROLLS = 20
MAXVROLL = 10.0
LAMBDA_RAD0 = np.deg2rad(191.05)

def _load_test_case(test_case_name):
    case_path = abspath(join(dirname(__file__), 'targets', test_case_name))
    assert isdir(case_path)
    arrs = (
        'x',
        'observable',
        'elongation_rad',
        'roll_rad',
        's_x', 's_y',
        'c1_x', 'c1_y',
        'c2_x', 'c2_y',
        'c3_x', 'c3_y',
        'scisize', 'sciscale', 'sciyangle',
        'n_x', 'n_y', 'e_x', 'e_y'
    )
    return (np.genfromtxt(join(case_path, '{}.csv'.format(n)), delimiter=',') for n in arrs)

def _compare_outputs(reference, computed):
    (
        x,
        observable,
        elongation_rad,
        roll_rad,
        s_x, s_y,
        c1_x, c1_y,
        c2_x, c2_y,
        c3_x, c3_y,
        scisize, sciscale, sciyangle,
        n_x, n_y, e_x, e_y
    ) = reference

    (
        t_x,
        t_observable,
        t_elongation_rad,
        t_roll_rad,
        t_s_x, t_s_y,
        t_c1_x, t_c1_y,
        t_c2_x, t_c2_y,
        t_c3_x, t_c3_y,
        t_scisize, t_sciscale, t_sciyangle,
        t_n_x, t_n_y, t_e_x, t_e_y
    ) = computed

    assert np.allclose(x, t_x)
    assert np.allclose(elongation_rad, t_elongation_rad)
    assert np.allclose(roll_rad, t_roll_rad, atol=2e-6)
    assert not np.any((observable == 1) ^ (t_observable == 1))

    assert np.allclose(scisize, t_scisize, atol=1e-7)
    assert np.allclose(sciscale, t_sciscale, atol=1e-7)
    assert np.allclose(sciyangle, t_sciyangle)

    assert np.allclose(s_x, t_s_x, atol=8e-4)  # TODO verify with cstark
    assert np.allclose(s_y, t_s_y, atol=8e-4)  # why this is large

    siaf_transform_epsilon = 13e-4
    # rationale: comparison of the SIAF transforms shows they should be
    # mathematically correct in both implementations, but numerical errors are
    # somehow being compounded to result in errors that are nevertheless small
    # relative to the size of a pixel (<< 0.01 px).

    # n.b. the residuals are larger in Y for this test case
    # see https://github.com/mperrin/jwxml/issues/4

    assert np.allclose(c1_x, t_c1_x, atol=siaf_transform_epsilon)
    assert np.allclose(c1_y, t_c1_y, atol=siaf_transform_epsilon)
    assert np.allclose(c2_x, t_c2_x, atol=siaf_transform_epsilon)
    assert np.allclose(c2_y, t_c2_y, atol=siaf_transform_epsilon)
    assert np.allclose(c3_x, t_c3_x, atol=siaf_transform_epsilon)
    assert np.allclose(c3_y, t_c3_y, atol=siaf_transform_epsilon)
    assert np.allclose( n_x,  t_n_x, atol=siaf_transform_epsilon)
    assert np.allclose( n_y,  t_n_y, atol=siaf_transform_epsilon)
    assert np.allclose( e_x,  t_e_x, atol=siaf_transform_epsilon)
    assert np.allclose( e_y,  t_e_y, atol=siaf_transform_epsilon)

def test_fomalhaut():
    reference = _load_test_case('Fomalhaut')
    computed = skyvec2ins.skyvec2ins(
        ra=344.41269,
        dec=-29.62224,
        pa1=325,
        pa2=0,
        pa3=0,
        separation_as1=10,
        separation_as2=0,
        separation_as3=0,
        instrname='NIRCam',
        apername='NRCA2_MASK210R',
        lambda_rad0=LAMBDA_RAD0,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)

def test_1RXSJ160929p1_210524():
    reference = _load_test_case('1RXSJ160929p1-210524')
    computed = skyvec2ins.skyvec2ins(
        ra=242.37628,
        dec=-21.08304,
        pa1=20,
        pa2=0,
        pa3=0,
        separation_as1=3,
        separation_as2=0,
        separation_as3=0,
        instrname='NIRCam',
        apername='NRCB3_MASKSWB',
        lambda_rad0=LAMBDA_RAD0,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)

def test_HR8799():
    reference = _load_test_case('HR8799')
    computed = skyvec2ins.skyvec2ins(
        ra=346.86965,
        dec=21.13425,
        pa1=45,
        separation_as1=1.7,
        pa2=325,
        separation_as2=1,
        pa3=190,
        separation_as3=0.65,
        instrname='MIRI',
        apername='MIRIM_MASK1065',
        lambda_rad0=LAMBDA_RAD0,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)

def test_NGC6543():
    reference = _load_test_case('NGC6543')
    computed = skyvec2ins.skyvec2ins(
        ra=269.63926,
        dec=66.63320,
        pa1=0,
        separation_as1=0,
        pa2=0,
        separation_as2=0,
        pa3=0,
        separation_as3=0,
        instrname='MIRI',
        apername='MIRIM_MASKLYOT',
        lambda_rad0=LAMBDA_RAD0,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)
