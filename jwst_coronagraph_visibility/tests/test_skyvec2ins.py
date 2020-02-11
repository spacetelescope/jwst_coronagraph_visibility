from os.path import abspath, dirname, join, isdir
import numpy as np
import datetime

from .. import skyvec2ins
from ..gui import get_aperture

TARGETS_DIR = abspath(join(dirname(__file__), 'targets'))
START_DATE = datetime.datetime(2018, 10, 1)
NPOINTS = 360
NROLLS = 20
MAXVROLL = 10.0

def _save_test_case(test_case_name, aperture,
                    ra, dec, pa1, pa2, pa3,
                    separation_as1, separation_as2, separation_as3):
    case_path = join(TARGETS_DIR, test_case_name)
    arrnames = (
        'x',
        'observable',
        'elongation_rad',
        'roll_rad',
        'c1_x', 'c1_y',
        'c2_x', 'c2_y',
        'c3_x', 'c3_y',
        'n_x', 'n_y',
        'e_x', 'e_y'
    )


    computed = skyvec2ins.skyvec2ins(
        ra=ra,
        dec=dec,
        pa1=pa1,
        separation_as1=separation_as1,
        pa2=pa2,
        separation_as2=separation_as2,
        pa3=pa3,
        separation_as3=separation_as3,
        aper=aperture,
        start_date=START_DATE,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    for name, arr in zip(arrnames, computed):
        outpath = join(case_path, '{}.csv'.format(name))
        np.savetxt(outpath, arr, delimiter=',')
        print('Saved', outpath)

def _generate_test_outputs():
    # Fomalhaut
    _save_test_case(
        'Fomalhaut',
        get_aperture('NIRCam', 'NRCA2_MASK210R'),
        ra=344.41269,
        dec=-29.62224,
        pa1=325,
        pa2=0,
        pa3=0,
        separation_as1=10,
        separation_as2=0,
        separation_as3=0,
    )
    # 1RXSJ160929p1-210524
    _save_test_case(
        '1RXSJ160929p1-210524',
        get_aperture('NIRCam', 'NRCB3_MASKSWB'),
        ra=242.37628,
        dec=-21.08304,
        pa1=20,
        pa2=0,
        pa3=0,
        separation_as1=3,
        separation_as2=0,
        separation_as3=0,
    )
    # HR8799
    _save_test_case(
        'HR8799',
        get_aperture('MIRI', 'MIRIM_MASK1065'),
        ra=346.86965,
        dec=21.13425,
        pa1=45,
        separation_as1=1.7,
        pa2=325,
        separation_as2=1,
        pa3=190,
        separation_as3=0.65,
    )
    # NGC 6543
    _save_test_case(
        'NGC6543',
        get_aperture('MIRI', 'MIRIM_MASKLYOT'),
        ra=269.63926,
        dec=66.63320,
        pa1=0,
        separation_as1=0,
        pa2=0,
        separation_as2=0,
        pa3=0,
        separation_as3=0,
    )

def _load_test_case(test_case_name):
    case_path = join(TARGETS_DIR, test_case_name)
    assert isdir(case_path)
    arrs = (
        'x',
        'observable',
        'elongation_rad',
        'roll_rad',
        'c1_x', 'c1_y',
        'c2_x', 'c2_y',
        'c3_x', 'c3_y',
        'n_x', 'n_y',
        'e_x', 'e_y'
    )
    return (np.genfromtxt(join(case_path, '{}.csv'.format(n)), delimiter=',') for n in arrs)

def _compare_outputs(reference, computed):
    (
        x,
        observable,
        elongation_rad,
        roll_rad,
        c1_x, c1_y,
        c2_x, c2_y,
        c3_x, c3_y,
        n_x, n_y,
        e_x, e_y
    ) = reference

    (
        t_x,
        t_observable,
        t_elongation_rad,
        t_roll_rad,
        t_c1_x, t_c1_y,
        t_c2_x, t_c2_y,
        t_c3_x, t_c3_y,
        t_n_x, t_n_y,
        t_e_x, t_e_y
    ) = computed

    assert np.allclose(x, t_x)
    assert np.allclose(elongation_rad, t_elongation_rad)
    assert np.allclose(roll_rad, t_roll_rad, atol=2e-6)
    assert not np.any((observable == 1) ^ (t_observable == 1))


    nircam_pixelscale = 0.0311 # for short-wavelen channels, SIAF PRDDEVSOC-D-012, 2016 April
    siaf_transform_epsilon = nircam_pixelscale / 100
    # rationale: comparison of the SIAF transforms shows they should be
    # mathematically correct in both implementations, but numerical errors are
    # somehow being compounded to result in errors that are nevertheless small
    # relative to the size of a pixel (<< 0.01 px). We set the tolerance at
    # 1/100 of a NIRCam pixel.

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
    aperture = get_aperture('NIRCam', 'NRCA2_MASK210R')
    computed = skyvec2ins.skyvec2ins(
        ra=344.41269,
        dec=-29.62224,
        pa1=325,
        pa2=0,
        pa3=0,
        separation_as1=10,
        separation_as2=0,
        separation_as3=0,
        aper=aperture,
        start_date=START_DATE,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)

def test_1RXSJ160929p1_210524():
    reference = _load_test_case('1RXSJ160929p1-210524')
    aperture = get_aperture('NIRCam', 'NRCB3_MASKSWB')
    computed = skyvec2ins.skyvec2ins(
        ra=242.37628,
        dec=-21.08304,
        pa1=20,
        pa2=0,
        pa3=0,
        separation_as1=3,
        separation_as2=0,
        separation_as3=0,
        aper=aperture,
        start_date=START_DATE,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)

def test_HR8799():
    reference = _load_test_case('HR8799')
    aperture = get_aperture('MIRI', 'MIRIM_MASK1065')
    computed = skyvec2ins.skyvec2ins(
        ra=346.86965,
        dec=21.13425,
        pa1=45,
        separation_as1=1.7,
        pa2=325,
        separation_as2=1,
        pa3=190,
        separation_as3=0.65,
        aper=aperture,
        start_date=START_DATE,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)

def test_NGC6543():
    reference = _load_test_case('NGC6543')
    aperture = get_aperture('MIRI', 'MIRIM_MASKLYOT')
    computed = skyvec2ins.skyvec2ins(
        ra=269.63926,
        dec=66.63320,
        pa1=0,
        separation_as1=0,
        pa2=0,
        separation_as2=0,
        pa3=0,
        separation_as3=0,
        aper=aperture,
        start_date=START_DATE,
        npoints=NPOINTS,
        nrolls=NROLLS,
        maxvroll=MAXVROLL,
    )
    _compare_outputs(reference, computed)
