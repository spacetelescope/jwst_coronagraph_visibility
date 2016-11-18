from .. import gui

def test_query_simbad():
    result = gui.query_simbad('HR 8799')
    eps = 0.5 / 60. / 60.   # 1/2 arcsec
    assert abs(result.ra - 346.86964613) < eps
    assert abs(result.dec - 21.13425148) < eps
    assert result.id == u'HD 218396'

    result = gui.query_simbad('sdfsdf')
    assert result is None
