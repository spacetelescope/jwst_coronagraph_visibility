#!/usr/bin/env python
# vim: set fileencoding=utf8 :
from __future__ import print_function, division
import sys
try:
    from tkinter import *
    from tkinter import ttk
except ImportError:
    from Tkinter import *
    import ttk
import os
import os.path
import datetime
import re
from collections import namedtuple
from contextlib import contextmanager

try:
    from urllib import quote
except ImportError:
    from urllib.parse import quote

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import patches
from matplotlib import pyplot as plt
plt.style.use('ggplot')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import numpy as np
import requests
import requests.exceptions
if getattr(sys, 'frozen', False):
    # we are running in a bundle
    bundle_dir = sys._MEIPASS
else:
    # we are running in a normal Python environment
    bundle_dir = os.path.dirname(os.path.abspath(__file__))

SimbadResult = namedtuple('SimbadResult', ['ra', 'dec', 'id'])

from jwxml import SIAF
import jwxml
from .skyvec2ins import skyvec2ins, ad2lb, lb2ad

from pprint import pprint

RED_GGPLOT = '#E24A33'
BLUE_GGPLOT = '#348ABD'
PURPLE_GGPLOT = '#988ED5'
GRAY_GGPLOT = '#777777'
YELLOW_GGPLOT = '#FBC15E'
GREEN_GGPLOT = '#8EBA42'
PINK_GGPLOT = '#FFB5B8'

QUERY_TIMEOUT_SEC = 1.0

DEFAULT_NPOINTS = 360
DEFAULT_NROLLS = 20

# Outlining the 'bad' areas of the NIRCam Module A coronagraphs requires some
# coordinate conversion gymnastics as there is an optical wedge in the pupil
# wheel that changes the angular to pixel transformation

def compute_v2v3_offset(aperture_a, aperture_b):
    '''
    For the same pixel coordinates, different V2, V3 coordinates are used
    depending on whether the coronagraph pupil wheel wedge is in the beam.
    The offset is computed by transforming the same pixel (Det) coordinates
    to V2, V3 in two different apertures and computing the difference in
    the resulting Tel frame coordinates
    '''
    x_a, y_a = aperture_a.Det2Tel(aperture_b.XDetRef,  aperture_b.YDetRef)
    x_b, y_b = aperture_b.Det2Tel(aperture_b.XDetRef,  aperture_b.YDetRef)
    return x_a - x_b, y_a - y_b

_NIRCAM_SIAF = SIAF('NIRCam')
_NIRCAM_CORON_OFFSET_TEL = compute_v2v3_offset(
    _NIRCAM_SIAF['NRCA5_MASKLWB'],
    _NIRCAM_SIAF['NRCA5_FULL']
)

# These bad areas were defined in raw detector pixel coordinates by
# John Stansberry using a backlit image of NIRCam A5 through the
# long-wavelength bar coronagraph pupil wedge. Colin Cox translated them to
# V2, V3 for A5

NIRCAM_CORON_BAD_AREAS = np.array(
       [[[  56.525 , -462.2092],
        [  56.5505, -456.9293],
        [  61.8066, -456.9541],
        [  61.7856, -462.2317]],

       [[  38.8981, -462.1289],
        [  38.9384, -456.8408],
        [  44.2083, -456.8827],
        [  44.1725, -462.1682]],

       [[  96.9128, -462.0046],
        [  96.9041, -456.7388],
        [ 102.1434, -456.7248],
        [ 102.1565, -461.9892]],

       [[ 117.1985, -461.9805],
        [ 117.1727, -456.7191],
        [ 122.4112, -456.6854],
        [ 122.4413, -461.9459]],

       [[  76.8615, -462.1433],
        [  76.8697, -456.8714],
        [  82.1149, -456.8766],
        [  82.1112, -462.1468]],

       [[ 134.8664, -462.0006],
        [ 134.8259, -456.7415],
        [ 140.0682, -456.6904],
        [ 140.113 , -461.949 ]],

       [[  38.9735, -444.1376],
        [  38.9862, -442.5039],
        [  41.2405, -442.5257],
        [  41.2284, -444.159 ]],

       [[  58.114 , -444.162 ],
        [  58.1216, -442.6566],
        [  60.3693, -442.6698],
        [  60.3623, -444.1749]],

       [[  78.2656, -444.1798],
        [  78.2685, -442.7392],
        [  80.5116, -442.7435],
        [  80.5091, -444.1839]],

       [[  98.3837, -443.9939],
        [  98.3826, -442.7426],
        [ 100.6232, -442.7381],
        [ 100.6248, -443.9892]],

       [[ 118.4879, -443.9169],
        [ 118.483 , -442.6667],
        [ 120.7234, -442.6532],
        [ 120.7286, -443.9034]],

       [[ 137.5412, -444.0201],
        [ 137.5313, -442.5203],
        [ 139.7735, -442.4982],
        [ 139.7839, -443.9979]],

       [[  21.83  , -463.4815],
        [  22.1856, -428.4652],
        [  38.844 , -428.6853],
        [  38.5725, -463.6381]],

       [[ 140.7495, -463.2577],
        [ 140.4967, -428.7508],
        [ 149.7032, -428.6478],
        [ 150.0012, -463.153 ]],

       [[  38.8609, -442.5026],
        [  38.9691, -428.6867],
        [ 140.4967, -428.7508],
        [ 140.5833, -442.49  ]],

       [[  38.1157, -465.8399],
        [  38.136 , -463.1937],
        [ 141.1852, -463.0652],
        [ 141.2094, -465.6963]]])
NIRCAM_CORON_BAD_AREAS[:,:,0] += _NIRCAM_CORON_OFFSET_TEL[0]
NIRCAM_CORON_BAD_AREAS[:,:,1] += _NIRCAM_CORON_OFFSET_TEL[1]
NIRCAM_CORON_BAD_AREAS.flags.writeable = False

def query_simbad(query_string):
    response = requests.get('http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oI?' + quote(query_string), timeout=QUERY_TIMEOUT_SEC)
    body = response.text
    ra = dec = canonical_id = None
    for line in body.split('\n'):
        if line[:2] == '%J' and ra is None:
            match = re.match('%J (\d+\.\d+) ([+\-]\d+\.\d+) .+', line)
            if match is None:
                return None
            ra, dec = map(float, match.groups())
        elif line[:4] == '%I.0' and canonical_id is None:
            match = re.match('%I.0 (.+)', line)
            if match is None:
                return None
            canonical_id = match.groups()[0]
    if ra is None or canonical_id is None:
        return None
    else:
        return SimbadResult(ra=ra, dec=dec, id=canonical_id)

def get_aperture(instrname, apername):
    # siaf_path = os.path.join(bundle_dir, 'data', '{}_SIAF.xml'.format(instrname))
    # assert os.path.exists(siaf_path), 'no SIAF for {} at {}'.format(instrname, siaf_path)
    siaf = SIAF(instr=instrname)
    return siaf[apername]

@contextmanager
def _busy_cursor(root):
    try:
        root.config(cursor='wait')
    except TclError:
        pass
    root.update()
    yield
    root.config(cursor='')
    root.update()

class VisibilityCalculation(object):
    def __init__(self, ra, dec, companions, aperture, start_date, npoints, nrolls):
        self.ra = ra
        self.dec = dec
        self.companions = companions
        self.aperture = aperture
        self.npoints = npoints
        self.nrolls = nrolls
        self.start_date = start_date

        # Outputs
        self.days = None
        self.observable = None
        self.elongation_rad = None
        self.roll_rad = None
        self.c1_x = None
        self.c1_y = None
        self.c2_x = None
        self.c2_y = None
        self.c3_x = None
        self.c3_y = None
        self.n_x = None
        self.n_y = None
        self.e_x = None
        self.e_y = None

    def calculate(self):
        (
            self.days,
            self.observable,
            self.elongation_rad,
            self.roll_rad,
            self.c1_x, self.c1_y,
            self.c2_x, self.c2_y,
            self.c3_x, self.c3_y,
            self.n_x, self.n_y,
            self.e_x, self.e_y
        ) = skyvec2ins(
            ra=self.ra,
            dec=self.dec,
            pa1=self.companions[0]['pa'],
            pa2=self.companions[1]['pa'],
            pa3=self.companions[2]['pa'],
            separation_as1=self.companions[0]['separation'],
            separation_as2=self.companions[1]['separation'],
            separation_as3=self.companions[2]['separation'],
            aper=self.aperture,
            start_date=self.start_date,
            npoints=self.npoints,
            nrolls=self.nrolls
        )

        # mask non-observable (roll, elongation) pairs from output data
        mask = self.observable == 0

        self.c1_x[mask] = np.nan
        self.c1_y[mask] = np.nan
        self.c2_x[mask] = np.nan
        self.c2_y[mask] = np.nan
        self.c3_x[mask] = np.nan
        self.c3_y[mask] = np.nan

        self.n_x[mask] = np.nan
        self.n_y[mask] = np.nan
        self.e_x[mask] = np.nan
        self.e_y[mask] = np.nan

class VisibilityCalculator(object):
    NIRCAM_A = 'NIRCam Channel A'
    NIRCAM_B = 'NIRCam Channel B'
    MIRI = 'MIRI'
    INSTRUMENTS = [NIRCAM_A, MIRI]
    NIRCAM_A_APERNAMES = [
        'NRCA2_MASK210R',
        'NRCA5_MASK335R',
        'NRCA5_MASK430R',
        'NRCA4_MASKSWB',
        'NRCA5_MASKLWB'
    ]
    NIRCAM_B_APERNAMES = [
        'NRCB1_MASK210R',
        'NRCB5_MASK335R',
        'NRCB5_MASK430R',
        'NRCB3_MASKSWB',
        'NRCB5_MASKLWB',
    ]
    MIRI_APERNAMES = [
        'MIRIM_CORON1065',
        'MIRIM_CORON1140',
        'MIRIM_CORON1550',
        'MIRIM_CORONLYOT'
    ]
    INSTRUMENT_TO_APERNAMES = {
        NIRCAM_A: NIRCAM_A_APERNAMES,
        MIRI: MIRI_APERNAMES
    }
    APERTURE_PA = 1
    V3_PA = 2
    USER_SUPPLIED_COORDS_MSG = '(User-supplied coordinates)'

    def __init__(self):
        self.root = Tk()
        self.root.title("JWST Coronagraph Visibility Tool")

        def close_app():
            self.root.quit()
            self.root.destroy()

        self.root.protocol("WM_DELETE_WINDOW", close_app)
        self.start_year = max(datetime.datetime.today().year, 2018)
        self._build()

    def start(self):
        self.root.lift()
        self.root.call('wm', 'attributes', '.', '-topmost', True)
        self.root.after_idle(self.root.call, 'wm', 'attributes', '.', '-topmost', False)
        self.root.mainloop()

    def error_modal(self, message, title="Error"):
        modal = Toplevel()
        modal.geometry('+400+400')
        modal.title(title)
        frame = ttk.Frame(modal, borderwidth=10)
        frame.grid(column=0, row=0, sticky=(N, S, E, W))
        msg = ttk.Label(frame, text=message)
        msg.grid(column=0, row=0)
        msg.grid_configure(padx=15, pady=15)
        button = ttk.Button(frame, command=modal.destroy, text="OK")
        button.grid(column=0, row=1)
        modal.transient(self.root)
        modal.grab_set()
        self.root.wait_window(modal)

    def show_about(self):
        self.error_modal(
            "The JWST Coronagraph Visibility tool provides approximate\n"
            "pointing restriction information for planning coronagraphic observations.\n\n"
            "For help, contact the helpdesk: help@stsci.edu",
            title="About"
        )

    def _build(self):
        # improve visual feedback for entries in 'disabled' state
        self.style = ttk.Style()
        self.style.map(
            'TEntry',
            background=[('disabled','#d9d9d9'),],
            foreground=[('disabled','#a3a3a3')]
        )
        self.root.minsize(width=1366, height=680)

        # ensure resizing happens:
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        self.main = ttk.Frame(self.root)
        self.main.grid(column=0, row=0, sticky=(N, W, E, S))

        menubar = Menu(self.root)
        appmenu = Menu(menubar, name='apple')
        examples_menu = Menu(menubar)
        self._build_examples_menu(examples_menu)
        menubar.add_cascade(menu=appmenu)
        menubar.add_cascade(menu=examples_menu, label='Examples')
        appmenu.add_command(label='About', command=self.show_about)
        appmenu.add_separator()
        self.root['menu'] = menubar

        # Target, companion, and detector controls
        self.controls_frame = ttk.Frame(self.main, width=240)
        self.controls_frame.grid(column=0, row=0, sticky=(N, W, E, S))
        self._build_controls(self.controls_frame)
        self.controls_frame.grid_propagate(False)

        # plot panel
        self.plot_frame = ttk.Frame(self.main)
        self.plot_frame.grid(column=1, row=0, sticky=(N, W, E, S), columnspan=2)
        self._build_plots(self.plot_frame)

        # massage the gui library a bit
        for child in self.main.winfo_children():
            child.grid_configure(padx=5, pady=5)

        self.main.columnconfigure(1, weight=1)

        self.main.rowconfigure(0, weight=1)

    def _build_controls(self, frame):
        # SIMBAD + RA/Dec
        simbad_frame = ttk.LabelFrame(frame, text="Target Location")
        self._build_simbad_lookup(simbad_frame)
        simbad_frame.grid(column=0, row=0, sticky=(N, W, E, S))

        # Companions
        companion_frame = ttk.LabelFrame(frame, text="Companions")
        self._build_companion_controls(companion_frame)
        companion_frame.grid(column=0, row=1, sticky=(N, W, E, S))
        companion_frame.grid_configure(pady=15)

        # Instrument/Mask selector
        instrument_mask_frame = ttk.Frame(frame)
        self._build_instrument_mask_controls(instrument_mask_frame)
        instrument_mask_frame.grid(column=0, row=2, sticky=(N, W, E, S))
        instrument_mask_frame.grid_configure(pady=15)

        # < > Aperture PA  < > V3 PA
        pa_control_frame = ttk.Frame(frame)
        self.pa_coords = IntVar(value=self.APERTURE_PA)
        aperture_pa_radio = ttk.Radiobutton(
            pa_control_frame,
            text='Aperture PA',
            value=self.APERTURE_PA,
            variable=self.pa_coords
        )
        aperture_pa_radio.grid(column=0, row=0)
        v3_pa_radio = ttk.Radiobutton(
            pa_control_frame,
            text='V3 PA',
            value=self.V3_PA,
            variable=self.pa_coords
        )
        v3_pa_radio.grid(column=1, row=0)
        pa_control_frame.grid(column=0, row=3)

        date_frame = ttk.LabelFrame(frame, text="Date and Sampling")
        self._build_date_controls(date_frame)
        date_frame.grid(column=0, row=4, sticky=(N, W, E, S))
        date_frame.grid_configure(pady=15)

        # Update Plot
        self.update_button = ttk.Button(frame, text="Update Plot", command=self.update_plot)
        self.update_button.grid(column=0, row=5, sticky=(E, W))
        self.progress = ttk.Progressbar(frame, orient='horizontal', mode='indeterminate')
        self.progress.grid(column=0, row=6, sticky=(E, W))
        #
        # examples_frame = ttk.LabelFrame(frame, text="Examples")
        # self._build_examples_frame(examples_frame)
        # examples_frame.grid(column=0, row=7, sticky=(W, E, S), pady=10)
        frame.columnconfigure(0, weight=1)

    def _build_examples_menu(self, menu):
        menu.add_command(label="Single companion, NIRCam 210R spot",
                         command=self._ex_single_companion)
        menu.add_command(label="Three companions, MIRI 4QPM", command=self._ex_three_companions)
        menu.add_command(label="North Ecliptic Pole, NIRCam long wavelength bar", command=self._ex_north_ecliptic)

    def _ex_single_companion(self):
        ra=344.41269
        dec=-29.62224
        pa1=325
        pa2=0
        pa3=0
        separation_as1=10
        separation_as2=0
        separation_as3=0
        apername='NRCA2_MASK210R'

        self.ra_value.set(ra)
        self.dec_value.set(dec)
        visible, pa, sep = self.companions[0]
        visible.set(True)
        pa.set(pa1)
        sep.set(separation_as1)
        visible, pa, sep = self.companions[1]
        visible.set(False)
        pa.set(0)
        sep.set(0)
        visible, pa, sep = self.companions[2]
        visible.set(False)
        pa.set(0)
        sep.set(0)
        self.instrument_value.set(self.NIRCAM_A)
        self.apername_value.set(apername)
        self.simbad_id.set("Example: Single companion")
        self.update_plot()

    def _ex_three_companions(self):
        ra=346.86965
        dec=21.13425
        pa1=45
        pa2=325
        pa3=190
        separation_as1=1.7
        separation_as2=1
        separation_as3=0.65
        apername='MIRIM_CORON1065'

        self.ra_value.set(ra)
        self.dec_value.set(dec)
        visible, pa, sep = self.companions[0]
        visible.set(True)
        pa.set(pa1)
        sep.set(separation_as1)
        visible, pa, sep = self.companions[1]
        visible.set(True)
        pa.set(pa2)
        sep.set(separation_as2)
        visible, pa, sep = self.companions[2]
        visible.set(True)
        pa.set(pa3)
        sep.set(separation_as3)

        self.instrument_value.set(self.MIRI)
        self.apername_value.set(apername)
        self.simbad_id.set("Example: Three companions")
        self.update_plot()

    def _ex_north_ecliptic(self):
        ra=270.0
        dec=66.5
        pa1=0
        pa2=120
        pa3=270
        separation_as1=3
        separation_as2=5
        separation_as3=10
        apername='NRCA5_MASKLWB'

        self.ra_value.set(ra)
        self.dec_value.set(dec)
        visible, pa, sep = self.companions[0]
        visible.set(True)
        pa.set(pa1)
        sep.set(separation_as1)
        visible, pa, sep = self.companions[1]
        visible.set(True)
        pa.set(pa2)
        sep.set(separation_as2)
        visible, pa, sep = self.companions[2]
        visible.set(True)
        pa.set(pa3)
        sep.set(separation_as3)

        self.instrument_value.set(self.NIRCAM_A)
        self.apername_value.set(apername)
        self.simbad_id.set("Example: North Ecliptic Pole")
        self.update_plot()

    def _build_date_controls(self, frame):
        date_label = ttk.Label(frame, text="Start date: October 1,")
        date_label.grid(column=0, row=0, sticky=(N, W))

        self.year_value = StringVar()
        self.year_value.set(self.start_year)
        year_label = ttk.Label(frame, textvariable=self.year_value, width=5)
        year_label.grid(column=1, row=0, sticky=(N, E, W))

        ttk.Label(frame, text="Timesteps per year:").grid(column=0, row=1, sticky=(N, W))
        self.npoints_value = StringVar()
        self.npoints_value.set(DEFAULT_NPOINTS)
        ttk.Entry(frame, textvariable=self.npoints_value, width=5).grid(column=1, row=1, sticky=(N, E, W))

        self.nrolls_value = StringVar()
        ttk.Label(frame, text="Rolls checked:").grid(column=0, row=2, sticky=(N, W))
        self.nrolls_value.set(DEFAULT_NROLLS)
        ttk.Entry(frame, textvariable=self.nrolls_value, width=5).grid(column=1, row=2, sticky=(N, E, W))

    def _build_simbad_lookup(self, frame):
        # SIMBAD lookup
        simbad_label = ttk.Label(frame, text="SIMBAD Target Resolver")
        simbad_label.grid(column=0, row=0, sticky=(N, W), columnspan=4)
        self.simbad_query = StringVar()
        simbad_entry = ttk.Entry(frame, textvariable=self.simbad_query)
        simbad_entry.grid(column=0, row=1, sticky=(N, W, E, S), columnspan=3)
        simbad_entry.bind('<Return>', lambda evt: self.do_simbad_lookup())
        simbad_button = ttk.Button(frame, text="Search", command=self.do_simbad_lookup)
        simbad_button.grid(column=3, row=1)

        # SIMBAD result status
        simbad_id_label = ttk.Label(frame, text="ID:")
        simbad_id_label.grid(column=0, row=2, sticky=(N, W))
        self.simbad_id = StringVar()
        simbad_id_value = ttk.Label(frame, textvariable=self.simbad_id)
        simbad_id_value.grid(column=1, row=2, sticky=(N, W), columnspan=3)

        # RA and Dec
        ra_label = ttk.Label(frame, text="RA:")
        ra_label.grid(column=0, row=3, sticky=(N, W), columnspan=3)
        self.ra_value = StringVar()
        ra_entry = ttk.Entry(frame, textvariable=self.ra_value)
        ra_entry.grid(column=1, row=3, sticky=(N, W, E), columnspan=2)
        ttk.Label(frame, text="º (decimal)").grid(column=3, row=3)

        dec_label = ttk.Label(frame, text="Dec:")
        dec_label.grid(column=0, row=4, sticky=(N, W))
        self.dec_value = StringVar()
        dec_entry = ttk.Entry(frame, textvariable=self.dec_value)
        dec_entry.grid(column=1, row=4, sticky=(N, W, E), columnspan=2)
        ttk.Label(frame, text="º (decimal)").grid(column=3, row=4)

        # Lambda and beta (ecliptic longitude and latitude)
        ecliptic_label = ttk.Label(frame, text="Ecliptic coordinates:")
        ecliptic_label.grid(column=0, row=5, sticky=(N, W), columnspan=4)
        self.ecliptic_value = StringVar()
        ecliptic_display = ttk.Label(frame, textvariable=self.ecliptic_value)
        ecliptic_display.grid(column=0, row=6, sticky=(N, W, E), columnspan=4)

        # Clear the SIMBAD ID when user edits RA or Dec
        def _clear_simbad_id(*_):
            self.simbad_id.set(self.USER_SUPPLIED_COORDS_MSG)

        def _update_ecliptic(*_):
            try:
                 ra, dec = float(self.ra_value.get()), float(self.dec_value.get())
            except ValueError:
                self.ecliptic_value.set('')
                return
            ecliptic_lambda, ecliptic_beta = ad2lb(np.deg2rad(ra), np.deg2rad(dec))
            ecliptic_display_val = '(l, b) = ({:1.4f}º, {:1.4f}º)'.format(
                np.rad2deg(ecliptic_lambda),
                np.rad2deg(ecliptic_beta)
            )
            self.ecliptic_value.set(ecliptic_display_val)

        for var in (self.ra_value, self.dec_value):
            var.trace('w', _clear_simbad_id)
            var.trace('w', _update_ecliptic)

        frame.columnconfigure(1, weight=1)

    def _build_companion_controls(self, frame):
        # (show?) PA deg   Sep arcsec
        ttk.Label(frame, text="PA (º)").grid(column=1, row=0)
        ttk.Label(frame, text="Sep (\")").grid(column=2, row=0)
        self.companions, self.companion_widgets = [], []
        for i in range(1, 4):
            # variables
            visible = BooleanVar(value=False)
            # ensure widgets are updated when `visible` changes:
            def _update_companions(*args):
                self.update_companions()
            visible.trace('w', _update_companions)
            pa = StringVar(value="0.00")
            sep = StringVar(value="0.00")
            self.companions.append((visible, pa, sep))

            # widgets
            check = ttk.Checkbutton(
                frame,
                variable=visible,
                onvalue=True,
                offvalue=False
            )
            check.grid(column=0, row=i)
            pa_entry = ttk.Entry(
                frame,
                textvariable=pa,
                state=DISABLED,
            )
            pa_entry.grid(column=1, row=i, sticky=(E, W))
            sep_entry = ttk.Entry(
                frame,
                textvariable=sep,
                state=DISABLED,
            )
            sep_entry.grid(column=2, row=i, sticky=(E, W))
            self.companion_widgets.append((check, pa_entry, sep_entry))
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(2, weight=1)

    def _build_instrument_mask_controls(self, frame):
        ttk.Label(frame, text="Instrument", anchor=E).grid(column=0, row=0)
        self.instrument_value = StringVar(value=self.NIRCAM_A)
        instrument_combo = ttk.Combobox(
            frame,
            textvariable=self.instrument_value,
            values=self.INSTRUMENTS,
            state='readonly'
        )
        instrument_combo.grid(
            column=0,
            row=1
        )

        initial_apernames = self.INSTRUMENT_TO_APERNAMES[self.NIRCAM_A]
        self.apername_value = StringVar(value=initial_apernames[0])
        ttk.Label(frame, text="Mask", anchor=E).grid(column=0, row=2)
        apername_combo = ttk.Combobox(
            frame,
            textvariable=self.apername_value,
            values=initial_apernames,
            state='readonly'
        )
        apername_combo.grid(
            column=0,
            row=3
        )

        # Hacks to prevent wonky looking text selection within readonly
        # combo boxes
        def _clear_selection_instr(evt):
            instrument_combo.selection_clear()
        instrument_combo.bind('<<ComboboxSelected>>', _clear_selection_instr)

        def _clear_selection_aper(evt):
            apername_combo.selection_clear()
        apername_combo.bind('<<ComboboxSelected>>', _clear_selection_aper)

        # Update apernames based on instrument
        def _update_apernames(*args):
            # throw away args, no useful info there
            values = self.INSTRUMENT_TO_APERNAMES[self.instrument_value.get()]
            apername_combo['values'] = values
            self.apername_value.set(values[0])
        self.instrument_value.trace('w', _update_apernames)

    def _build_plots(self, frame):
        self.figure = Figure(figsize=(8, 8), dpi=72)

        # initialized when the plot is updated:
        self._pick_event_handler_id = None
        self._plot_overlay_elements = []
        self._mask_artists = []

        obs_axes = (0.1, 0.3, 0.35, 0.6)  # (left, bottom, width, height)
        self.observability_ax = self.figure.add_axes(obs_axes)

        detector_axes = (0.55, 0.3, 0.4, 0.6)
        self.detector_ax = self.figure.add_axes(detector_axes)
        self.detector_ax.set_aspect('equal', anchor='SE')

        # companion legend markers
        self.companion_legend_markers = []
        self.companion_legend_labels = []
        self.companion_info = []
        v_pos = 0.2
        self.observable_pa = matplotlib.text.Text(x=0.55, y=v_pos, text="PA:", transform=self.figure.transFigure, figure=self.figure)
        self.figure.texts.append(self.observable_pa)

        line_height = 0.04
        for i, color in enumerate((RED_GGPLOT, BLUE_GGPLOT, PURPLE_GGPLOT)):
            v_pos -= line_height
            marker = matplotlib.patches.Rectangle((0.55, v_pos), width=0.01, height=0.015, facecolor=color, transform=self.figure.transFigure, figure=self.figure)
            self.figure.patches.append(marker)
            self.companion_legend_markers.append(marker)

            label = matplotlib.text.Text(x=0.57, y=v_pos, text="Companion {}".format(i + 1), transform=self.figure.transFigure, figure=self.figure)
            self.figure.texts.append(label)
            self.companion_legend_labels.append(label)

            info = matplotlib.text.Text(x=0.57, y=v_pos - line_height / 2, text="# arcsec @ # deg", transform=self.figure.transFigure, figure=self.figure)
            self.figure.texts.append(info)
            self.companion_info.append(info)

        self._canvas = FigureCanvasTkAgg(self.figure, master=frame)
        self._canvas.show()
        self._canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self._toolbar = NavigationToolbar2TkAgg(self._canvas, frame)
        self._toolbar.update()
        self._canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        def on_key_event(event):
            key_press_handler(event, self._canvas, self._toolbar)

        self._canvas.mpl_connect('key_press_event', on_key_event)

    def do_simbad_lookup(self):
        search_string = self.simbad_query.get()
        if not len(search_string.strip()) > 0:
            self.error_modal("Search query for SIMBAD must not be empty")
            return

        with _busy_cursor(self.root):
            try:
                result = query_simbad(search_string.strip())
            except requests.exceptions.Timeout:
                self.error_modal("Cannot reach SIMBAD! Check your network connection, or see if perhaps SIMBAD is down...")
                return

            if result is None:
                self.error_modal("No object found for this identifier! Try a different query, or supply RA and Dec manually.")
                return

            self.ra_value.set(str(result.ra))
            self.dec_value.set(str(result.dec))
            self.simbad_id.set(result.id)


    def update_companions(self):
        # handle disabling / enabling entries
        for comp, widg in zip(self.companions, self.companion_widgets):
            visible, pa, sep = comp
            check, pa_entry, sep_entry = widg
            if visible.get():
                pa_entry.config(state="normal")
                sep_entry.config(state="normal")
            else:
                pa_entry.config(state="disabled")
                sep_entry.config(state="disabled")

    def update_plot(self):
        try:
            ra = float(self.ra_value.get())
            dec = float(self.dec_value.get())
        except ValueError:
            self.error_modal("RA and Declination must be given in decimal degrees")
            return
        if ra > 360 or ra < 0:
            self.error_modal("RA must be between 0 and 360 degrees")
            return
        if dec > 90 or dec < -90:
            self.error_modal("Declination must be between -90 and 90 degrees")
            return
        try:
            npoints = int(self.npoints_value.get())
            nrolls = int(self.nrolls_value.get())
        except ValueError:
            self.error_modal("Number of points and roll angle sampling must be integers")
            return

        try:
            start_year = int(self.year_value.get())
            if start_year < 2000:
                raise ValueError("sun_ecliptic_longitude works for years after 2000 only")
            start_date = datetime.datetime(start_year, 10, 1)
        except ValueError:
            self.error_modal("Supply a four-digit year after 2000")

        # ugly loop unroll for the 3 companions
        shown, pa, sep = self.companions[0]
        if shown.get():
            pa1 = float(pa.get())
            separation_as1 = float(sep.get())
        else:
            pa1 = 0.0
            separation_as1 = 0.0
        shown, pa, sep = self.companions[1]
        if shown.get():
            pa2 = float(pa.get())
            separation_as2 = float(sep.get())
        else:
            pa2 = 0.0
            separation_as2 = 0.0
        shown, pa, sep = self.companions[2]
        if shown.get():
            pa3 = float(pa.get())
            separation_as3 = float(sep.get())
        else:
            pa3 = 0.0
            separation_as3 = 0.0

        instrument = self.instrument_value.get()
        if instrument == self.NIRCAM_A:
            instrname = 'NIRCam'
        elif instrument == self.MIRI:
            instrname = 'MIRI'
        else:
            raise Exception("Unsupported instrument!")
        apername = self.apername_value.get()

        # busy cursor start
        self.update_button.config(state='disabled')
        self.progress.start()
        with _busy_cursor(self.root):
            aper = get_aperture(instrname, apername)
            self.result = VisibilityCalculation(
                ra,
                dec,
                [
                    {'pa': pa1, 'separation': separation_as1},
                    {'pa': pa2, 'separation': separation_as2},
                    {'pa': pa3, 'separation': separation_as3},
                ],
                aper,
                start_date,
                npoints,
                nrolls
            )
            self.result.calculate()

            self._clear_plot_overlay()
            self._update_observability()
            if self._pick_event_handler_id is None:
                self._pick_event_handler_id = self.figure.canvas.mpl_connect('pick_event', self._on_pick)
            self._update_detector()
            self._canvas.show()

        self.progress.stop()
        self.update_button.config(state='normal')

    def _update_observability(self):
        days = self.result.days
        elongation_rad = self.result.elongation_rad
        roll_rad = self.result.roll_rad
        observable = self.result.observable
        ax = self.observability_ax
        ax.clear()
        if self.simbad_id.get() == self.USER_SUPPLIED_COORDS_MSG:
            ax.set_title('Observability of\nRA: {:3.5f} deg Dec: {:+3.5f} deg'.format(self.result.ra, self.result.dec))
        else:
            ax.set_title('Observability of {}'.format(self.simbad_id.get()))

        (elongation_line,) = ax.plot(days, np.rad2deg(elongation_rad[0]), color='black', label='Solar elongation')  # same for all 20 roll angles?? pick first

        collapsed_mask = np.any(observable, axis=0)
        observable_series = ax.scatter(
            days[collapsed_mask],
            np.rad2deg(elongation_rad[0])[collapsed_mask],
            color='none',
            marker='o',
            edgecolor=RED_GGPLOT,
            s=40,
            label='Observable Elongations'
        )

        self._last_plotted_pa = self.pa_coords.get()
        if self._last_plotted_pa == self.APERTURE_PA:
            # Aperture PA
            pa_label = 'Aperture PA'
            pa_color = BLUE_GGPLOT

            # n.b. sciyangle is incorrect in PRDSOCDEV 012
            theta = np.rad2deg(np.arctan2(self.result.n_x, self.result.n_y))
        else:
            # v3 PA
            pa_label = 'V3 PA'
            pa_color = PURPLE_GGPLOT
            theta = np.rad2deg(roll_rad)

        theta %= 360

        mask = observable != 0
        # there might be a better way to get a 'days' the right shape
        days_for_all_rolls = np.repeat(days[np.newaxis,:], self.result.nrolls, axis=0)
        days_for_all_rolls[self.result.observable == 0] = np.nan
        theta[self.result.observable == 0] = np.nan
        # TODO there should be a more elegant way to hold on to the actual plotted arrays
        # for later interactivity
        self._days_for_all_rolls, self._theta = days_for_all_rolls, theta

        self._pa_series = ax.scatter(days_for_all_rolls, theta, color=pa_color, label=pa_label, picker=True)

        ax.set_xlim(0, 366)
        ax.set_xlabel('Days since Oct 1 {}'.format(self.result.start_date.year))
        legend = ax.legend(
            (elongation_line, observable_series, self._pa_series),
            ('Solar elongation', 'Observable elongations', pa_label),
            bbox_to_anchor=(0.1, 0.1, 0.36, .102),
            bbox_transform=self.figure.transFigure,
            mode="expand", borderaxespad=0.,
            framealpha=0.0,
        )

        ax.set_ylim(0, 400)
        ax.set_ylabel('Degrees')

    def work_backwards(self, x_array, y_array, xdata, ydata):
        dist = (x_array - xdata)**2 + (y_array - ydata)**2
        dist[self.result.observable == 0] = np.nan
        y, x = np.unravel_index(np.nanargmin(dist), dist.shape)
        return y, x

    def _on_pick(self, event):
        self._clear_plot_overlay()
        if event.artist.axes == self.detector_ax:
            self._on_detector_pick(event)
        elif event.artist.axes == self.observability_ax:
            self._on_observability_pick(event)

    def _on_observability_pick(self, event):
        yidx, xidx = self.work_backwards(self._days_for_all_rolls, self._theta, event.mouseevent.xdata, event.mouseevent.ydata)
        self._add_plot_overlay(yidx, xidx)

    def _on_detector_pick(self, event):
        companions = (
            (self.c1_plot_group, (self.result.c1_x, self.result.c1_y)),
            (self.c2_plot_group, (self.result.c2_x, self.result.c2_y)),
            (self.c3_plot_group, (self.result.c3_x, self.result.c3_y)),
        )
        for idx, (artist, (xarr, yarr)) in enumerate(companions):
            if self.result.companions[idx]['separation'] == 0.0:
                continue
            if artist == event.artist:
                yidx, xidx = self.work_backwards(xarr, yarr, event.mouseevent.xdata, event.mouseevent.ydata)
                self._add_plot_overlay(yidx, xidx)
                return

    def _clear_plot_overlay(self):
        while len(self._plot_overlay_elements):
            elem = self._plot_overlay_elements.pop()
            elem.remove()
        for text in self.companion_info:
            text.set_text('')

    def _add_plot_overlay(self, yidx, xidx):
        obs_highlight = self.observability_ax.scatter(self._days_for_all_rolls[yidx, xidx], self._theta[yidx, xidx], color='white', edgecolor='black', s=100)
        self._plot_overlay_elements.append(obs_highlight)
        obs_vline = self.observability_ax.axvline(self._days_for_all_rolls[yidx, xidx], color=BLUE_GGPLOT)
        self._plot_overlay_elements.append(obs_vline)
        obs_hline = self.observability_ax.axhline(self._theta[yidx, xidx], color=BLUE_GGPLOT)
        self._plot_overlay_elements.append(obs_hline)
        if self._last_plotted_pa == self.APERTURE_PA:
            pa_label = 'Aperture PA'
        else:
            pa_label = 'V3 PA'

        self.observable_pa.set_text("{pa_label} = {pa:.2f} deg".format(
            pa_label=pa_label,
            pa=self._theta[yidx, xidx],
        ))

        for idx, companion in enumerate(self.result.companions):
            if companion['separation'] == 0:
                continue
            c_x = getattr(self.result, 'c{}_x'.format(idx + 1))
            c_y = getattr(self.result, 'c{}_y'.format(idx + 1))
            x, y = c_x[yidx, xidx], c_y[yidx, xidx]
            highlight = self.detector_ax.scatter(x, y, color='white', edgecolor='black', s=100)
            self.companion_info[idx].set_text('{dist:.2f} arcsec @ {angle:.2f} deg'.format(
                dist=companion['separation'],
                angle=np.rad2deg(np.arctan2(-x, y))
            ))
            self._plot_overlay_elements.append(highlight)

        separations = [self.result.companions[i]['separation'] for i in range(3)]
        if np.max(separations) > 0:
            scale_factor = 1.1 * np.max(separations)
        else:
            arcsec_per_pixel = np.average([self.result.aperture.XSciScale, self.result.aperture.YSciScale])
            scale_factor = self.result.aperture.XSciSize * arcsec_per_pixel / 2.0

        n_x_temp = self.result.n_x[yidx, xidx] / np.sqrt(self.result.n_x[yidx, xidx]**2 + self.result.n_y[yidx, xidx]**2)
        n_y_temp = self.result.n_y[yidx, xidx] / np.sqrt(self.result.n_x[yidx, xidx]**2 + self.result.n_y[yidx, xidx]**2)

        u = np.array([0, 1])
        north_line, = self.detector_ax.plot(scale_factor * n_x_temp * u, scale_factor * n_y_temp * u, color=RED_GGPLOT)
        self._plot_overlay_elements.append(north_line)

        north_label = self.detector_ax.text(scale_factor / 2 * n_x_temp, scale_factor / 2 * n_y_temp, "N")
        self._plot_overlay_elements.append(north_label)

        e_x_temp = self.result.e_x[yidx, xidx] / np.sqrt(self.result.e_x[yidx, xidx]**2 + self.result.e_y[yidx, xidx]**2)
        e_y_temp = self.result.e_y[yidx, xidx] / np.sqrt(self.result.e_x[yidx, xidx]**2 + self.result.e_y[yidx, xidx]**2)

        east_line, = self.detector_ax.plot(scale_factor * e_x_temp * u, scale_factor * e_y_temp * u, color=YELLOW_GGPLOT)
        self._plot_overlay_elements.append(east_line)
        east_label = self.detector_ax.text(scale_factor / 2 * e_x_temp, scale_factor / 2 * e_y_temp, "E")
        self._plot_overlay_elements.append(east_label)

        self._canvas.show()

    def _update_detector(self):
        ax = self.detector_ax
        ax.clear()
        aperture = self.result.aperture
        arcsec_per_pixel = np.average([aperture.XSciScale, aperture.YSciScale])
        ax.set_title('{name}\n({x_size:.0f} x {y_size:.0f} pixels, {scale:1.4f} arcsec/pixel)'.format(
            name=aperture.AperName,
            x_size=aperture.XSciSize,
            y_size=aperture.YSciSize,
            scale=arcsec_per_pixel,
        ))
        self._mask_artists = []
        ax.set_aspect('equal')

        aper_corners_x, aper_corners_y = aperture.corners(frame='Idl')
        verts = np.concatenate([aper_corners_x[:,np.newaxis], aper_corners_y[:,np.newaxis]], axis=1)
        patch = patches.Polygon(verts, facecolor='none', edgecolor='red', alpha=0.5, linestyle='--', linewidth=3)
        ax.add_artist(patch)

        self.c1_plot_group = ax.scatter(self.result.c1_x, self.result.c1_y, picker=True, color=RED_GGPLOT)
        self.c2_plot_group = ax.scatter(self.result.c2_x, self.result.c2_y, picker=True, color=BLUE_GGPLOT)
        self.c3_plot_group = ax.scatter(self.result.c3_x, self.result.c3_y, picker=True, color=PURPLE_GGPLOT)

        ax.set_xlim(np.min(aper_corners_x) - 5, np.max(aper_corners_x) + 5)
        ax.set_ylim(np.min(aper_corners_y) - 5, np.max(aper_corners_y) + 5)
        ax.set_xlabel('x (arcsec, ideal frame)')
        ax.set_ylabel('y (arcsec, ideal frame)')

        self._overlay_mask()

    def _overlay_mask(self):
        while self._mask_artists:
            artist = self._mask_artists.pop()
            artist.remove()

        aperture = self.result.aperture
        aperture_name = aperture.AperName
        arcsec_per_pixel = np.average([aperture.XSciScale, aperture.YSciScale])
        x_sci_size, y_sci_size = aperture.XSciSize, aperture.YSciSize
        mask_patches = []

        if 'NRC' in aperture_name:
            for quad_verts in NIRCAM_CORON_BAD_AREAS:
                v2, v3 = quad_verts[:,0], quad_verts[:,1]
                xidl, yidl = aperture.Tel2Idl(v2, v3)
                idl_verts = np.concatenate([xidl[:,np.newaxis], yidl[:,np.newaxis]], axis=1)
                patch = patches.Polygon(idl_verts, facecolor='red', edgecolor='none', alpha=0.5)
                mask_patches.append(patch)
            if aperture_name[-1] == 'R':
                if '210R' in aperture_name:
                    radius_arcsec = 0.40
                elif '335R' in aperture_name:
                    radius_arcsec = 0.64
                elif '430R' in aperture_name:
                    radius_arcsec = 0.82
                else:
                    raise RuntimeError("Invalid mask!")
                # make a circle
                mask_patches.append(patches.Circle((0, 0), radius=radius_arcsec, alpha=0.5))
            else:
                if 'LWB' in aperture_name:
                    thin_extent_arcsec = 0.58 * (2 / 4)
                    thick_extent_arcsec = 0.58 * (6 / 4)
                elif 'SWB' in aperture_name:
                    thin_extent_arcsec = 0.27 * (2 / 4)
                    thick_extent_arcsec = 0.27 * (6 / 4)
                else:
                    raise RuntimeError("Invalid mask!")

                x_verts = x_sci_size / 2 * np.array([-1, 1, 1, -1])
                y_verts = np.array([
                    thin_extent_arcsec / arcsec_per_pixel,
                    thick_extent_arcsec / arcsec_per_pixel,
                    -thick_extent_arcsec / arcsec_per_pixel,
                    -thin_extent_arcsec / arcsec_per_pixel
                ])
                x_idl_verts, y_idl_verts = aperture.Sci2Idl(x_verts + aperture.XSciRef, y_verts + aperture.YSciRef)
                verts = np.concatenate([x_idl_verts[:,np.newaxis], y_idl_verts[:,np.newaxis]], axis=1)
                patch = patches.Polygon(verts, alpha=0.5)
                mask_patches.append(patch)
                # self._mask_artists = self.detector_ax.add_artist(patch)
        elif 'MIRI' in aperture_name:
            y_angle = np.deg2rad(aperture.V3IdlYAngle)
            corners_x, corners_y = aperture.corners(frame='Idl')
            min_x, min_y = np.min(corners_x), np.min(corners_y)
            max_x, max_y = np.max(corners_x), np.max(corners_y)

            if 'LYOT' in aperture_name:
                # David Law, personal communication, May 2016:
                # The clear-aperture area for the Lyot is 272x272 pixels

                width_arcsec = 0.72
                x_verts = width_arcsec * np.array([-1, -1, 1, 1]) / 2
                y_verts = np.array([min_y, max_y, max_y, min_y])
                x_verts = np.cos(y_angle) * x_verts + np.sin(y_angle) * y_verts
                y_verts = -np.sin(y_angle) * x_verts + np.cos(y_angle) * y_verts

                verts = np.concatenate([x_verts[:,np.newaxis], y_verts[:,np.newaxis]], axis=1)
                rectangular_part = patches.Polygon(verts)
                # already in Idl coords
                radius_arcsec = 2.16
                circular_part = patches.Circle((0, 0), radius=radius_arcsec)
                mask_collection = PatchCollection([rectangular_part, circular_part], alpha=0.5)
                mask_patches.append(mask_collection)
            elif '1065' in aperture_name or '1140' in aperture_name or '1550' in aperture_name:
                width_arcsec = 0.33
                # David Law, personal communication, May 2016:
                # The clear-aperture area for the 4QPM is 216x216 pixels
                x_verts = np.array([
                    min_x,
                    -width_arcsec,
                    -width_arcsec,
                    width_arcsec,
                    width_arcsec,
                    max_x,
                    max_x,
                    width_arcsec,
                    width_arcsec,
                    -width_arcsec,
                    -width_arcsec,
                    min_x
                ])
                y_verts = np.array([
                    width_arcsec,
                    width_arcsec,
                    max_y,
                    max_y,
                    width_arcsec,
                    width_arcsec,
                    -width_arcsec,
                    -width_arcsec,
                    min_y,
                    min_y,
                    -width_arcsec,
                    -width_arcsec
                ])
                x_verts = np.cos(y_angle) * x_verts + np.sin(y_angle) * y_verts
                y_verts = -np.sin(y_angle) * x_verts + np.cos(y_angle) * y_verts

                verts = np.concatenate([x_verts[:, np.newaxis], y_verts[:, np.newaxis]], axis=1)
                mask_patches.append(patches.Polygon(verts, alpha=0.5))
            else:
                raise RuntimeError("Invalid mask!")

        for patch in mask_patches:
            self._mask_artists.append(self.detector_ax.add_artist(patch))


def run():
    app = VisibilityCalculator()
    app.start()
