#!/usr/bin/env python
# vim: set fileencoding=utf8 :
from __future__ import print_function, division

try:
    from tkinter import *
    from tkinter import ttk
except ImportError:
    from Tkinter import *
    import ttk
import os.path

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

import astropy.coordinates as coord
import astropy.units as u
from astroquery.simbad import Simbad

from .jwxml import SIAF
from .skyvec2ins import skyvec2ins

from pprint import pprint

RED_GGPLOT = '#E24A33'
BLUE_GGPLOT = '#348ABD'
PURPLE_GGPLOT = '#988ED5'
GRAY_GGPLOT = '#777777'
YELLOW_GGPLOT = '#FBC15E'
GREEN_GGPLOT = '#8EBA42'
PINK_GGPLOT = '#FFB5B8'

class VisibilityCalculation(object):
    def __init__(self, ra, dec, companions, aperture, npoints, nrolls):
        self.ra = ra
        self.dec = dec
        self.companions = companions
        self.aperture = aperture
        self.npoints = npoints
        self.nrolls = nrolls

        # Outputs
        self.x = None
        self.observable = None
        self.elongation_rad = None
        self.roll_rad = None

        self.scisize = None
        self.sciscale = None
        self.sciyangle = None

        self.s_x = None
        self.s_y = None
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
            self.x,
            self.observable,
            self.elongation_rad,
            self.roll_rad,
            self.s_x, self.s_y,
            self.c1_x, self.c1_y,
            self.c2_x, self.c2_y,
            self.c3_x, self.c3_y,
            self.scisize, self.sciscale, self.sciyangle,
            self.n_x, self.n_y, self.e_x, self.e_y
        ) = skyvec2ins(
            ra=self.ra,
            dec=self.dec,
            pa1=self.companions[0]['pa'],
            pa2=self.companions[1]['pa'],
            pa3=self.companions[2]['pa'],
            separation_as1=self.companions[0]['separation'],
            separation_as2=self.companions[1]['separation'],
            separation_as3=self.companions[2]['separation'],
            # TODO pass Aperture object instead of instrument + apername
            instrname=self.aperture.instrument,
            apername=self.aperture.AperName,
            npoints=self.npoints,
            nrolls=self.nrolls
        )

        self.days = self.x
        # mask non-observable (roll, elongation) pairs from output data
        mask = self.observable == 0
        self.s_x[mask] = np.nan
        self.s_y[mask] = np.nan

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
    INSTRUMENTS = [NIRCAM_A, NIRCAM_B, MIRI]
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
        'MIRIM_MASK1065',
        'MIRIM_MASK1140',
        'MIRIM_MASK1550',
        'MIRIM_MASKLYOT'
    ]
    INSTRUMENT_TO_APERNAMES = {
        NIRCAM_A: NIRCAM_A_APERNAMES,
        NIRCAM_B: NIRCAM_B_APERNAMES,
        MIRI: MIRI_APERNAMES
    }
    DETECTOR_PA = 1
    V3_PA = 2
    USER_SUPPLIED_COORDS_MSG = '(User-supplied coordinates)'
    def __init__(self):
        self.root = Tk()
        self.root.title("JWST Visibility Calculator")

        def close_app():
            self.root.quit()
            self.root.destroy()

        self.root.protocol("WM_DELETE_WINDOW", close_app)
        self._build()

    def start(self):
        self.root.lift()
        self.root.call('wm', 'attributes', '.', '-topmost', True)
        self.root.after_idle(self.root.call, 'wm', 'attributes', '.', '-topmost', False)
        self.root.mainloop()

    def error_modal(self, message):
        modal = Toplevel()
        modal.geometry('+400+400')
        modal.title("Error")
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

    def _build(self):
        # improve visual feedback for entries in 'disabled' state
        self.style = ttk.Style()
        self.style.map(
            'TEntry',
            background=[('disabled','#d9d9d9'),],
            foreground=[('disabled','#a3a3a3')]
        )
        self.root.minsize(width=900, height=500)
        # ensure resizing happens:
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        self.main = ttk.Frame(self.root)
        self.main.grid(column=0, row=0, sticky=(N, W, E, S))

        # Target, companion, and detector controls
        self.controls_frame = ttk.Frame(self.main)
        self.controls_frame.grid(column=0, row=0, sticky=(N, W, E, S))
        self._build_controls(self.controls_frame)

        # plot panel
        self.plot_frame = ttk.Frame(self.main)
        self.plot_frame.grid(column=1, row=0, sticky=(N, W, E, S), columnspan=2)
        self._build_plots(self.plot_frame)

        # massage the gui library a bit
        for child in self.main.winfo_children():
            child.grid_configure(padx=5, pady=5)

        self.main.columnconfigure(1, weight=1)
        self.main.columnconfigure(2, weight=1)
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
        instrument_mask_frame = ttk.LabelFrame(frame, text="Instrument/Mask Selector")
        self._build_instrument_mask_controls(instrument_mask_frame)
        instrument_mask_frame.grid(column=0, row=2, sticky=(N, W, E, S))
        instrument_mask_frame.grid_configure(pady=15)

        # < > Detector PA  < > V3 PA
        pa_control_frame = ttk.Frame(frame)
        self.pa_coords = IntVar(value=self.DETECTOR_PA)
        detector_pa_radio = ttk.Radiobutton(
            pa_control_frame,
            text='Detector PA',
            value=self.DETECTOR_PA,
            variable=self.pa_coords
        )
        detector_pa_radio.grid(column=0, row=0)
        v3_pa_radio = ttk.Radiobutton(
            pa_control_frame,
            text='V3 PA',
            value=self.V3_PA,
            variable=self.pa_coords
        )
        v3_pa_radio.grid(column=1, row=0)
        pa_control_frame.grid(column=0, row=3)
        # Update Plot
        update_button = ttk.Button(frame, text="Update Plot", command=self.update_plot)
        update_button.grid(column=0, row=4, sticky=(E, W))

        # TODO make this data driven or get rid of it
        examples_frame = ttk.LabelFrame(frame, text="Examples")
        self._build_examples_frame(examples_frame)
        examples_frame.grid(column=0, row=5, sticky=(N, W, E, S))
        frame.columnconfigure(0, weight=1)

    def _build_examples_frame(self, frame):
        ttk.Button(
            frame,
            text="Fomalhaut on NIRCam MASK210R",
            command=self._ex_fomalhaut
        ).grid(column=0, row=0, sticky=(E, W))
        ttk.Button(
            frame,
            text="1RXS J160929.1-210524 on NIRCam B MASKSWB",
            command=self._ex_1RXSJ160929p1_210524
        ).grid(column=0, row=1, sticky=(E, W))
        ttk.Button(
            frame,
            text="HR 8799 on MIRI MASK1065",
            command=self._ex_HR8799
        ).grid(column=0, row=2, sticky=(E, W))

    def _ex_fomalhaut(self):
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
        self.update_plot()

    def _ex_1RXSJ160929p1_210524(self):
        ra=242.37628
        dec=-21.08304
        pa1=20
        pa2=0
        pa3=0
        separation_as1=3
        separation_as2=0
        separation_as3=0
        apername='NRCB3_MASKSWB'

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
        self.instrument_value.set(self.NIRCAM_B)
        self.apername_value.set(apername)
        self.update_plot()

    def _ex_HR8799(self):
        ra=346.86965
        dec=21.13425
        pa1=45
        pa2=325
        pa3=190
        separation_as1=1.7
        separation_as2=1
        separation_as3=0.65
        apername='MIRIM_MASK1065'

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
        self.update_plot()

    def _build_simbad_lookup(self, frame):
        # SIMBAD lookup
        simbad_label = ttk.Label(frame, text="SIMBAD Target Resolver")
        simbad_label.grid(column=0, row=0, sticky=(N, W), columnspan=4)
        self.simbad_query = StringVar()
        simbad_entry = ttk.Entry(frame, textvariable=self.simbad_query)
        simbad_entry.grid(column=0, row=1, sticky=(N, W, E, S), columnspan=3)
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
        ra_label.grid(column=0, row=3, sticky=(N, W))
        self.ra_value = StringVar()
        ra_entry = ttk.Entry(frame, textvariable=self.ra_value)
        ra_entry.grid(column=1, row=3, sticky=(N, W, E), columnspan=2)
        ttk.Label(frame, text="degrees (decimal)").grid(column=3, row=3)

        dec_label = ttk.Label(frame, text="Dec:")
        dec_label.grid(column=0, row=4, sticky=(N, W))
        self.dec_value = StringVar()
        dec_entry = ttk.Entry(frame, textvariable=self.dec_value)
        dec_entry.grid(column=1, row=4, sticky=(N, W, E), columnspan=2)
        ttk.Label(frame, text="degrees (decimal)").grid(column=3, row=4)

        # Clear the SIMBAD ID when user edits RA or Dec
        def _user_edited_coords(*args):
            self.simbad_id.set(self.USER_SUPPLIED_COORDS_MSG)
        self.ra_value.trace('w', _user_edited_coords)

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
        ttk.Label(frame, text="Instrument").grid(column=0, row=0)
        self.instrument_value = StringVar(value=self.NIRCAM_A)
        instrument_combo = ttk.Combobox(
            frame,
            textvariable=self.instrument_value,
            values=self.INSTRUMENTS,
            state='readonly'
        )
        instrument_combo.grid(
            column=1,
            row=0
        )

        initial_apernames = self.INSTRUMENT_TO_APERNAMES[self.NIRCAM_A]
        self.apername_value = StringVar(value=initial_apernames[0])
        ttk.Label(frame, text="Mask").grid(column=0, row=1)
        apername_combo = ttk.Combobox(
            frame,
            textvariable=self.apername_value,
            values=initial_apernames,
            state='readonly'
        )
        apername_combo.grid(
            column=1,
            row=1
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

    def _build_observability(self, frame):
        self.observability_figure = Figure(figsize=(8, 8), dpi=72)
        self.observability_ax = self.observability_figure.add_subplot(1, 1, 1)

        self._observability_canvas = FigureCanvasTkAgg(self.observability_figure, master=frame)
        self._observability_canvas.show()
        self._observability_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self._observability_toolbar = NavigationToolbar2TkAgg(self._observability_canvas, frame)
        self._observability_toolbar.update()
        self._observability_canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        def on_key_event(event):
            key_press_handler(event, self._observability_canvas, self._observability_toolbar)

        self._observability_canvas.mpl_connect('key_press_event', on_key_event)

    def _build_detector(self, frame):
        self.detector_figure = Figure(figsize=(8, 8), dpi=72)
        self.detector_ax = self.detector_figure.add_subplot(1, 1, 1)
        self.detector_ax.set_aspect('equal')

        self._detector_canvas = FigureCanvasTkAgg(self.detector_figure, master=frame)
        self._detector_canvas.show()
        self._detector_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self._detector_toolbar = NavigationToolbar2TkAgg(self._detector_canvas, frame)
        self._detector_toolbar.update()
        self._detector_canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        def on_key_event(event):
            # print('you pressed %s' % event.key)
            key_press_handler(event, self._detector_canvas, self._detector_toolbar)

        self._detector_canvas.mpl_connect('key_press_event', on_key_event)

    def _build_plots(self, frame):
        self.figure = Figure(figsize=(8, 8), dpi=72)

        # initialized when the plot is updated:
        self._pick_event_handler_id = None
        self._plot_overlay_elements = []
        self._mask_artist = None

        self.observability_ax = self.figure.add_subplot(1, 2, 1)
        self.detector_ax = self.figure.add_subplot(1, 2, 2)
        self.detector_ax.set_aspect('equal')

        self._canvas = FigureCanvasTkAgg(self.figure, master=frame)
        self._canvas.show()
        self._canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self._toolbar = NavigationToolbar2TkAgg(self._canvas, frame)
        self._toolbar.update()
        self._canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        def on_key_event(event):
            # print('you pressed %s' % event.key)
            key_press_handler(event, self._canvas, self._toolbar)

        self._canvas.mpl_connect('key_press_event', on_key_event)

    def do_simbad_lookup(self):
        search_string = self.simbad_query.get()
        if not len(search_string.strip()) > 0:
            self.error_modal("Search query for SIMBAD must not be empty")
            return
        self.root.config(cursor='wait')
        self.root.update()
        result_table = Simbad.query_object(search_string)
        if len(result_table) > 1:
            self.error_modal("More than one object found for this identifier! Try a more specific query, or supply RA and Dec manually.")
            return
        if len(result_table) == 0:
            self.error_modal("No objects found for this identifier! Try a different query, or supply RA and Dec manually.")
            return
        ra_str, dec_str = result_table[0]['RA'], result_table[0]['DEC']
        ra_deg = coord.Angle(ra_str, unit=u.hour).degree
        dec_deg = coord.Angle(dec_str, unit=u.degree).degree
        self.ra_value.set(str(ra_deg))
        self.dec_value.set(str(dec_deg))
        self.simbad_id.set(result_table[0]['MAIN_ID'])
        self.root.config(cursor='')
        self.root.update()

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
        ra = float(self.ra_value.get())
        dec = float(self.dec_value.get())

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
        if instrument in (self.NIRCAM_A, self.NIRCAM_B):
            instrname = 'NIRCam'
        else:
            instrname = 'MIRI'
        apername = self.apername_value.get()

        # TODO we want to eventually pass jwxml.Aperture instances into the calculation
        siaf_path = os.path.join(os.path.dirname(__file__), 'data', '{}_SIAF.xml'.format(instrname))
        assert os.path.exists(siaf_path), 'no SIAF for {}'.format(instrname)
        siaf = SIAF(instr=instrname, filename=siaf_path)
        aper = siaf[apername]

        npoints = 360
        nrolls = 20

        self.result = VisibilityCalculation(
            ra,
            dec,
            [
                {'pa': pa1, 'separation': separation_as1},
                {'pa': pa2, 'separation': separation_as2},
                {'pa': pa3, 'separation': separation_as3},
            ],
            aper,
            npoints,
            nrolls
        )
        # busy cursor start
        self.root.config(cursor='wait')
        self.root.update()
        self.result.calculate()

        if self.simbad_id.get() == self.USER_SUPPLIED_COORDS_MSG:
            self.figure.suptitle('Target at RA: {:3.5f} deg Dec: {:+3.5f} deg'.format(ra, dec))
        else:
            self.figure.suptitle('{} at RA: {:3.5f} deg Dec: {:+3.5f} deg'.format(self.simbad_id.get(), ra, dec))

        self._clear_plot_overlay()
        self._update_observability()
        if self._pick_event_handler_id is None:
            self._pick_event_handler_id = self.figure.canvas.mpl_connect('pick_event', self._on_pick)
        self._update_detector()
        self._canvas.show()
        # busy cursor end
        self.root.config(cursor='')
        self.root.update()

    def _update_observability(self):
        days = self.result.days
        elongation_rad = self.result.elongation_rad
        roll_rad = self.result.roll_rad
        observable = self.result.observable
        sciyangle = self.result.sciyangle

        ax = self.observability_ax
        ax.clear()
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

        if self.pa_coords.get() == self.DETECTOR_PA:
            # detector PA
            pa_label = 'Detector PA'
            pa_color = BLUE_GGPLOT
            theta = np.rad2deg(roll_rad + sciyangle)
        else:
            # v3 PA
            pa_label = 'V3 PA'
            pa_color = PURPLE_GGPLOT
            theta = np.rad2deg(roll_rad)

        theta %= 360

        mask = observable != 0
        # there might be a better way to get a 'days' the right shape
        days_for_all_rolls = np.repeat(days[np.newaxis,:], 20, axis=0)
        days_for_all_rolls[self.result.observable == 0] = np.nan
        theta[self.result.observable == 0] = np.nan
        # TODO there should be a more elegant way to hold on to the actual plotted arrays
        # for later interactivity
        self._days_for_all_rolls, self._theta = days_for_all_rolls, theta

        self._pa_series = ax.scatter(days_for_all_rolls, theta, color=pa_color, label=pa_label, picker=True)

        ax.set_xlim(0, 366)
        ax.set_xlabel('Days since Oct 1 2018')
        ax.legend(
            (elongation_line, observable_series, self._pa_series),
            ('Solar elongation', 'Observable elongations', pa_label)
        )

        ax.set_ylim(0, 360)
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
        for artist, (xarr, yarr) in companions:
            if artist == event.artist:
                yidx, xidx = self.work_backwards(xarr, yarr, event.mouseevent.xdata, event.mouseevent.ydata)
                self._add_plot_overlay(yidx, xidx)
                return

    def _clear_plot_overlay(self):
        while len(self._plot_overlay_elements):
            elem = self._plot_overlay_elements.pop()
            elem.remove()

    def _add_plot_overlay(self, yidx, xidx):
        obs_highlight = self.observability_ax.scatter(self._days_for_all_rolls[yidx, xidx], self._theta[yidx, xidx], color='white', edgecolor='black', s=100)
        self._plot_overlay_elements.append(obs_highlight)
        obs_vline = self.observability_ax.axvline(self._days_for_all_rolls[yidx, xidx])
        self._plot_overlay_elements.append(obs_vline)
        obs_hline = self.observability_ax.axhline(self._theta[yidx, xidx])
        self._plot_overlay_elements.append(obs_hline)

        if self.result.companions[0]['separation'] != 0.0:
            c1_highlight = self.detector_ax.scatter(self.result.c1_x[yidx, xidx], self.result.c1_y[yidx, xidx], color='white', edgecolor='black', s=100)
            self._plot_overlay_elements.append(c1_highlight)
        if self.result.companions[1]['separation'] != 0.0:
            c2_highlight = self.detector_ax.scatter(self.result.c2_x[yidx, xidx], self.result.c2_y[yidx, xidx], color='white', edgecolor='black', s=100)
            self._plot_overlay_elements.append(c2_highlight)
        if self.result.companions[2]['separation'] != 0.0:
            c3_highlight = self.detector_ax.scatter(self.result.c3_x[yidx, xidx], self.result.c3_y[yidx, xidx], color='white', edgecolor='black', s=100)
            self._plot_overlay_elements.append(c3_highlight)

        scale_factor = (1/5) * np.average(self.result.scisize)


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
        self._mask_artist = None
        ax.set_aspect('equal')

        self.c1_plot_group = ax.scatter(self.result.c1_x, self.result.c1_y, picker=True, color=RED_GGPLOT)
        self.c2_plot_group = ax.scatter(self.result.c2_x, self.result.c2_y, picker=True, color=BLUE_GGPLOT)
        self.c3_plot_group = ax.scatter(self.result.c3_x, self.result.c3_y, picker=True, color=PURPLE_GGPLOT)

        ax.legend(
            (self.c1_plot_group, self.c2_plot_group, self.c3_plot_group),
            ('Companion 1', 'Companion 2', 'Companion 3')
        )

        x_sci_size, y_sci_size = self.result.scisize
        ax.set_xlim(-x_sci_size / 2, x_sci_size / 2)
        ax.set_ylim(-y_sci_size / 2, y_sci_size / 2)
        ax.set_xlabel('x (pixels from center)')
        ax.set_ylabel('y (pixels from center)')
        ax.scatter(self.result.s_x, self.result.s_y, facecolor='yellow', edgecolor='black', marker='*', s=100)

        self._overlay_mask(self.result.aperture.AperName)

    def _overlay_mask(self, apername):
        if self._mask_artist is not None:
            self._mask_artist.remove()
        arcsec_per_pixel = np.average(self.result.sciscale)  # avg of scale in x and y
        x_sci_size, y_sci_size = self.result.scisize
        if 'NRC' in apername and apername[-1] == 'R':
            if '210R' in apername:
                radius_arcsec = 0.40
            elif '335R' in apername:
                radius_arcsec = 0.64
            elif '430R' in apername:
                radius_arcsec = 0.82
            else:
                raise RuntimeError("Invalid mask!")
            # make a circle
            self._mask_artist = self.detector_ax.add_artist(patches.Circle((0, 0), radius=radius_arcsec / arcsec_per_pixel, alpha=0.5))
        elif 'NRC' in apername:
            if 'LWB' in apername:
                thin_extent_arcsec = 0.58 * (2 / 4)
                thick_extent_arcsec = 0.58 * (6 / 4)
            elif 'SWB' in apername:
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
            verts = np.concatenate([x_verts[:,np.newaxis], y_verts[:,np.newaxis]], axis=1)
            patch = patches.Polygon(verts, alpha=0.5)
            self._mask_artist = self.detector_ax.add_artist(patch)
        elif 'MIRI' in apername:
            y_angle = float(self.result.sciyangle)
            if 'LYOT' in apername:
                width_arcsec = 0.72
                radius_arcsec = 2.16
                circular_part = patches.Circle((0, 0), radius=radius_arcsec / arcsec_per_pixel)
                x_verts = width_arcsec / arcsec_per_pixel * np.array([-1, -1, 1, 1])
                y_verts = y_sci_size / 2 * np.array([-1, 1, 1, -1])
                x_verts_rot = np.cos(y_angle) * x_verts + np.sin(y_angle) * y_verts
                y_verts_rot = -np.sin(y_angle) * x_verts + np.cos(y_angle) * y_verts
                verts = np.concatenate([x_verts_rot[:,np.newaxis], y_verts_rot[:,np.newaxis]], axis=1)
                rectangular_part = patches.Polygon(verts)
                mask_collection = PatchCollection([rectangular_part, circular_part], alpha=0.5)
                self._mask_artist = self.detector_ax.add_artist(mask_collection)
            elif '1065' in apername or '1140' in apername or '1550' in apername:
                width_arcsec = 0.33
                x = x_sci_size / 2
                y = y_sci_size / 2
                maskwidth = width_arcsec / arcsec_per_pixel
                x_verts = np.array([-x, -maskwidth, -maskwidth, maskwidth, maskwidth, x, x, maskwidth, maskwidth, -maskwidth, -maskwidth, -x])
                y_verts = np.array([maskwidth, maskwidth, y, y, maskwidth, maskwidth, -maskwidth, -maskwidth, -y, -y, -maskwidth, -maskwidth])
                x_verts_rot = np.cos(y_angle) * x_verts + np.sin(y_angle) * y_verts
                y_verts_rot = -np.sin(y_angle) * x_verts + np.cos(y_angle) * y_verts
                verts = np.concatenate([x_verts_rot[:,np.newaxis], y_verts_rot[:,np.newaxis]], axis=1)
                mask = patches.Polygon(verts, alpha=0.5)
                self._mask_artist = self.detector_ax.add_artist(mask)
            else:
                raise RuntimeError("Invalid mask!")


def run():
    app = VisibilityCalculator()
    app.start()
