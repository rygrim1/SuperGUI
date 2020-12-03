import tkinter as tk
from tkinter import StringVar, IntVar, ttk
from collections import namedtuple
import re, glob, yaml
from os.path import expanduser 
import supermods.plotbands as plot_mods                                                                                                                                                                  
import supermods.Wyckoff as wyck
import supermods.col_key as col
from math import copysign
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from sys import platform as sys_pf
if sys_pf == 'darwin':
	matplotlib.use("TkAgg")

# Add K-path finder script that reads QE input file and determines k-path 
# plot of total energy convergence in real time from log file
# integrate into automation system to check congergence in real time, etc

class CheckPrimitive:
	def __init__(self, params):
		self.cell = params[1]
		self.space_group = params[0][1]
		self.param = namedtuple('Wyckoff_Parameter', ['atom', 'letter', 'x', 'y', 'z'])
#--------------------------------------------------------------------------------------------------		
	def constructPrim(self, validate):
		if validate=='cell':
			return self.cell
		else:
			wyckoff_params = [self.param(*data) for data in self.cell]
			self.cell.clear()
			[wyck.get_wyckoff(self.cell, self.space_group, p.letter, p.atom, p.x, p.y, p.z) for p in wyckoff_params]
			return self.cell
class SuperCell:
	def __init__(self, from_yaml, link, unit_cell, output, cell_dm, sort_keys, shifts, charges):
		self.primitive = unit_cell[1]
		self.name=unit_cell[0]
		self.format=output
		self.sort_by=sort_keys
		self.cell_dm=cell_dm
		self.shifts=shifts
		self.charges = charges
		self.invert = False
		self.invert_key = False
		self.origin = False
		self.zero_key = False
		self.from_yaml=from_yaml
		self.bilbao = link
#--------------------------------------------------------------------------------------------------
	def constructCell(self, X, Y, Z, x_shift, y_shift, z_shift):
		""" Input is self.primitive, dimensions X, Y, and Z, and the shifts for each dimension
		Returns the unsorted bulk cell
		"""
		iterated_cell = []
		# Core algorithm to produce supercell coordinate sets"""
		for i in range(len(self.primitive)):
			# Iterates individually upon values (atoms and coordinates) in self.primitive
			prim = self.primitive[i]
			atom = prim[0]
			# Creation of supercell basis in XYZ
			basis = [
				prim[1] / X,
				prim[2] / Y,
				prim[3] / Z, 
			]
			# Inner comprehension of all reduced coordinates in X
			X_red = [ 
				(
					basis[0] + i / X,
					basis[1],
					basis[2], 
				)
				for i in range(X)
			]
			# Inner comprehension of all reduced coordinates in Y   
			XY_red = [
				(
					xred[0],
					xred[1] + i / Y,
					xred[2],
				)
				for xred in X_red
				for i in range(Y)
			]
			# Final comprehension, adds reduced Z coordinates and any shift parameters
			XYZ_red = [
				(
					atom,
					round(xred[0] + x_shift, 7),
					round(xred[1] + y_shift, 7),
					round(xred[2] + i / Z + z_shift, 7),
				)
				for xred in XY_red for i in range(Z)
			]
			for xred in XYZ_red:
			# Updates iterated_cell with expanded coordinate sets for each individual atom from self.primitive until the entire
			# bulk cell is returned
				iterated_cell.append(xred)
		return iterated_cell
#--------------------------------------------------------------------------------------------------
	def layerCell(self, *args):
		""" Layer cell """
		if len(args)==0:
			bulk_cell = self.constructCell(*self.cell_dm, *self.shifts)
		if len(args)==1:
			bulk_cell=args[0]
		layered_cell = sorted(bulk_cell,
			key = lambda x: tuple(x[self.sort_by[i]] for i in range(3)))
		for idx, pos in enumerate(layered_cell):
			if idx < len(layered_cell)-1 and ((list(layered_cell[idx+1])[self.sort_by[0]] == 1-pos[self.sort_by[0]] or list(layered_cell[idx+1]) == [pos[0], -pos[1], -pos[2], -pos[3]])):
					self.invert = True
					self.invert_key = idx+1
					print('ok')
			elif idx < len(layered_cell)-1 and list(layered_cell[idx+1])[self.sort_by[0]] == -pos[self.sort_by[0]]:
					self.invert = True
					self.invert_key = idx+1
			if pos[1:4] == [0]*3:
				self.origin = True
				self.zero_key = idx+1
		return layered_cell
#--------------------------------------------------------------------------------------------------
	def displayParam(self, disp_id, X, Y, Z):
		name=[AnsiiCodes.cyan+str(s) for s in self.name]
		if disp_id ==1:
			spc=55
		else:
			spc=0
		if X == Y == Z == 1: 
			print('\n'+'*'*spc+'\n', StringFormats.cell_name.format(*name), '\n')
		else: 
			print('\n'+'*'*spc+'\n', StringFormats.supercell_name.format(*self.cell_dm, *name), '\n')
		if self.bilbao:
			print('Cell constructed from Wyckoff positions using standard ITA settings:\n%s' %(self.bilbao))
#--------------------------------------------------------------------------------------------------
	def displayCell(self, *args):
		if len(args)==0:
			layered_supercell = self.layerCell()
		if len(args)==1:
			layered_supercell = self.layerCell(args[0])
		self.displayParam(1, *self.cell_dm)
		print('\n', StringFormats.header.format('X','Y','Z'))
		for idx, position in enumerate(layered_supercell, 1):
			if self.invert and idx in [self.invert_key, self.invert_key + 1]:
				print(StringFormats.invert.format(idx, *position))
			if not self.invert or (self.invert and idx not in [self.invert_key, self.invert_key + 1]):
				print(StringFormats.main.format(idx, *position))
		if self.invert: 
			print(StringFormats.invert_notice.format(self.invert_key, self.invert_key + 1))
		if self.origin: 
			print(StringFormats.invert_notice.format(self.zero_key))
#--------------------------------------------------------------------------------------------------
	def cleaveCell(self, surface, *args):
		if len(args)==0:
			layers = self.layerCell()
		if len(args)==1:
			layers = args[0]
		if surface==False:
			surface=[1, len(layers)]
		if max(surface) > len(layers):
			print('Invalid index')
			return None
		self.displayParam(2, *self.cell_dm)
		cleaved_surface = layers[surface[0]-1:surface[1]]
		quantum_sort=sorted(cleaved_surface,
			key = lambda x: tuple(x[self.sort_by[i]] for i in range(3)))
		abi_sort=sorted(cleaved_surface,
			key = lambda x:(x[0],x[self.sort_by[0]]))
		for i in range(len(cleaved_surface)):
			if self.format == 'ABINIT':
				print(StringFormats.abinit.format(*abi_sort[i][1:4], abi_sort[i][0]))
			if self.format == 'Quantum ESP.':
				print(StringFormats.quantum.format(*quantum_sort[i]))
		try:
			net_charge = sum([self.charges.get(value[0], 'Null') for value in cleaved_surface])
			if net_charge > 0:
				print(StringFormats.charge.format(AnsiiCodes.bold, AnsiiCodes.green, '+' + str(net_charge), AnsiiCodes.end))
			elif net_charge == 0:
				print(StringFormats.charge.format(AnsiiCodes.bold, AnsiiCodes.cyan, net_charge, AnsiiCodes.end))
			else:
				print(StringFormats.charge.format(AnsiiCodes.bold, AnsiiCodes.red, net_charge, AnsiiCodes.end))
		except TypeError:
			print('\nSurface charge not available. Update atomic charges in cells.yaml.')
			pass
#--------------------------------------------------------------------------------------------------
class DropDown(ttk.OptionMenu):
	def __init__(self, parent, options: list, initial_value: str=None, style: str=None):
		self.var = tk.StringVar(parent)
		self.var.set(initial_value if initial_value else options[0])
		self.option_menu = ttk.OptionMenu.__init__(self, parent, self.var, initial_value, *options, style=style)
		self.callback = None
		
	def add_callback(self, callback: callable):
		def internal_callback(*args):
			callback()
		self.var.trace("w", internal_callback)
	def get(self):
		return self.var.get()
	def set(self, value: str):
		self.var.set(value)
#--------------------------------------------------------------------------------------------------
class SuperGui(tk.Tk):
	def __init__(self):
		tk.Tk.__init__(self)
		self.yaml_path=path=str(Path(__file__).parent.absolute())+'/cells.yaml'
		with open(self.yaml_path) as yaml_file:
			self.yaml_data = yaml.load(yaml_file, Loader = yaml.FullLoader)
		self.charges = self.yaml_data['Charges']
		self.struct = [key for key in self.yaml_data if key.lower() != 'charges']
		self.key_map = {'x': 1, 'y': 2, 'z': 3 }
		self.cell = []
		#self.tk.call('tk', 'scaling', 4.0)
		self.title("SuperGUI.py v1.a")
		self.geometry('875x625')
		self.config(bg='gray93')
		self.style = ttk.Style()
		self.note=ttk.Notebook(self)
		self.main=ttk.Frame(self.note)
		self.adsorb=ttk.Frame(self.note)
		self.plots=ttk.Frame(self.note)
		self.exit_gui=ttk.Frame(self.note)
		self.pdos_plots=ttk.Frame(self.note)
		self.note.pack(expand=1, fill="both")
		self.note.add(self.main, text='Main')
		self.note.add(self.adsorb, text='Adsorbates')
		self.note.add(self.plots, text='Bands & PDOS')
		self.note.add(self.pdos_plots, text='View plots')
		self.note.add(self.exit_gui, text='Quit')
		self.note.grid(row=0, column=0, columnspan=6, rowspan=10, sticky='NWES')

		self.style.theme_create( "yummy", parent='alt', settings={
        "TNotebook": {
        	"configure": {"background": "lavender", "tabmargins": [2, 5, 2, 0],}},
        "TNotebook.Tab": {
            "configure": {"padding": [5, 1], "background": "gray88", "font": ["serif", 12, "bold"]},
            "map":       {"background": [("selected", "gray93")], 
                          "expand": [("selected", [1, 1, 1, 0])], "foreground": [("selected", "blue2")] }},
        "TFrame": {
        	"configure": {"fieldbackground": "gray93", "background": "gray93"}},
        "TLabel": {
        	"configure": {"fieldbackground": "gray93", "background": "gray93", "font": ["serif", 12, "bold"]}},
        "Taki.TLabel": {
        	"configure": {"background": "gray93", "fieldbackground": "gray93", "font": ["serif", 12, "bold"]}},
        "TLabelframe.Label": {
        	"configure": {"background": "gray93", "fieldbackground": "gray93", "borderwidth": 2, 
        					"bordercolor": "black", "font": ["serif", 12, "bold"] }},
        "TEntry": {
        	"configure": {"fieldbackground": "ghost white", "border": 2, "bordercolor": "black", "padding": 2}},
        "DOS.TEntry": {
        	"configure": {"fieldbackground": "ghostwhite"}, "font": ["serif", 12, "bold"]},
        "XY.TEntry": {
        	"configure": {"fieldbackground": "ghostwhite"}, "font": ["serif", 12, "bold"]},
        "TButton": {
        	"configure": {"relief": "raised", "background": "lavender", "foreground": "navy", "anchor": "center", 
        					"font": ["serif", 12, "bold"], "borderwidth": 4, "bordercolor": "black"}},
        "Make.TButton": {
        	"configure": {"relief": "raised", "background": "pink", 
        				  "foreground": "navy", "anchor": "center", "font": ["serif", 12, "bold"], "borderwidth": 2}},
        "Col.TButton": {
        	"configure": {"relief": "raised", 
        				  "foreground": "black", "anchor": "center", "font": ["serif", 12, "bold"], "borderwidth": 4}},
        "Toolbutton": {
        	"configure": {"font": ["serif", 12, "bold"], "borderwidth": 4, "relief": "raised",
        				  "background": "lavender", "anchor": "center", "bordercolor": "navy"}},
        "Invert.Toolbutton": {
        	"configure": {"relief": "raised", "background": "lavender", "borderwidth": 4, "font": ["serif", 12, "bold"], 
        					"padding": 2, "anchor": "center"}},
        "Plot.TButton": {
        	"configure": {"font": ["serif", 12, "bold"], "relief": "raised",
        				  "background": "lavender", "foreground": "red2", "borderwidth": 4, "anchor": "center"}},
        "Save.TButton": {
        	"configure": {"font": ["serif", 12, "bold"], "relief": "raised",
        				  "background": "lavender", "foreground": "red2", "borderwidth": 4, "anchor": "center"}},
        "TMenubutton": {
        	"configure": {"relief": "raised", "background": "lavender",  
        				  "width": 12, "borderwidth": 4, "anchor": "center", "font": ["serif", 12, 'bold']}},
        "DPI.TMenubutton": {
        	"configure": {"relief": "raised", "background": "lavender",  
        				  "width": 8, "anchor": "center", "font": ["serif", 12, "bold"]}},
        "Format.TMenubutton": {
        	"configure": {"relief": "raised", "background": "lavender",  
        				  "width": 10, "anchor": "center", "font": ["serif", 12, "bold"]}},
        "Bravais.TMenubutton": {
        	"configure": {"relief": "raised", "background": "lavender",  
        				  "width": 11, "anchor": "center", "font": ["serif", 12, "bold"]}},
        "Palette.TMenubutton": {
        	"configure": {"relief": "raised", "background": "lavender",  
        				  "width": 13, "borderwidth": 4, "anchor": "center", "font": ["serif", 12, 'bold']}},
        "Dim.TMenubutton": {
        	"configure": {"foreground": "black", "font": ["serif", 12, "bold"], "width": 10}},
        "Treeview": {
        	"configure": {"fieldbackground": "mint cream", "background": "mint cream", "font": ["menlo", 13]}, 
        				  "map": {"background": [("selected", "cyan")]}},
        "Treeview.Heading": {
        	"configure":{"font": ["serif", 12, "bold"]}
        			}})
		self.style.layout("Dim.TEntry",
                   [('Entry.plain.field', {'children': [(
                       'Entry.background', {'children': [(
                           'Entry.padding', {'children': [(
                               'Entry.textarea', {'sticky': 'nswe'})],
                      'sticky': 'nswe'})], 'sticky': 'nswe'})],
                      'border':'2', 'sticky': 'nswe'})])
		
		self.style.layout("Adsorb.TEntry",
                   [('Entry.plain.field', {'children': [(
                       'Entry.background', {'children': [(
                           'Entry.padding', {'children': [(
                               'Entry.textarea', {'sticky': 'nswe'})],
                      'sticky': 'nswe'})], 'sticky': 'nswe'})],
                      'border':'2', 'sticky': 'nswe'})])

		self.style.layout("Cleave.TEntry",
                   [('Entry.plain.field', {'children': [(
                       'Entry.background', {'children': [(
                           'Entry.padding', {'children': [(
                               'Entry.textarea', {'sticky': 'nswe'})],
                      'sticky': 'nswe'})], 'sticky': 'nswe'})],
                      'border':'2', 'sticky': 'nswe'})])
		self.style.layout("Path.TEntry",
                   [('Entry.plain.field', {'children': [(
                       'Entry.background', {'children': [(
                           'Entry.padding', {'children': [(
                               'Entry.textarea', {'sticky': 'nswe'})],
                      'sticky': 'nswe'})], 'sticky': 'nswe'})],
                      'border':'2', 'sticky': 'nswe'})])
		self.style.theme_use("yummy")

		# Frame configuration
		for i in range(1):
			self.columnconfigure(i, weight=1)
		for i in range(1):
			self.rowconfigure(i, weight=1)
		for i in range(0, 5):
			self.main.columnconfigure(i, weight=1)
			self.adsorb.columnconfigure(i, weight=1)
		for i in range(0,7):
			self.plots.columnconfigure(i, weight=1)
		for i in range(1, 12):
			self.main.rowconfigure(i, weight=1)
			self.exit_gui.rowconfigure(i, weight=1)
			self.pdos_plots.columnconfigure(i, weight=1)
			self.pdos_plots.rowconfigure(i, weight=1)
			self.adsorb.rowconfigure(i, weight=1)
			self.exit_gui.columnconfigure(i, weight=1)
			self.plots.rowconfigure(i, weight=1)

		# main
		self.X = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.X.grid(row=2, column=1, sticky='nwes', padx=2, pady=4)
		self.X.insert(0, 1)
		self.Y = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.Y.grid(row=2, column=2, sticky='nsew', padx=2, pady=4)
		self.Y.insert(0, 1)
		self.Z = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.Z.grid(row=2, column=3, sticky='nsew', padx=2, pady=4)
		self.Z.insert(0, 1)

		self.x_shift = ttk.Entry(self.main, width=4, justify="center", style='Dim.TEntry')
		self.x_shift.grid(row=3, column=1, sticky='nsew', padx=2, pady=4)
		self.x_shift.insert(0, 0.0)
		self.y_shift = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.y_shift.grid(row=3, column=2, sticky='nsew', padx=2, pady=4)
		self.y_shift.insert(0, 0.0)
		self.z_shift = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.z_shift.grid(row=3, column=3, sticky='nsew', padx=2, pady=4)
		self.z_shift.insert(0, 0.0)

		self.x_sort = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.x_sort.grid(row=4, column=1, sticky='nsew', padx=2, pady=4)
		self.x_sort.insert(0, 'x')
		self.y_sort = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.y_sort.grid(row=4, column=2, sticky='nsew', padx=2, pady=4)
		self.y_sort.insert(0, 'y')
		self.z_sort = ttk.Entry(self.main, width=4,  justify="center", style='Dim.TEntry')
		self.z_sort.grid(row=4, column=3, sticky='nsew', padx=2, pady=4)
		self.z_sort.insert(0, 'z')

		self.cleave_a = ttk.Entry(self.main, width=4,  justify="center", style='Cleave.TEntry')
		self.cleave_a.grid(row=5,column=1,sticky='nsew',padx=2, pady=4)
		self.cleave_b = ttk.Entry(self.main, width=4,  justify="center", style='Cleave.TEntry')
		self.cleave_b.grid(row=5,column=3,sticky='nsew',padx=2, pady=4)

		self.path= ttk.Entry (self.main, justify="center", style='Path.TEntry')
		self.path.grid(row=9, column=1, sticky='nsew', padx=2, pady=4, columnspan=3)
		self.path.insert(0,'')

		# main labels
		self.xyz = ["X", "Y", "Z"]
		self.ad_labels= ["Atom", "Distance (\u212B)", "xpa"]
		for i in range(1,4):
			self.test = ttk.Label(self.main, text=self.xyz[i-1], font=("serif", 14, 'bold')).grid(row=1, column=i)
			self.test = ttk.Label(self.adsorb, text=self.ad_labels[i-1], font=("serif", 12, 'bold')).grid(row=8, column=i-1)
		self.dim_label=ttk.Label(self.main, text='Supercell dimensions').grid(row=2, column=0)
		self.shifts_label =ttk.Label(self.main, text='Origin shifts').grid(row=3, column=0)
		self.sort_label =ttk.Label(self.main, text='Sorting priority', style='Taki.TLabel').grid(row=4, column=0)
		self.cleave_atoms=ttk.Label(self.main, text='Cleave surface', style='Taki.TLabel').grid(row=5, column=0)
		self.cleave_arrow=ttk.Label(self.main, text="\u2192", font=("serif", 20, "bold")).grid(row=5, column=2)
		self.lattice_params=ttk.Label(self.main, text='Lattice constants',font=("serif",12, "bold")).grid(row=6, column=0)
		self.angles=ttk.Label(self.main, text='\u03B1, \u03B2, \u03B3',font=("serif",12, "bold")).grid(row=7, column=0)
		#self.miller_indices=ttk.Label(self.main, text='Miller indices',font=("serif",12, "bold")).grid(row=8, column=0)
		self.N = StringVar(self.adsorb, value='N atoms: 0')
		self.n_atoms=ttk.Label(self.adsorb, textvariable=self.N, width=10, foreground='red2', font=("serif", 12, 'bold'), justify="center")
		self.n_atoms.grid(row=8, column=3)

		# main buttons
		self.create_cell = ttk.Button(self.main, text='Construct cell', 
			width=12, style='TButton', command=lambda: self.validate_params(1)).grid(row=2, column=4, ipady=4)
		self.create_cell = ttk.Button(self.adsorb, text='Construct cell', 
			width=12, style='TButton', command=lambda: self.validate_params(4)).grid(row=9, column=5, ipady=4)
		
		self.exit_button=ttk.Button(self.main, text='Quit', width=10, style="TButton", command=self.quit)
		self.exit_button.grid(row=4, column=4, ipady=4)

		self.create_file=ttk.Button(self.main, text='Create pw.x file', 
			style = 'Make.TButton', width=12).grid(row=9, column=4, ipady=4)
		self.cleave_button=ttk.Button(self.main, text='Cleave surface', 
			style = 'TButton', command=lambda: self.validate_params(2), width=12).grid(row=3, column=4, ipady=4)
		self.add_button=ttk.Button(self.adsorb, text='Add adsorbate(s)',
			command=lambda: self.validate_params(5), style = 'TButton', width=13).grid(row = 10, column=2, ipady=4)
		self.add_button=ttk.Button(self.adsorb, text='Print cell',
			command=lambda: self.validate_params(1, self.cell), style = 'TButton', width=10).grid(row=10, column=3, ipady=4)
		self.clear_menu=ttk.Button(self.adsorb, text='Clear menu',
			command=lambda: self.validate_params(1, self.cell), style = 'TButton', width=12).grid(row = 9, column=4, ipady=4)
		# main checkbuttons
		self.check_var = IntVar()
		self.pwx_inputs=ttk.Checkbutton(self.main, width=12, command=self.validate_path,
			var=self.check_var, text='Load pw.x file', style="Toolbutton").grid(row=9, column=0, ipady=4)
		self.inv_var = IntVar()
		self.inversion_pairs=ttk.Checkbutton(self.adsorb, text='Inversion sym', command=self.validate_path, 
			var=self.inv_var, style="Invert.Toolbutton", width=11, ).grid(row=10, column=1, ipady=2)
		self.alat = ttk.Entry(self.main, width=10, justify="center", style='Adsorb.TEntry')
		self.alat.grid(row=6, column=1, sticky='nsew', padx=2, pady=4)
		self.alat.insert(0, 'a')
		self.blat = ttk.Entry(self.main, width=10, 
			justify="center", style='Adsorb.TEntry')
		self.blat.grid(row=6, column=2, sticky='nsew', padx=2, pady=4)
		self.blat.insert(0, 'b')

		self.clat = ttk.Entry(self.main, width=10, 
			justify="center", style='Adsorb.TEntry')
		self.clat.grid(row=6, column=3, sticky='nsew', padx=2, pady=4)
		self.clat.insert(0, 'c')
		self.alpha = ttk.Entry(self.main, width=10, 
			justify="center", style='Adsorb.TEntry')
		self.alpha.insert(0,'\u03B1')
		self.alpha.grid(row=7, column=1, sticky='nsew', padx=2, pady=4)
		self.beta = ttk.Entry(self.main, width=10, 
			justify="center", style='Adsorb.TEntry')
		self.beta.grid(row=7, column=2, sticky='nsew', padx=2, pady=4)
		self.beta.insert(0, '\u03B2')
		self.gamma = ttk.Entry(self.main, width=10, 
			justify="center", style='Adsorb.TEntry')
		self.gamma.grid(row=7, column=3, sticky='nsew', padx=2, pady=4)
		self.gamma.insert(0,'\u03B3')
		self.sorted_dos=[]
		self.plot_fermi = IntVar(value=0)
		self.ticks_var = IntVar(value=0)
		self.plot_kpath = IntVar()
		self.plot_bg = IntVar(value=0)
		self.fill_pdos = IntVar(value=1)
		self.plot_shift = IntVar(value=0)

		self.pdos_ticks=tk.Checkbutton(self.plots, variable=self.ticks_var, highlightthickness=0, bd=4, relief='flat', text='PDOS ticks', justify="center", width=11, bg='gray93', font=("serif", 12), 
			command=lambda: self.plot_params(1))
		self.pdos_ticks.grid(row=3, column=5)
		self.fermi=tk.Checkbutton(self.plots, variable=self.plot_fermi, text='Fermi level', justify="center", width=12, bg='gray93', font=("serif", 12), 
			command=lambda: self.plot_params(2))
		self.fermi.grid(row=2, column=6, columnspan=1, padx=(0,0))
		self.fermi_shift=tk.Checkbutton(self.plots, variable=self.plot_shift, text='Fermi shift', justify="center", width=12, bg='gray93', font=("serif", 12), 
			command=lambda: self.plot_params(8))
		self.fermi_shift.grid(row=3, column=6, columnspan=1, padx=(0,0))
		self.bandgap=tk.Checkbutton(self.plots, variable=self.plot_bg, text='Band gap', justify="center", width=11, bg='gray93', font=("serif", 12), 
			command=lambda: self.plot_params(3))
		self.bandgap.grid(row=2, column=4, padx=(0,0))
		self.kpath=tk.Checkbutton(self.plots, variable=self.plot_kpath, text='K-path', justify="center", width=11, bg='gray93', font=("serif", 12), 
			command=lambda: self.plot_params(4))
		self.kpath.grid(row=3, column=4, padx=(0,0))
		self.filldos=tk.Checkbutton(self.plots, variable=self.fill_pdos, text='Fill PDOS', justify="center",  width=11, bg='gray93', font=("serif", 12),  fg="green", 
			command=lambda: self.plot_params(5))
		self.filldos.grid(row=2, column=5, padx=(0,0))

		#self.adsorbates=tk.Checkbutton(self.plots, variable=self.plot_fermi, text='Adsorbates', justify="center", width=12, bg='gray93', font=("serif", 12), 
		#	command=lambda: self.plot_params(2))
		#self.adsorbates.grid(row=6, column=0, columnspan=1, padx=(0,0))
		#self.sel_ads=ttk.Button(self.plots, text='Select adsorbate', style = 'Col.TButton', width=14, 
		#	command=self.change_col).grid(row=6, rowspan=1, column=1, columnspan=1, ipady=2, padx=(0,0), sticky='w')
		#self.clear_dos=ttk.Button(self.plots, width=10, text='Clear items', command=self.clear_pdos, style="TButton")
		#self.clear_dos.grid(row=5, column=4, ipadx=4, ipady=3)

		#Plot labels
		self.eV_range=ttk.Label(self.plots, width=15, text='Energy range (eV)',
			font=("serif", 12)).grid(row=3, column=0, ipadx=0, columnspan=1, sticky='e')
		self.x_range=ttk.Label(self.plots, width=15, text='PDOS range (arb)',
			font=("serif", 12)).grid(row=2, column=0, ipadx=0, sticky='e')

		self.plot_var=tk.IntVar(value=4)
		self.select_bands=tk.Radiobutton(self.plots, text='Band structure', variable=self.plot_var, value=2, bg='gray93', selectcolor='green', font=("serif", 12))
		self.select_bands.grid(row=1, column=4, columnspan=1, ipadx=0, padx=(0,0))
		self.select_PDOS=tk.Radiobutton(self.plots, text='PDOS', variable=self.plot_var, value=3, bg='gray93', selectcolor='green', font=("serif", 12))
		self.select_PDOS.grid(row=1, column=5, ipadx=0, columnspan=1)
		self.select_both=tk.Radiobutton(self.plots, text='Combined', variable=self.plot_var, value=4, selectcolor='green', bg='gray93', font=("serif", 12))
		self.select_both.grid(row=1, column=6, columnspan=1, ipadx=0, padx=(0,0))

		self.plot_text = tk.StringVar()
		self.plot_text.set('Create plot')

		self.create_plot=ttk.Button(self.plots, command=lambda: self.plot_params(6), 
			style="Plot.TButton", width=10, textvariable=self.plot_text)
		self.create_plot.grid(row=9, column=6, ipady=2, columnspan=1, padx=(0,10))

		# PDOS plot
		self.figure=plt.figure(figsize=(8,6))
		self.figure.patch.set_facecolor('whitesmoke')
		self.plt_status=True
		# Bands plot
		self.bandplot = FigureCanvasTkAgg(self.figure, self.pdos_plots)
		self.bandplot.get_tk_widget().grid(row=0, column=0, rowspan=10, columnspan=10, sticky='nsew')
		self.bandplot.draw()
		self.ax=plot_mods.clear_plot(self.figure)

		#self.plt_status=True
		self.save_plt=ttk.Button(self.plots, width=10, command=lambda: self.plot_params(6,7), text='Save plot', style="Save.TButton")
		self.save_plt.grid(row=10, column=6, ipady=2, padx=(0,10))
		self.plt_dpi= DropDown(self.plots, [300, 600, 1200], 'DPI', style='DPI.TMenubutton')
		self.plt_dpi.grid(row=10, column=5, ipady=2, pady=(0,0), padx=(0,0), sticky='w')
		self.plt_fmt= DropDown(self.plots, ['png', 'pdf', 'ps', 'eps'], 'File format', 'Format.TMenubutton')
		self.plt_fmt.grid(row=10, rowspan=1, ipady=2, column=4, padx=(0,5), sticky='e')
		self.eV_a = ttk.Entry(self.plots, width=8, justify="center", style='Cleave.TEntry')
		self.eV_a.grid(row=3, column=1, sticky='')
		self.eV_a.insert(0, -10)
		self.eV_b = ttk.Entry(self.plots, width=8, justify="center", style='Cleave.TEntry')
		self.eV_b.grid(row=3, column=2, padx=(0,0), sticky='')
		self.eV_b.insert(0, 10)
		self.dos_a = ttk.Entry(self.plots, width=8, justify="center", style='Cleave.TEntry')
		self.dos_a.grid(row=2, column=1, sticky='', ipadx=0)
		self.dos_b = ttk.Entry(self.plots, width=8, justify="center", style='Cleave.TEntry')
		self.dos_b.grid(row=2, column=2, sticky='', ipadx=0)
		self.dos_a.insert(0, 0)
		self.dos_b.insert(0, 2)

		self.plotpathlabel=ttk.Label(self.plots, width=12, text='Directory path',
			font=("serif", 12)).grid(row=0, column=0, padx=(0,0), sticky='e')
		self.plotpath= ttk.Entry (self.plots, width=35, font=("serif", 13), 
			justify="left", style='Cleave.TEntry')
		self.plotpath.grid(row=0, column=1, ipady=2, columnspan=4, pady=(5,0), padx=(15,0), sticky='w')
		self.plotpath.insert(0,' ~/SuperGUI/kaol.bulk.wfc/')

		self.plot_label =ttk.Label(self.plots, text='Plot title', width=8, font=("serif", 12)).grid(row=1, column=0, padx=(0,0), sticky='e')
		self.plot_title_entry= ttk.Entry (self.plots, font=("serif", 12), width=15, justify="left", style='Cleave.TEntry')
		self.plot_title_entry.grid(row=1, column=1, columnspan=2, ipady=2, padx=(15,0), ipadx=0, sticky='w')
		self.plot_title_entry.insert(0, ' ')

		self.pdos_tree = ttk.Treeview(self.plots, columns=('natom', 'Orbital'), selectmode='none')
		self.pdos_tree.column('natom', stretch='no', minwidth=60, width=60, anchor='center')
		self.pdos_tree.column('Orbital', stretch='no', minwidth=60, width=90, anchor='center')
		self.pdos_tree.heading('natom', text='pw.x #')
		self.pdos_tree.heading('Orbital', text='Orbital')

		self.pdos_tree.grid(row=5, column=1, columnspan=2, rowspan=4, pady=(0,0), padx=(10,0), sticky='nes')
		self.pdos_tree['show'] = 'headings'
		self.pdos_tree.bind("<ButtonRelease-1>", self.pdos_select)
		self.ads_tree = ttk.Treeview(self.plots, columns=('Adsorbate'), selectmode='none')
		self.ads_tree.column('Adsorbate', stretch='no', minwidth=60, width=110, anchor='center')
		self.ads_tree.heading('Adsorbate', text='Adsorbates')
		self.ads_tree.grid(row=5, column=6, columnspan=1, rowspan=4, pady=(0,0), padx=(0,10), ipadx=5, sticky='news')
		self.ads_tree['show'] = 'headings'
		self.ads_tree.bind("<ButtonRelease-1>", self.pdos_select)

		self.label_tree = ttk.Treeview(self.plots, columns=('Labels'), selectmode='none')
		self.label_tree.column('Labels', stretch='no', minwidth=60, width=135, anchor='center')
		self.label_tree.heading('Labels', text='Legend label')
		self.label_tree.grid(row=5, column=5, columnspan=1, rowspan=4, pady=(0,0), padx=(0,0), ipadx=0, sticky='news')
		self.label_tree['show'] = 'headings'
		self.label_tree.bind("<ButtonRelease-1>", self.label_select)
		# ADSORBATE MENU
		self.format_menu = DropDown(self.main, ['Quantum ESP.', 'ABINIT'])
		self.format_menu.grid(row = 1, column=4, ipady=2)
		self.unit_menu = DropDown(self.main, ['Angstroms', 'Bohr'], 'Angstroms', 'Dim.TMenubutton')
		self.unit_menu.grid(row = 6, column=4, ipady=2)
		#self.sym_menu = DropDown(self.plots, ['Bravais lattice', 'Tetragonal'], 'Bravais lattice', 'Palette.TMenubutton')
		#self.sym_menu.grid(row = 8 , column=0, columnspan=1,  padx=(15,0), ipady=2, pady=(0,0), sticky='n')
		self.struct_menu = DropDown(self.main, self.struct, 'Structure')
		self.struct_menu.grid(row = 1, column=0, ipady=4, ipadx=2)
		self.adsorb_dim= DropDown(self.adsorb, ['x', 'y', 'z'], 'Direction', 'Dim.TMenubutton')
		self.adsorb_dim.grid(row=10, column=0, ipady=2)
		self.bands_lw_label=ttk.Label(self.plots, width=16, anchor="center", text='Band linewidth', style='Taki.TLabel',
			font=("serif", 12)).grid(row=9, column=0, columnspan=1, padx=(5,0), sticky='e')
		self.dos_lw_label=ttk.Label(self.plots, width=16, anchor="center", text='DOS linewidth',
			font=("serif", 12)).grid(row=9, column=3, columnspan=1, padx=(5,0), sticky='e')
		self.lw_bands=tk.Scale(self.plots, from_=0.0, to=2.0, digits = 3, tickinterval=0.50, borderwidth=4, 
			resolution = 0.05, width=20, length=185, orient='horizontal', bg='gray93', troughcolor="ghostwhite", 
			font=("serif", 12, 'bold'))
		self.lw_dos=tk.Scale(self.plots, length=185, from_=0.0, to=2.0, digits = 3, tickinterval=0.50, resolution = 0.05, 
			borderwidth=4, width=20, orient='horizontal', bg='gray93', troughcolor="ghostwhite", 
				font=("serif", 12, 'bold'))
		self.lw_dos.grid(row=9, column=4, columnspan=2, padx=(0,5))
		self.lw_dos.set(1.25)
		self.lw_bands.grid(row=9, column=1, columnspan=2, padx=(10,5))
		self.lw_bands.set(0.50)

		self.color_choice= DropDown(self.plots, col.color_key(1), 'Colors', 'Palette.TMenubutton')
		self.color_choice.grid(row = 5, rowspan=1, column=0, columnspan=1, padx=(15,0), ipady=2, sticky='')

		self.color_palettes= DropDown(self.plots, ['Default'], 'Color palettes', 'Palette.TMenubutton')
		self.color_palettes.grid(row = 7, rowspan=1, column=0, columnspan=1, ipadx=1, padx=(15,0), ipady=2, sticky='')
		self.change_color=ttk.Button(self.plots, text='Change orbital color', style = 'Col.TButton', width=14, 
			command=self.change_col).grid(row=6, rowspan=1, column=0, columnspan=1, ipady=2, padx=(15,0), sticky='nw')


		self.load_palette=ttk.Button(self.plots, text='Load color palette', style = 'Col.TButton', width=14, 
			command=self.change_col).grid(row=8, rowspan=1, column=0, columnspan=1, ipady=2, padx=(15,0), sticky='n')
		self.save_palette=ttk.Button(self.plots, text='Save palette', style = 'Col.TButton', width=11, 
			command=self.change_col).grid(row=10, rowspan=1, column=2, columnspan=1, ipady=2, padx=(5,0), sticky='w')
		self.tree = ttk.Treeview(self.adsorb, columns=('natom', 'Atom','x','y','z'), selectmode='none')
		self.new_palette_label=ttk.Label(self.plots, width=15, anchor="center", text='Save palette as',
			font=("serif", 12)).grid(row=10, column=0, columnspan=1, padx=(0,0), sticky='e')
		self.new_palette_entry = ttk.Entry(self.plots, width=12, font=("serif", 13), justify="left", style='Cleave.TEntry')
		self.new_palette_entry.grid(row=10, column=1, ipady=0, padx=(0,0), sticky='w')
		self.new_palette_entry.insert(0, ' new palette')

		self.leg_label =ttk.Label(self.plots, text='Add legend label', width=13, font=("serif", 13)).grid(row=4, column=0, padx=(0,0), sticky='e')
		self.leg_entry = ttk.Entry(self.plots, width=10, font=("serif", 13), justify="left", style='Cleave.TEntry')
		self.leg_entry.grid(row=4, column=1, padx=(13,0), sticky='w')

		self.change_label=ttk.Button(self.plots, text='Swap label', style = 'Col.TButton', width=8, 
			command=self.change_label).grid(row=4, rowspan=1, column=2, columnspan=1, ipady=2, padx=(5,0), sticky='w')

		self.tree.column('natom', stretch='yes', minwidth=0, width=5, anchor='center')
		self.tree.column('Atom', stretch='yes', minwidth=0, width=10, anchor='center')
		self.tree.column('x', stretch='yes', minwidth=0, width=5, anchor='center')
		self.tree.column('y', stretch='yes', minwidth=0, width=5, anchor='center')
		self.tree.column('z', stretch='yes', minwidth=0, width=5, anchor='center')
		self.tree.heading('natom', text='#')
		self.tree.heading('Atom', text='Atom')
		self.tree.heading('x', text='x')
		self.tree.heading('y', text='y')
		self.tree.heading('z', text='z')
		self.tree.grid(row=0, column=0, columnspan=6, sticky="NSEW", padx=15, rowspan=8, pady=10)
		self.tree['show'] = 'headings'
		self.tree.bind("<ButtonRelease-1>", self.adsorb_select)

		self.color_tree = ttk.Treeview(self.plots, columns=('Atom', 'DOS color'), selectmode='none')
		self.color_tree.column('Atom', stretch='no', width=90, anchor='center')
		self.color_tree.column('DOS color', stretch='no', width=165, anchor='center')
		self.color_tree.heading('Atom', text='Orbital')
		self.color_tree.heading('DOS color', text='DOS color')
		self.color_tree.grid(row=5, column=3, columnspan=2, rowspan=4, padx=(0,0), pady=(0,0), sticky='news')
		self.color_tree['show'] = 'headings'
		self.color_tree.bind("<ButtonRelease-1>", self.file_select)

		self.ad_deselect=ttk.Button(self.adsorb, text='Deselect all',
			command=self.ad_deselect, style = 'TButton', width=12)
		self.ad_deselect.grid(row=10, column=4, ipady=4)

		self.find_dir=ttk.Button(self.plots, width=10, 
			command=lambda: self.insert_dosmenu(), text='Load wfc files', style="Col.TButton")
		self.find_dir.grid(row=0, column=5, ipady=2, columnspan=1, padx=(0,5), pady=(5,0), sticky='e')
		self.clear_fil=ttk.Button(self.plots, width=10, 
			command=lambda: self.insert_dosmenu(), text='Clear files', style="Col.TButton")
		self.clear_fil.grid(row=0, column=6, ipady=2, columnspan=1, padx=(0,5), pady=(5,0))

		self.create_file=ttk.Button(self.adsorb, text='Create pw.x file', 
			style = 'Make.TButton', width=12).grid(row=10, column=5, ipady=4)

		self.del_atom=ttk.Button(self.adsorb, text='Delete atom(s)',
			command=self.ad_delete, style ='TButton', width=11)
		self.del_atom.grid(row=9, column=3, ipady=4)

		self.dim = [self.X, self.Y, self.Z]
		self.shifts = [self.x_shift, self.y_shift, self.z_shift]
		self.sort = [self.x_sort, self.y_sort, self.z_sort]
		self.celldm=[self.alat.get(), self.blat.get(), self.clat.get()]
		self.enter_adsorb = ttk.Entry(self.adsorb, width=8, font=("serif", 12), justify="center", style='TEntry')

		self.enter_adsorb.grid(row=9, column=0, ipady=4)
		self.enter_adsorb.insert(0, 'H')
		self.angs = ttk.Entry(self.adsorb, width=8, font=("serif", 12), justify="center", style='TEntry')

		self.angs.grid(row=9, column=1, ipady=4)
		self.angs.insert(0, 1.0)
		self.xpa = ttk.Entry(self.adsorb, width=8, font=("serif", 12), justify="center", style='TEntry')
		self.ad_index=[]
		self.xpa.grid(row=9, column=2, ipady=4)
		self.xpa.insert(0, 'auto')
		self.retain=[]
		self.QE_sort=[]
		self.dos_col={}
		self.dos_dict={}
		self.validate_path()
		self.style.configure("Dim.TEntry", padding=2, fieldbackground="SlateGray1", 
			fg='black', anchor='left')
		self.style.configure("Adsorb.TEntry", padding=2, fieldbackground="lavender", 
			bordercolor="black", anchor='center')
		self.style.configure("Cleave.TEntry", padding=2, fieldbackground="ghost white", 
			bordercolor="black", anchor='center')
		self.style.configure("Checkbutton", relief="flat", offrelief='flat')
		self.save=False
#--------------------------------------------------------------------------------------------------
	def adsorb_select(self, event):
		self.tree.selection_toggle(self.tree.focus())
		self.ad_index=max([self.tree.item(sel)['values'][0] for sel in self.tree.get_children()])
	def ad_deselect(self):
		for sel in self.tree.selection():
			self.tree.selection_remove(sel)
	def ad_delete(self):
		for sel in self.tree.selection():
			atom=self.tree.item(sel)['values'][1:]
			atom=[atom[0], float(atom[1]), float(atom[2]), float(atom[3])]
			self.tree.delete(sel)
			self.cell.remove(atom)
			if atom in self.retain:
				self.retain.remove(atom)
		self.N.set('N atoms: '+str(len(self.cell)))
	def insert_admenu(self, cell):
		self.cell=[[v[0], round(v[1],6), round(v[2],6), round(v[3],6)] for v in cell]
		for atom in self.cell:
			self.tree.insert("", 'end', iid=self.cell.index(atom), text=atom, values=(self.cell.index(atom)+1, *atom))
	def clear_pdos(self):
		self.color_tree.delete(*self.color_tree.get_children())
		self.pdos_tree.delete(*self.pdos_tree.get_children())
	def label_select(self, event):
		self.label_tree.selection_toggle(self.label_tree.focus())
		index=self.label_tree.focus()
		label=self.label_tree.item(index)['values'][0]
	def pdos_select(self, event):
		self.pdos_tree.selection_toggle(self.pdos_tree.focus())
		index=self.pdos_tree.focus()
		pdos_vals=self.pdos_tree.item(index)['values']
		if index not in self.color_tree.get_children():
			self.dos_dict.update({index: self.dos_col[index]})
			color=self.dos_dict[index][1]
			col_name=col.color_key(5, color)
			self.color_tree.insert("", 'end', iid=index, text=pdos_vals[2], values=(index, col_name, pdos_vals[2]))
			self.color_tree.item(index, tags=index.replace(" ",""))
			self.color_tree.tag_configure(index.replace(" ",""), foreground=(color))
			self.label_tree.insert("", 'end', iid=index, text=index, values=(index, index))
		for orb in self.color_tree.get_children():
			iid=self.color_tree.item(orb)['values'][0]
			if iid not in self.pdos_tree.selection():
				self.color_tree.delete(iid)
				del self.dos_dict[iid]
		for orb in self.label_tree.get_children():
			iid=self.label_tree.item(orb)['values'][1]
			if iid not in self.pdos_tree.selection():
				self.label_tree.delete(iid)
	def pdos_deselect(self):
		for sel in self.pdos_tree.selection():
			self.pdos_tree.selection_remove(sel)
	def label_deselect(self):
		for sel in self.label_tree.selection():
			self.label_tree.selection_remove(sel)
	def change_col(self):
		for sel in self.color_tree.selection():
			color=self.color_choice.get()
			hex_col=col.color_key(4, color)
			atom=self.color_tree.item(sel)['values'][2]
			iid=self.color_tree.item(sel)['values'][0]
			self.color_tree.tag_configure(iid.replace(" ",""), foreground=(hex_col))
			self.color_tree.item(iid, values=(iid, color, atom))
			self.dos_dict[iid][1]=hex_col
			self.color_tree.selection_remove(sel)
	def change_label(self):
		for sel in self.label_tree.selection():
			label=self.leg_entry.get()
			iid=self.label_tree.item(sel)['values'][1]
			self.label_tree.item(iid, values=(label, iid))
			self.label_tree.selection_remove(sel)
			self.dos_dict[iid][4]=label

	def insert_dosmenu(self):
		self.dos_col.clear()
		self.sorted_dos.clear()
		self.pdos_tree.delete(*self.pdos_tree.get_children())
		self.color_tree.delete(*self.color_tree.get_children())
		self.label_tree.delete(*self.label_tree.get_children())
		DOS_files = glob.glob(expanduser(self.plotpath.get().strip(" ")) + "/*wfc*")
		if not DOS_files:
			return None
		scf=glob.glob(expanduser(self.plotpath.get().strip(" ")) + "/*scf*")
		if not scf:
			self.plot_fermi.set(0)
			self.plot_bg.set(0)
			self.plot_shift.set(0)
		else:
			self.fermi.configure(fg="green")
			self.plot_fermi.set(1)
		for file in DOS_files:
			try:
				atom = file.split("atm")[1]
				index=re.findall(r'%s(\d+)' % '#', atom)
				n=index[1]
				st_atom=re.sub('[()]', '', re.findall('\(.*?\)',atom)[0])
				orbital=st_atom+' '+'#'+n+' '+re.findall('\(.*?\)',atom)[1][1]
				color=col.color_key(2, st_atom)
				orb_id=(index[0]+' '+orbital).replace(" ","")
				self.dos_col.update({orb_id: [st_atom, color, file, orb_id, orb_id]})
				self.sorted_dos.append([int(index[0]), orbital, (index[0]+' '+orbital).replace(" ",""), st_atom, file, color])
			except IndexError:
				continue
		self.sorted_dos=sorted(self.sorted_dos)
		for orb in self.sorted_dos:
			self.pdos_tree.insert("", 'end', iid=orb[2], text=orb[1], values=(orb[0], orb[1], orb[3]))
		corrected_path=expanduser(self.plotpath.get().strip(" "))
		self.color_tree.delete(*self.color_tree.get_children())
		self.QE_files = glob.glob(corrected_path + "/*bands*")
		file_key = {'gnu': 0, 'out': 1, '.in': 2}
		self.QE_sort=[file_key.get(key[-3:], 'Null') for key in self.QE_files]
		if self.QE_sort:
			self.style.configure("Plot.TButton", foreground='green')
			self.style.configure("Save.TButton", foreground='green')
			self.plot_text.set('Create plot')
	def file_select(self, event):
		self.color_tree.selection_toggle(self.color_tree.focus())
	def file_deselect(self):
		for sel in self.color_tree.selection():
			self.color_tree.selection_remove(sel)
	def plot_params(self, *params):
		try:
			scf=glob.glob(expanduser(self.plotpath.get().strip(" ")) + "/*scf*")
			with open(scf[0], 'r+') as scf:
				fermi_level=[float(line.split()[4]) for line in scf if "Fermi" in line][0]
		except:
			fermi_level=False
			self.plot_fermi.set(0)
			self.plot_bg.set(0)
			self.plot_shift.set(0)
		def set_plot(mode, text, relief):
			self.bandplot.draw()
			self.style.configure("Plot.TButton", relief=relief)
			self.plot_text.set(text)
			self.plt_status=mode
		def color_check(key, params, check_var, plot_var):
			if key in params and check_var.get():
				plot_var.configure(fg="green")
			if key in params and not check_var.get():
				plot_var.configure(fg="black")
		color_check(1, params, self.ticks_var, self.pdos_ticks)
		color_check(2, params, self.plot_fermi, self.fermi)
		color_check(3, params, self.plot_bg, self.bandgap)
		color_check(4, params, self.plot_kpath, self.kpath)
		color_check(5, params, self.fill_pdos, self.filldos)
		color_check(8, params, self.plot_shift, self.fermi_shift)
		if 7 in params:
			figure=plt.figure(figsize=(8, 6), dpi=1200)
			self.save=True
		else:
			figure=self.figure
		if 6 in params or 7 in params:
			if len(self.QE_files)==3:
				if 7 not in params:
					self.figure.clf()
				if not self.plt_status and 7 not in params:
					plot_mods.clear_plot(self.figure)
					set_plot(True, "Create plot", "raised")
					return None
				else:
					try:
						y_range=[float(self.eV_a.get()), float(self.eV_b.get())]
					except ValueError:
						return None
					if self.plot_var.get()==2:
						self.bands(y_range, fermi_level, figure, self.save)
						set_plot(False, "Clear plot", "sunken")
						return None
					try:
						x_range=[float(self.dos_a.get()), float(self.dos_b.get())]
					except ValueError:
						return None
					if self.plot_var.get() in [3, 4] and not self.pdos_tree.selection():
						print('At least one atomic orbital must be selected to plot PDOS')
						return None
					if self.plot_var.get()==3:
						self.pdos(x_range, y_range, fermi_level, figure, self.save)
					if self.plot_var.get()==4:
						self.combined(x_range, y_range, fermi_level, figure, self.save)
				if 7 not in params:
					set_plot(False, "Clear plot", "sunken")
	def validate_path(self):
		if not self.check_var.get():
			self.style.configure("Toolbutton", relief="raised", background="lavender")
			self.style.configure("Path.TEntry", fieldbackground="gray80")
			self.style.configure("Make.TButton", background="pink")
			self.path.delete(0, '')
			if self.xpa.get()=='':
				self.xpa.insert(0,'auto')
		else:
			self.style.configure("Toolbutton", relief="sunken", background="seagreen2")
			self.style.configure("Path.TEntry", fieldbackground="ghost white")
			if self.path.get()=='':
				self.path.insert(0,'~/')
			self.style.configure("Make.TButton", background="gray88")
			self.xpa.delete(0,'')

		if not self.inv_var.get():
			self.style.configure("Invert.Toolbutton", relief="raised", background="lavender")
		else:
			self.style.configure("Invert.Toolbutton", relief="sunken", background="seagreen2")
#--------------------------------------------------------------------------------------------------
	def from_file(self, path: str=None):
		expanded_path = expanduser(path)
		cell_from_file=[[path, 0], []]
		try:
			with open(expanded_path, 'r+') as file:
				natom=[int(line.split()[2].replace(",","")) for line in file if "nat" in line][0]
				file.seek(0)
				positions=[file.readline().split() for line in file for i in range(natom) if "ATOMIC_POSITIONS" in line]
				for pos in positions:
					cell_from_file[1].append([pos[0], *[float(p) for p in pos[1:4]]])
		except IndexError:
			print('Error parsing file')
			return None
		return cell_from_file
#--------------------------------------------------------------------------------------------------
		#Add option to create VASP/QE/xsf files?
	#<<NotebookTabChanged>>
	def adsorbate(self, cell, prim_axis, ads, angs, xpa, sort_keys):
		prim_axis=self.key_map.get(self.adsorb_dim.get()[0], 'Null')
		if prim_axis=='Null':
			self.style.configure("Dim.TMenubutton", foreground="red2")
			return None
		if not self.check_var.get():
			xpa=xpa[prim_axis-1]
			self.xpa.delete(0, 'end')
			self.xpa.insert(0, round(xpa,4))
		self.cell=[list(atom) for atom in cell]
		def add_to_cell(prim_axis, adsorb, xpa, angs):
			def ad_index():
				return max([self.tree.item(sel)['values'][0] for sel in self.tree.get_children()])+1
			for sel in self.tree.selection():
					atom=self.tree.item(sel)['values'][1:]
					ad=[adsorb, float(atom[1]), float(atom[2]), float(atom[3])]
					ad[prim_axis]+=copysign(round(xpa*angs, 3), ad[prim_axis])
					ad[prim_axis]=round(ad[prim_axis], 5)
					if ad[prim_axis]<0:
						slot='0'
					else:
						slot='end'
					self.tree.item(sel, tags="red")
					self.tree.tag_configure("red", foreground=('red2'))
					if self.inv_var.get():
						for inv_sel in self.tree.get_children():
							inv_val = self.tree.item(inv_sel)['values']
							if (-float(inv_val[1:5][sort_keys[0]])==float(atom[sort_keys[0]]) and (-float(inv_val[1:5][sort_keys[1]])==float(atom[sort_keys[1]]))) or (-float(inv_val[1:5][sort_keys[0]])==float(atom[sort_keys[0]]) and (-float(inv_val[1:5][sort_keys[2]])==float(atom[sort_keys[2]]))):
								invert=True
								inv_val=[adsorb, *[float(v) for v in inv_val[2:5]]]
								inv_val[prim_axis]+=copysign(round(xpa*angs, 3), inv_val[prim_axis])
								inv_val[prim_axis]=round(inv_val[prim_axis], 5)
								self.tree.item(inv_sel, tags="red")
								self.tree.tag_configure("red", foreground=('red2'))
								if inv_val not in self.cell and inv_val not in self.retain and invert==True:
									self.retain.append(inv_val)
									self.cell.append(inv_val)
									if inv_val[prim_axis]<0:
										inv_slot='0'
									else:
										inv_slot='end'
									self.tree.insert("", inv_slot, iid=ad_index(), text=inv_val, values=(ad_index(), *inv_val), tags=('blue'))
									self.tree.tag_configure("blue", foreground=('blue'))
					self.tree.selection_remove(sel)
					if ad not in self.cell and ad not in self.retain:
						self.retain.append(ad)
						self.cell.append(ad)
						self.tree.insert("", slot, iid=ad_index(), text=ad, values=(ad_index(), *ad), tags=('blue'))
						self.tree.tag_configure("blue", foreground=('blue'))
		add_to_cell(prim_axis, ads, xpa, angs)
#--------------------------------------------------------------------------------------------------
	def validate_params(self, button_id, *args):
		def set_dim(dim, sort, shift):
			try:
				dim = [int(dim.get()) for dim in dim]
				sort_by = [sort.get() for sort in sort]
				shifts= [float(shift.get()) for shift in shift]
				sort_keys = [self.key_map.get(key, 'Null') for key in sort_by]
				output_format=self.format_menu.get()
			except ValueError:
				return None
			return output_format, dim, sort_keys, shifts
		params = set_dim(self.dim, self.sort, self.shifts)
		if self.check_var.get():
			unit_cell=self.from_file(self.path.get())
			source=False
			bilbao_link=False
		else:
			structure=self.struct_menu.get()
			source = True
			if structure=='Structure':
				print('Select a structure or input file')
				return None

			spg = re.split(' |,', structure)[1]
			unparsed_cell=(re.split(' |,', self.yaml_data[structure]))
			while '' in unparsed_cell:
				unparsed_cell.remove('')

			if True in [s.isalpha() for s in unparsed_cell[1]]:
				parser=[5, 'wyckoff']
				bilbao_link=AnsiiCodes.bold+"https://www.cryst.ehu.es/cgi-bin/"\
				"cryst/programs/nph-wp-list?gnum="+str(spg)+'\n'+AnsiiCodes.end
			else:
				parser=[4, 'cell']
				bilbao_link=False
			split_cell = [
				unparsed_cell[i:i + parser[0]] for i in range(0, len(unparsed_cell), parser[0])]
			unit_cell = [
			[structure.split()[0], int(structure.split()[1])],
			[[*data[:-3], *(float(v) for v in data[-3::])] for data in split_cell]]
			unit_cell[1]=CheckPrimitive(unit_cell).constructPrim(parser[1])
		if button_id==1 and len(args)==0:
			SuperCell(source, bilbao_link, unit_cell, *params, self.charges).displayCell()
		if button_id==1 and len(args)==1:
			SuperCell(source, bilbao_link, self.cell, *params, self.charges).displayCell(self.cell)
			SuperCell(source, bilbao_link, self.cell, *params, self.charges).cleaveCell(False, self.cell)
		if button_id==4:# Construct cell for adsorbate menu
			self.retain.clear()
			self.cell.clear()
			self.tree.delete(*self.tree.get_children())
			self.xpa.delete(0, '')
			self.xpa.insert(0,'auto')
			self.insert_admenu(SuperCell(source, bilbao_link, unit_cell, *params, self.charges).layerCell())
		if button_id==2:#cleave cell
			if self.cleave_a.get()==self.cleave_b.get()=='':
				surface=False
			else:
				try:
					surface = [int(self.cleave_a.get()), int(self.cleave_b.get())]
				except ValueError:
					return None
				if surface[0] > surface[1]:
					return None
			SuperCell(source, bilbao_link, unit_cell, *params, self.charges).cleaveCell(surface)
		if button_id==5: # Add to cell
		# Tab for adding unit cells? (editing yml file)
			if self.check_var.get():
				try:
					xpa=float(self.xpa.get())
				except ValueError:
					return None
			if self.enter_adsorb.get().isalpha():
				try: 
					angs = float(self.angs.get())
					ads = self.enter_adsorb.get().upper()
					celldm=[float(self.alat.get()), float(self.blat.get()), float(self.clat.get())]
				except ValueError:
					return None
			prim_axis=self.key_map.get(self.adsorb_dim.get()[0], 'Null')
			if prim_axis=='Null':
				self.style.configure("Dim.TMenubutton", foreground="red2")
				return None
			cell = SuperCell(source, bilbao_link, unit_cell, *params, self.charges).layerCell()
			prim_axis=self.key_map.get(self.adsorb_dim.get()[0], 'Null')
			def xpa_xyz(cell,i, dim):
				xsort=sorted(cell, key = lambda x:(x[i]))
				return (abs(xsort[0][i])+xsort[-1][i])/(float(celldm[i-1])*dim[i-1])
			if not self.check_var.get():
				xpa=[xpa_xyz(cell,i,params[1]) for i in range(1,4)]
			if self.retain==[]:
				self.adsorbate(SuperCell(source, bilbao_link, unit_cell, 
					*params, self.charges).layerCell(), prim_axis, ads, angs, xpa, params[2])
			else:
				self.adsorbate(self.cell, prim_axis, ads, angs, xpa, params[2])
		self.N.set('N atoms: '+str(len(self.cell)))
	def bands(self, y_range, fermi_level, figure, save_key):
		plot_mods.plotbands(*[self.QE_files[self.QE_sort.index(i)] for i in range(0,3)], self.plot_bg.get(),
			fermi_level, self.plot_fermi.get(), self.plot_shift.get(), *y_range, figure, save_key, self.plot_title_entry.get().strip(' '), self.lw_bands.get(), True, False)
		self.bandplot.draw()
	def pdos(self, x_range, y_range, fermi_level, figure, save_key): 
		dos_files=[[self.dos_dict[iid][3], self.dos_dict[iid][1], self.dos_dict[iid][2], self.dos_dict[iid][4]] for iid in self.pdos_tree.selection()]
		plot_mods.dosplot(dos_files, figure, save_key, self.plot_title_entry.get().strip(' '), fermi_level, y_range, x_range, 
			self.dos_dict, self.plot_fermi.get(), self.plot_shift.get(), True, False, self.fill_pdos.get(), self.lw_dos.get())
	def combined(self, x_range, y_range, fermi_level, figure, save_key):
		dos_files=[[self.dos_dict[iid][3], self.dos_dict[iid][1], self.dos_dict[iid][2], self.dos_dict[iid][4]] for iid in self.pdos_tree.selection()]
		plot_mods.bands_and_pdos(*[self.QE_files[self.QE_sort.index(i)] for i in range(0,3)], fermi_level, *y_range, 
			figure, save_key, self.plot_title_entry.get().strip(' '), self.lw_bands.get(), self.lw_dos.get(), self.plot_bg.get(), self.plot_fermi.get(), 
			self.plot_shift.get(), dos_files, x_range, self.dos_dict, self.fill_pdos.get())

#--------------------------------------------------------------------------------------------------
class AnsiiCodes:
    cyan = '\033[96m'
    bold = '\033[1m'
    end = '\033[0m'
    highlight= '\033[01;97;105m'
    pink = '\033[95m'
    green = '\033[32m'
    red = '\033[31m'
class StringFormats:
    main = '{:^4}:  {:2}  | {:12.9f} | {:12.9f} | {:12.9f}'
    abinit = '  {:12.9f}  {:12.9f}  {:12.9f}  #{:3}'
    quantum = '{:2}     {:12.9f}  {:12.9f}  {:12.9f}'
    charge = '\nNet charge of surface: {}{}{}{}\n'
    invert = AnsiiCodes.highlight + '{:^4}:  {:2}  | {:12.9f} | {:12.9f} | {:12.9f}' + AnsiiCodes.end
    invert_notice = AnsiiCodes.bold + '\n***Possible inversion center found at {}-{}***\n' + AnsiiCodes.end
    header = '             {:^12}   {:^12}   {:^12}\n'
    cleaved_header ='      {:^10}    {:^10}    {:^10}\n'
    supercell_name = AnsiiCodes.cyan+AnsiiCodes.bold+'\n{} x {} x {} {} #{} supercell'+AnsiiCodes.end
    cell_name= AnsiiCodes.bold+ '\n{} #{} cell' +AnsiiCodes.end
#--------------------------------------------------------------------------------------------------
app=SuperGui()
app.mainloop()