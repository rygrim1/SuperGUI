from collections import namedtuple
import sys, re, glob, yaml
from os.path import expanduser 
import supermods.Wyckoff as wyck

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
	def __init__(self, link, unit_cell, output, cell_dm, sort_keys, shifts, charges):
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
		bulk_cell = self.constructCell(*self.cell_dm, *self.shifts)
		layered_cell = sorted(bulk_cell,
			key = lambda x: tuple(x[self.sort_by[i]] for i in range(3)))
		for idx, pos in enumerate(layered_cell):
			if idx < len(layered_cell)-1 and ((list(layered_cell[idx+1])[self.sort_by[0]] == 1-pos[self.sort_by[0]] or list(layered_cell[idx+1]) == [pos[0], -pos[1], -pos[2], -pos[3]])):
					self.invert = True
					self.invert_key = idx+1
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
	def displayCell(self):
		layered_supercell = self.layerCell()
		print(StringFormats.header.format('','','x','y','z'))
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

		while True:
			surface = input(StringFormats.cleave_surface).split('-')
			print('\n')
			if surface != [''] and (len(surface) == 2): 
				cleaved_surface = layered_supercell[int(surface[0])-1:int(surface[1])]
				quantum_sort=sorted(cleaved_surface,
					key = lambda x: tuple(x[self.sort_by[i]] for i in range(3)))
				abi_sort=sorted(cleaved_surface,
					key = lambda x:(x[0],x[self.sort_by[0]]))
			else: 
				print('Exiting...\n')
				raise SystemExit
			for i in range(len(cleaved_surface)):
				if self.format == 'a':
					print(StringFormats.abinit.format(*abi_sort[i][1:4], abi_sort[i][0]))
				if self.format == 'q':
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
				print('\nSurface charge not available. Update atomic charges in cells.yaml')
				pass
			print('Press enter to exit')
			self.cleaved = {}
#--------------------------------------------------------------------------------------------------
class Yaml_Retrieve:
	def __init__(self):
		self.yaml_path=path=str(sys.argv[1])
		with open(self.yaml_path) as yaml_file:
			self.yaml_data = yaml.load(yaml_file, Loader = yaml.FullLoader)
		self.charges = self.yaml_data['Charges']
		self.struct = [key for key in self.yaml_data if key.lower() != 'charges']
		self.key_map = {'x': 1, 'y': 2, 'z': 3 }
	def validate_params(self):
		def set_dim():
			try:
				dim = [int(v) for v in self.yaml_data['Cell dimensions'].values()]
				sort_by = self.yaml_data['Sort priority'].strip(" ").split()
				shifts = [float(v) for v in self.yaml_data['Origin shifts'].values()]
				sort_keys = [self.key_map.get(key, 'Null') for key in sort_by]
				output_format=self.yaml_data['Output format'][0].lower()
				struct = self.yaml_data['Structure']
			except ValueError:
				print('Parameter error. Check cells.yaml')
				return None
			return output_format, dim, sort_keys, shifts
		struct=self.yaml_data['Structure']
		fmtd_struct=struct.strip(" ").split()
		name = fmtd_struct
		struct_a = self.yaml_data['Structures'][" ".join(fmtd_struct)]
		params = set_dim()
		spg = int(name[1])
		unparsed_cell=(re.split(' |,', struct_a))
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
		[name[0], int(name[1])],
		[[*data[:-3], *(float(v) for v in data[-3::])] for data in split_cell]]
		unit_cell[1]=CheckPrimitive(unit_cell).constructPrim(parser[1])
		SuperCell(bilbao_link, unit_cell, *params, self.charges).displayCell()
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
    cleave_surface = '\nEnter two numbers e.g. "10-20" corresponding to the atomic'\
    				' positions you would like to cleave a surface between:\n\n'

if __name__ == '__main__':
	Yaml_Retrieve().validate_params()
