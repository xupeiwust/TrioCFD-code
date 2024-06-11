import os
import numpy as np

import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", help="Input lata file file")
parser.add_argument("output", help="Output column text file")
parser.add_argument("--static", action="store_true", help="Indicates a case without bubble motion")
args = parser.parse_args()

f = np.fromfile(args.input, dtype=np.float32)

length = 0.005
n = int(np.round(np.shape(f)[0]**(1./3.)))

x = np.arange(n)*length/n - 0.002

m = np.moveaxis(np.array(np.meshgrid(x,x,x,indexing='ij')), 0, -1)
f = np.reshape(f, np.shape(m)[:3], order='F')

final_y_location = 0
final_z_location = 0 if args.static else .5e-3
jslice = np.argwhere(m[0,:,0,1] >= final_y_location).flatten()[0]
kslice = np.argwhere(m[0,0,:,2] >= final_z_location).flatten()[0]

x_line = m[:,jslice,kslice,0]
f_line = f[:,jslice,kslice]

np.savetxt(args.output, np.vstack((x_line, f_line)).T)
