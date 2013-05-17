from SimpleMD import SimpleMD
from pylab import *

md = SimpleMD(100, 3, steps=50000, exeDir="grtest")
md.setup_positions(density=0.7)
md.setup_masses(1)
md.log_positions(period=5)
md.log_output()
md.execute()
(r, n, gr) =  md.calc_gofr(bin_width=0.025)


plot(r, gr, linewidth=1.5)

savefig('gr.png', bbox_inches=0)

md.clean_files()
