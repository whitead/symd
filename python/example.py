from SimpleMD import SimpleMD
import numpy as np

md = SimpleMD(100, 3, steps=500, exeDir="example")
md.setup_positions(density=0.7)
md.log_positions(period=5)
md.log_output()
md.execute()
print(f'Average Potential Energy: {np.mean(md.pe):.3f}')
