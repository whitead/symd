#!/usr/bin/python
import numpy as np

from scipy.integrate import quad
from math import pi

def read_run_params(fname):
    params = {}
    with open(fname, 'r') as f:
        for line in f.readlines():
            sline = line.split()
            if(len(sline) == 2):
                params[sline[0]] = sline[1]
    return params

def pair_distance(x, y, ndims, img, PBC):
    dsum = 0
    dist = 0
    for j in range(ndims):
        dist = (x[j] - y[j])
        if(PBC):
            dist = dist - round(dist / img[j]) * img[j]
        dsum += dist ** 2
    return np.sqrt(dsum)

def bin_atom(i, distances, positions, bin_edges, ndims, img, PBC):
    for j in range(i+1, len(positions)):
        distances[j] = pair_distance(positions[i], positions[j], ndims, img, PBC)
    return np.histogram(distances[i+1:], bin_edges)

#square regions
def region_1(rho, d):
    return 0.5 * (d**2 * np.sqrt(rho**2 - 2 * d ** 2))

def region_1_2_theta(rho, d):
    return np.arccos(d / np.sqrt(rho ** 2 - d ** 2))

#curved square regions
def region_2_integrand(theta, rho, d):
    return np.sqrt(np.cos(theta) ** 2 - (d / rho)**2) / (np.cos(theta) ** 3)

def region_2(rho, d):
    i4 = d**3 / 6. * (rho**2 / d **2 - 1) * (pi / 4 - region_1_2_theta(rho, d))
    i3 = d ** 2 * rho / 3. * quad(region_2_integrand, region_1_2_theta(rho,d), pi / 4, args=(rho, d))[0] 
    return i4 + i3


#spherical region
def region_3_integrand(theta, rho, d):
    return np.sqrt(np.cos(theta) ** 2 - (d / rho)**2) / np.cos(theta)

def region_3(rho, d):
    return rho ** 3 / 3. * (d / rho * (pi / 4 - region_1_2_theta(rho, d)) - quad(region_3_integrand, region_1_2_theta(rho, d), pi / 4, args=(rho, d))[0])

def calc_volume(rho, d):

    alpha = rho / d

    if(alpha <= 1):
        return 4./3 * pi * (rho) ** 3

    if(alpha <= np.sqrt(2)):
        return 4. / 3 * pi * (rho) ** 3 - 6. * pi * (2 * (rho)**3 - 3 * d * (rho)**2 + d**3) / 3.

    if(alpha < np.sqrt(3)):
        return 16. * (region_1(rho,d) + region_2(rho, d) + region_3(rho,d))

    return 8. * d ** 3




def calc_gofr(xyz_file, box_size, bin_width=0.1, ndims=3, PBC=True): 

    if(ndims != 3):
        raise Exception("Have not implemented n-dimensional g(r)")
    #if we don't have a cube, I don't yet have the g(r) cube-sphere intersection implemented
    if(not(box_size[0] == box_size[1] and box_size[1] == box_size[2])):
        bin_edges = np.arange(0, min(box_size) / 2, bin_width)
    else:
        bin_edges = np.arange(0, np.sqrt(3) * box_size[0] / 2, bin_width)
    (hist, bin_edges) = np.histogram([-1], bin_edges)    

    with open(xyz_file, 'r') as f:

        lines_read = 2
        frames_read = 0
        natoms = 0
        try:
            natoms = int(f.readline())
        except ValueError:
            raise Exception("Doesn't look like %s in an xyz file" % xyz_file)

        if(natoms == 0):
            return None

        positions = [0 for x in range(natoms)]
        distances = [0 for x in range(natoms)]
        
        print "Processing frame...0",
        while(1):
            try:
                #skip over frames and number of atoms
                line = f.readline() 
                while(line != "" and line.split(":")[0] != "Frame"):
                    line = f.readline() 
                for i in range(natoms):
                    tokens = f.readline().split()
                    if(len(tokens) != ndims + 1):
                        raise EOFError()
                    positions[i] = [float(x) for x in tokens[1:]]
                for i in range(natoms):
                    (newhist, bin_edges) = bin_atom(i,distances, positions, bin_edges, ndims, box_size, PBC)
                    hist = np.add(hist, newhist)
                frames_read += 1
                print "\rProcessing frame...%d" % frames_read,
            except EOFError:
                break

        print "\nCalculating g(r)..."
        volume = 1.
        for x in box_size:
            volume = volume * x
        bulk_density = sum(hist) / volume 

    r = np.empty(len(bin_edges) - 1)
    n = np.empty(len(bin_edges) - 1)
    gofr = np.empty(len(bin_edges) - 1)
    
    for i in range(len(bin_edges) - 1) :
        v = calc_volume(bin_edges[i+1], box_size[0] / 2.) - calc_volume(bin_edges[i], box_size[0] / 2.)
        ideal = v * bulk_density
        actual = hist[i]
        r[i] = (bin_edges[i + 1]  + bin_edges[i]) / 2.
        n[i] = actual
        gofr[i] = actual / ideal
    return (r,n,gofr)

if __name__ == "__main__":
    import sys
    print "Note: This program assumes your box is a cube. It does not check"

    if(len(sys.argv) != 5):
        print "usage: [gofr.py] [input xyz file] [output file] [run_params_file] [bin_width]"
        exit()
    
    bin_width = float(sys.argv[4])
    params = read_run_params(sys.argv[3])
    ndims = int(params['n_dims'])
    (r, n, gofr) = calc_gofr(sys.argv[1],
                   box_size=[int(params['box_%d_size' % x]) for x in range(1,ndims+1)], 
                   bin_width=bin_width,
                   ndims=ndims)
    
    with open(sys.argv[2], 'w') as f:
        for (ri, ni, gofri) in zip(r,n,gofr):
            f.write("%g %g %g\n" % (ri, ni, gofri))

