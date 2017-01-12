# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:19:55 2016

Concatenate files for velocity data

@author: hkarimi
"""

import numpy as np
# Generate file names to concatenate all b
b_array = np.linspace(0,1,21)

uxsdxfiles = []; uxsdyfiles = []; uxsdzfiles = [];
uysdxfiles = []; uysdyfiles = []; uysdzfiles = [];
uzsdxfiles = []; uzsdyfiles = []; uzsdzfiles = [];
for b in b_array:
    uxsdxfiles.append('uxb'+'{:.2f}'.format(b) + 'sdx.txt')
    uxsdyfiles.append('uxb'+'{:.2f}'.format(b) + 'sdy.txt')
    uxsdzfiles.append('uxb'+'{:.2f}'.format(b) + 'sdz.txt')
    uysdxfiles.append('uyb'+'{:.2f}'.format(b) + 'sdx.txt')
    uysdyfiles.append('uyb'+'{:.2f}'.format(b) + 'sdy.txt')
    uysdzfiles.append('uyb'+'{:.2f}'.format(b) + 'sdz.txt')
    uzsdxfiles.append('uzb'+'{:.2f}'.format(b) + 'sdx.txt')
    uzsdyfiles.append('uzb'+'{:.2f}'.format(b) + 'sdy.txt')
    uzsdzfiles.append('uzb'+'{:.2f}'.format(b) + 'sdz.txt')

# concatenate all b
for sd in ['x','y','z']:
    for vd in ['x','y','z']:
        filenames = eval('u' + vd + 'sd' + sd + 'files')
        output = 'u' + vd + 'sd' + sd + '.txt'
        with open('./concatenated/' + output, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())
                    
# concatenate all sd
uxfiles = ['uxsdx.txt','uxsdy.txt','uxsdz.txt']
uyfiles = ['uysdx.txt','uysdy.txt','uysdz.txt']
uzfiles = ['uzsdx.txt','uzsdy.txt','uzsdz.txt']

for vd in ['x','y','z']:
    filenames = eval('u' + vd + 'files')
    output = 'u' + vd + '.txt'
    with open('./concatenated/'+output,'w') as outfile:
        for fname in filenames:
            with open('./concatenated/'+fname) as infile:
                outfile.write(infile.read())
                
numlinesx = sum(1 for line in open('./concatenated/ux.txt'))
numlinesy = sum(1 for line in open('./concatenated/uy.txt'))
numlinesz = sum(1 for line in open('./concatenated/uz.txt'))