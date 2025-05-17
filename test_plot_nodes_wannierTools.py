from plotting_class import plottingTools
from styling_class import styler

import matplotlib.pyplot as plt

import time

path_read1 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_71'
path_read2 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_72'
path_read3 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_73'
path_read4 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_74'
path_read5 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_75'
path_read6 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_76'
path_read7 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_77'
path_read8 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_78'
path_read9 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_79'
path_read10 = '/home/carlos/Documents/scriptsCalculations/projects/RPt2B/LaPt2B/PBEsol/SG180/k997/200bands/150MLWFs/nodes_wannierTools/find_nodes_80'
fig, axs = plt.subplots(4, 1, figsize=[12, 10]) # 15, 6 # 5, 7

t0 = time.time()
plotter = plottingTools(main_linewidth=1.8, k_linewidth=0.8)
# =========================================================
# =========================================================
kpath = [[0,   0,   0],  # G
         [0,   1/2, 0],  # M
         [1/3, 1/3, 0],  # K
         [0,   0,   0],  # G

         [0,   0,   1/2],  # A
         [1/2, 0,   1/2],  # L
         [1/3, 1/3, 1/2],  # H
         [0,   0,   1/2],  # A

         [1/2, 0,   1/2],  # L
         [1/2, 0,   0  ],  # M
         
         [1/3, 1/3, 1/2],  # H
         [1/3, 1/3, 0  ],  # K
         ]
kpath = None

nodes, E, gap, path = plotter.separate_nodes(path_read1, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read1+'/nodes_path_separated_k_E.txt')
nodes, E, gap, path = plotter.separate_nodes(path_read2, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read2+'/nodes_path_separated_k_E.txt')
nodes, E, gap, path = plotter.separate_nodes(path_read3, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read3+'/nodes_path_separated_k_E.txt')
nodes, E, gap, path = plotter.separate_nodes(path_read4, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read4+'/nodes_path_separated_k_E.txt')
nodes, E, gap, path = plotter.separate_nodes(path_read5, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read5+'/nodes_path_separated_k_E.txt')
nodes, E, gap, path = plotter.separate_nodes(path_read6, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read6+'/nodes_path_separated_k_E.txt')
nodes, E, gap, path = plotter.separate_nodes(path_read7, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read7+'/nodes_path_separated_k_E.txt')
nodes, E, gap, path = plotter.separate_nodes(path_read8, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read8+'/nodes_path_separated_k_E.txt')
# nodes, E, gap, path = plotter.separate_nodes(path_read9, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read9+'/nodes_path_separated_k_E.txt')
# nodes, E, gap, path = plotter.separate_nodes(path_read10, kpath, [-0.1, 0.1], ktol=0.0001, savefile=path_read10+'/nodes_path_separated_k_E.txt')