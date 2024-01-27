# Generic Simulation Script

import numpy as np
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
# from rod_derivative_bank import *
# from triangle_derivative_bank_h3 import *
# from triangle_derivative_bank_s3 import *
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, h3distb2, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, r4distb2, killingvecs3, rot2transh3pos, rot2transh3vel, boostxh3, rotzh3
# from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3, h3rads2tridae, h3rads3tridae, s3rads2tridae, s3rads3tridae


# Time arrays based on time step
t_max = 10      # Total simulation time

dt1 = .1    
t_arr1 = np.arange(0.,t_max+dt1,dt1)

dt01 = .01    
t_arr01 = np.arange(0.,t_max+dt01,dt01)

dt001 = .001    
t_arr001 = np.arange(0.,t_max+dt001,dt001)

dt0001 = .0001    
t_arr0001 = np.arange(0.,t_max+dt0001,dt0001)

dt00001 = .00001    
t_arr00001 = np.arange(0.,t_max+dt00001,dt00001)

# # # Rod Spring Data

#     # H3

# # mathematicadata = np.genfromtxt("./data_files/spring_potential/rod/h3_gauss3s_h0001_t10_boo.dat")

# datarh3gtdt00001gs3  = np.load("./data_files/spring_potential/rod/h3_r_gausss3_anal_tmax10_dt00001.npy")
# # datarh3gtdt000005gs3 = np.load("./data_files/spring_potential/rod/h3_r_gausss3_gt_tmax10_dt000005.npy")

# datarh3dt1gs1     = np.load("./data_files/spring_potential/rod/h3_r_gausss1_sim_tmax10_dt1.npy")
# datarh3dt1gs2     = np.load("./data_files/spring_potential/rod/h3_r_gausss2_sim_tmax10_dt1.npy")
# datarh3dt1gs3     = np.load("./data_files/spring_potential/rod/h3_r_gausss3_sim_tmax10_dt1.npy")

# # datarh3dt1gs1m    = np.genfromtxt("./data_files/mathematica_data/gauss_error_analysis_data/h3_gauss1s_h1_t10.dat")
# # datarh3dt1gs2m    = np.genfromtxt("./data_files/mathematica_data/gauss_error_analysis_data/h3_gauss2s_h1_t10.dat")
# # datarh3dt1gs3m    = np.genfromtxt("./data_files/mathematica_data/gauss_error_analysis_data/h3_gauss3s_h1_t10.dat")

# datarh3dt01gs1    = np.load("./data_files/spring_potential/rod/h3_r_gausss1_sim_tmax10_dt01.npy")
# datarh3dt01gs2    = np.load("./data_files/spring_potential/rod/h3_r_gausss2_sim_tmax10_dt01.npy")
# datarh3dt01gs3    = np.load("./data_files/spring_potential/rod/h3_r_gausss3_sim_tmax10_dt01.npy")

# datarh3dt001gs1   = np.load("./data_files/spring_potential/rod/h3_r_gausss1_sim_tmax10_dt001.npy")
# datarh3dt001gs2   = np.load("./data_files/spring_potential/rod/h3_r_gausss2_sim_tmax10_dt001.npy")
# datarh3dt001gs3   = np.load("./data_files/spring_potential/rod/h3_r_gausss3_sim_tmax10_dt001.npy")

# # datarh3dt001gs1m    = np.genfromtxt("./data_files/mathematica_data/gauss_error_analysis_data/h3_gauss1s_h001_t10.dat")
# # datarh3dt001gs2m    = np.genfromtxt("./data_files/mathematica_data/gauss_error_analysis_data/h3_gauss2s_h001_t10.dat")
# # datarh3dt001gs3m    = np.genfromtxt("./data_files/mathematica_data/gauss_error_analysis_data/h3_gauss3s_h001_t10.dat")

# datarh3dt0001gs1  = np.load("./data_files/spring_potential/rod/h3_r_gausss1_sim_tmax10_dt0001.npy")
# datarh3dt0001gs2  = np.load("./data_files/spring_potential/rod/h3_r_gausss2_sim_tmax10_dt0001.npy")
# datarh3dt0001gs3  = np.load("./data_files/spring_potential/rod/h3_r_gausss3_sim_tmax10_dt0001.npy")

# datarh3dt00001gs1 = np.load("./data_files/spring_potential/rod/h3_r_gausss1_sim_tmax10_dt00001.npy")
# datarh3dt00001gs2 = np.load("./data_files/spring_potential/rod/h3_r_gausss2_sim_tmax10_dt00001.npy")
# datarh3dt00001gs3 = np.load("./data_files/spring_potential/rod/h3_r_gausss3_sim_tmax10_dt00001.npy")

#     # S3

# datars3gtdt00001gs3  = np.load("./data_files/spring_potential/rod/s3_r_gausss3_anal_tmax10_dt00001.npy")
# # datars3gtdt000005gs3 = np.load("./data_files/spring_potential/rod/s3_r_gausss3_gt_tmax10_dt000005.npy")

# datars3dt1gs1     = np.load("./data_files/spring_potential/rod/s3_r_gausss1_sim_tmax10_dt1.npy")
# datars3dt1gs2     = np.load("./data_files/spring_potential/rod/s3_r_gausss2_sim_tmax10_dt1.npy")
# datars3dt1gs3     = np.load("./data_files/spring_potential/rod/s3_r_gausss3_sim_tmax10_dt1.npy")

# datars3dt01gs1    = np.load("./data_files/spring_potential/rod/s3_r_gausss1_sim_tmax10_dt01.npy")
# datars3dt01gs2    = np.load("./data_files/spring_potential/rod/s3_r_gausss2_sim_tmax10_dt01.npy")
# datars3dt01gs3    = np.load("./data_files/spring_potential/rod/s3_r_gausss3_sim_tmax10_dt01.npy")

# datars3dt001gs1   = np.load("./data_files/spring_potential/rod/s3_r_gausss1_sim_tmax10_dt001.npy")
# datars3dt001gs2   = np.load("./data_files/spring_potential/rod/s3_r_gausss2_sim_tmax10_dt001.npy")
# datars3dt001gs3   = np.load("./data_files/spring_potential/rod/s3_r_gausss3_sim_tmax10_dt001.npy")

# datars3dt0001gs1  = np.load("./data_files/spring_potential/rod/s3_r_gausss1_sim_tmax10_dt0001.npy")
# datars3dt0001gs2  = np.load("./data_files/spring_potential/rod/s3_r_gausss2_sim_tmax10_dt0001.npy")
# datars3dt0001gs3  = np.load("./data_files/spring_potential/rod/s3_r_gausss3_sim_tmax10_dt0001.npy")

# datars3dt00001gs1 = np.load("./data_files/spring_potential/rod/s3_r_gausss1_sim_tmax10_dt00001.npy")
# datars3dt00001gs2 = np.load("./data_files/spring_potential/rod/s3_r_gausss2_sim_tmax10_dt00001.npy")
# datars3dt00001gs3 = np.load("./data_files/spring_potential/rod/s3_r_gausss3_sim_tmax10_dt00001.npy")

# Triangle Spring Data

    # H3

datath3gtdt00001gs3  = np.load("./data_files/spring_potential/triangle/h3_t_gausss3_anal_tmax10_dt00001.npy")
# datath3gtdt000005gs3 = np.load("./data_files/spring_potential/triangle/h3_t_gausss3_gt_tmax10_dt000005.npy")

datath3dt1gs1 = np.load("./data_files/spring_potential/triangle/h3_t_gausss1_sim_tmax10_dt1.npy")
datath3dt1gs2 = np.load("./data_files/spring_potential/triangle/h3_t_gausss2_sim_tmax10_dt1.npy")
datath3dt1gs3 = np.load("./data_files/spring_potential/triangle/h3_t_gausss3_sim_tmax10_dt1.npy")

datath3dt01gs1 = np.load("./data_files/spring_potential/triangle/h3_t_gausss1_sim_tmax10_dt01.npy")
datath3dt01gs2 = np.load("./data_files/spring_potential/triangle/h3_t_gausss2_sim_tmax10_dt01.npy")
datath3dt01gs3 = np.load("./data_files/spring_potential/triangle/h3_t_gausss3_sim_tmax10_dt01.npy")

datath3dt001gs1 = np.load("./data_files/spring_potential/triangle/h3_t_gausss1_sim_tmax10_dt001.npy")
datath3dt001gs2 = np.load("./data_files/spring_potential/triangle/h3_t_gausss2_sim_tmax10_dt001.npy")
datath3dt001gs3 = np.load("./data_files/spring_potential/triangle/h3_t_gausss3_sim_tmax10_dt001.npy")

datath3dt0001gs1 = np.load("./data_files/spring_potential/triangle/h3_t_gausss1_sim_tmax10_dt0001.npy")
datath3dt0001gs2 = np.load("./data_files/spring_potential/triangle/h3_t_gausss2_sim_tmax10_dt0001.npy")
datath3dt0001gs3 = np.load("./data_files/spring_potential/triangle/h3_t_gausss3_sim_tmax10_dt0001.npy")

datath3dt00001gs1 = np.load("./data_files/spring_potential/triangle/h3_t_gausss1_sim_tmax10_dt00001.npy")
datath3dt00001gs2 = np.load("./data_files/spring_potential/triangle/h3_t_gausss2_sim_tmax10_dt00001.npy")
datath3dt00001gs3 = np.load("./data_files/spring_potential/triangle/h3_t_gausss3_sim_tmax10_dt00001.npy")

    # S3

datats3gtdt00001gs3  = np.load("./data_files/spring_potential/triangle/s3_t_gausss3_anal_tmax10_dt00001.npy")
# datats3gtdt000005gs3 = np.load("./data_files/spring_potential/triangle/s3_t_gausss3_gt_tmax10_dt000005.npy")

datats3dt1gs1 = np.load("./data_files/spring_potential/triangle/s3_t_gausss1_sim_tmax10_dt1.npy")
datats3dt1gs2 = np.load("./data_files/spring_potential/triangle/s3_t_gausss2_sim_tmax10_dt1.npy")
datats3dt1gs3 = np.load("./data_files/spring_potential/triangle/s3_t_gausss3_sim_tmax10_dt1.npy")

datats3dt01gs1 = np.load("./data_files/spring_potential/triangle/s3_t_gausss1_sim_tmax10_dt01.npy")
datats3dt01gs2 = np.load("./data_files/spring_potential/triangle/s3_t_gausss2_sim_tmax10_dt01.npy")
datats3dt01gs3 = np.load("./data_files/spring_potential/triangle/s3_t_gausss3_sim_tmax10_dt01.npy")

datats3dt001gs1 = np.load("./data_files/spring_potential/triangle/s3_t_gausss1_sim_tmax10_dt001.npy")
datats3dt001gs2 = np.load("./data_files/spring_potential/triangle/s3_t_gausss2_sim_tmax10_dt001.npy")
datats3dt001gs3 = np.load("./data_files/spring_potential/triangle/s3_t_gausss3_sim_tmax10_dt001.npy")

datats3dt0001gs1 = np.load("./data_files/spring_potential/triangle/s3_t_gausss1_sim_tmax10_dt0001.npy")
datats3dt0001gs2 = np.load("./data_files/spring_potential/triangle/s3_t_gausss2_sim_tmax10_dt0001.npy")
datats3dt0001gs3 = np.load("./data_files/spring_potential/triangle/s3_t_gausss3_sim_tmax10_dt0001.npy")

datats3dt00001gs1 = np.load("./data_files/spring_potential/triangle/s3_t_gausss1_sim_tmax10_dt00001.npy")
datats3dt00001gs2 = np.load("./data_files/spring_potential/triangle/s3_t_gausss2_sim_tmax10_dt00001.npy")
datats3dt00001gs3 = np.load("./data_files/spring_potential/triangle/s3_t_gausss3_sim_tmax10_dt00001.npy")

# Rod Rigid Data

#     # H3

# # datarigrh3dt1rs2 = np.load("./data_files/rigid_rod/rod/h3_rig_r_rads2_sim_tmax10_dt1.npy")
# # datarigrh3dt1rs3 = np.load("./data_files/rigid_rod/rod/h3_rig_r_rads3_sim_tmax10_dt1.npy")

# # datarigrh3dt01rs2 = np.load("./data_files/rigid_rod/rod/h3_rig_r_rads2_sim_tmax10_dt01.npy")
# # datarigrh3dt01rs3 = np.load("./data_files/rigid_rod/rod/h3_rig_r_rads3_sim_tmax10_dt01.npy")

# # datarigrh3dt001rs2 = np.load("./data_files/rigid_rod/rod/h3_rig_r_rads2_sim_tmax10_dt001.npy")
# # datarigrh3dt001rs3 = np.load("./data_files/rigid_rod/rod/h3_rig_r_rads3_sim_tmax10_dt001.npy")

# # thesis data

# datarigrh3dt1rs2 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_r_rads2_sim_tmax10_dt1.npy")
# datarigrh3dt1rs3 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_r_rads3_sim_tmax10_dt1.npy")

# datarigrh3dt01rs2 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_r_rads2_sim_tmax10_dt01.npy")
# datarigrh3dt01rs3 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_r_rads3_sim_tmax10_dt01.npy")

# datarigrh3dt001rs2 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_r_rads2_sim_tmax10_dt001.npy")
# datarigrh3dt001rs3 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_r_rads3_sim_tmax10_dt001.npy")

#     # S3

# # datarigrs3dt1rs2 = np.load("./data_files/rigid_rod/rod/s3_rig_r_rads2_sim_tmax10_dt1.npy")
# # datarigrs3dt1rs3 = np.load("./data_files/rigid_rod/rod/s3_rig_r_rads3_sim_tmax10_dt1.npy")

# # datarigrs3dt01rs2 = np.load("./data_files/rigid_rod/rod/s3_rig_r_rads2_sim_tmax10_dt01.npy")
# # datarigrs3dt01rs3 = np.load("./data_files/rigid_rod/rod/s3_rig_r_rads3_sim_tmax10_dt01.npy")

# # datarigrs3dt001rs2 = np.load("./data_files/rigid_rod/rod/s3_rig_r_rads2_sim_tmax10_dt001.npy")
# # datarigrs3dt001rs3 = np.load("./data_files/rigid_rod/rod/s3_rig_r_rads3_sim_tmax10_dt001.npy")

# # thesis data

# datarigrs3dt1rs2 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_r_rads2_sim_tmax10_dt1.npy")
# datarigrs3dt1rs3 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_r_rads3_sim_tmax10_dt1.npy")

# datarigrs3dt01rs2 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_r_rads2_sim_tmax10_dt01.npy")
# datarigrs3dt01rs3 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_r_rads3_sim_tmax10_dt01.npy")

# datarigrs3dt001rs2 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_r_rads2_sim_tmax10_dt001.npy")
# datarigrs3dt001rs3 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_r_rads3_sim_tmax10_dt001.npy")

# # Triangle Rigid Data

# #     # H3

# # datarigth3dt1rs2 = np.load("./data_files/rigid_rod/triangle/h3_rig_t_rads2_sim_tmax10_dt1.npy")
# # datarigth3dt1rs3 = np.load("./data_files/rigid_rod/triangle/h3_rig_t_rads3_sim_tmax10_dt1.npy")

# # datarigth3dt01rs2 = np.load("./data_files/rigid_rod/triangle/h3_rig_t_rads2_sim_tmax10_dt01.npy")
# # datarigth3dt01rs3 = np.load("./data_files/rigid_rod/triangle/h3_rig_t_rads3_sim_tmax10_dt01.npy")

# # datarigth3dt001rs2 = np.load("./data_files/rigid_rod/triangle/h3_rig_t_rads2_sim_tmax10_dt001.npy")
# # datarigth3dt001rs3 = np.load("./data_files/rigid_rod/triangle/h3_rig_t_rads3_sim_tmax10_dt001.npy")

# # thesis data

# datarigth3dt1rs2 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_t_rads2_sim_tmax10_dt1.npy")
# datarigth3dt1rs3 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_t_rads3_sim_tmax10_dt1.npy")

# datarigth3dt01rs2 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_t_rads2_sim_tmax10_dt01.npy")
# datarigth3dt01rs3 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_t_rads3_sim_tmax10_dt01.npy")

# datarigth3dt001rs2 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_t_rads2_sim_tmax10_dt001.npy")
# datarigth3dt001rs3 = np.load("./data_files/rigid_rod/thesis_data/h3_rig_t_rads3_sim_tmax10_dt001.npy")

# #     # S3

# # datarigts3dt1rs2 = np.load("./data_files/rigid_rod/triangle/s3_rig_t_rads2_sim_tmax10_dt1.npy")
# # datarigts3dt1rs3 = np.load("./data_files/rigid_rod/triangle/s3_rig_t_rads3_sim_tmax10_dt1.npy")

# # datarigts3dt01rs2 = np.load("./data_files/rigid_rod/triangle/s3_rig_t_rads2_sim_tmax10_dt01.npy")
# # datarigts3dt01rs3 = np.load("./data_files/rigid_rod/triangle/s3_rig_t_rads3_sim_tmax10_dt01.npy")

# # datarigts3dt001rs2 = np.load("./data_files/rigid_rod/triangle/s3_rig_t_rads2_sim_tmax10_dt001.npy")
# # datarigts3dt001rs3 = np.load("./data_files/rigid_rod/triangle/s3_rig_t_rads3_sim_tmax10_dt001.npy")

# # thesis data

# datarigts3dt1rs2 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_t_rads2_sim_tmax10_dt1.npy")
# datarigts3dt1rs3 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_t_rads3_sim_tmax10_dt1.npy")

# datarigts3dt01rs2 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_t_rads2_sim_tmax10_dt01.npy")
# datarigts3dt01rs3 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_t_rads3_sim_tmax10_dt01.npy")

# datarigts3dt001rs2 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_t_rads2_sim_tmax10_dt001.npy")
# datarigts3dt001rs3 = np.load("./data_files/rigid_rod/thesis_data/s3_rig_t_rads3_sim_tmax10_dt001.npy")


# Rigid Triangle H2xE data

# # Length data

# datarigth2edt01s12 = np.loadtxt("./data_files/forbidden_isometry/h2eriglen3sep1data.dat")
# datarigth2edt01s13 = np.loadtxt("./data_files/forbidden_isometry/h2eriglen3sep2data.dat")
# datarigth2edt01s23 = np.loadtxt("./data_files/forbidden_isometry/h2eriglen3sep3data.dat")

# # Angle data

# datarigth2edt01ang1 = np.loadtxt("./data_files/forbidden_isometry/h2eriglen3ang1data.dat")
# datarigth2edt01ang2 = np.loadtxt("./data_files/forbidden_isometry/h2eriglen3ang2data.dat")
# datarigth2edt01ang3 = np.loadtxt("./data_files/forbidden_isometry/h2eriglen3ang3data.dat")

# Curvature Deteector Data

# databumpcurvdet = np.load("./data_files/curvature_detector/bumpgeo_gausss3_sim_tmax10_dtp01.npy")


#########
# Plots #
#########

# activate latex text rendering
rc('text', usetex=True)
rc('axes', linewidth=2)
rc('font', weight='bold')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']


###################################################

######### Spring Rod plots for H3 and S3 ##########

###################################################


#####
# Spring Potential Rod Separation plot for H3 and S3
#####


# fig,ax=plt.subplots(1,1,figsize=((8,6)))

# h3distdatarspgs3 = np.zeros(np.shape(datarh3dt0001gs3)[0])
# s3distdatarspgs3 = np.zeros(np.shape(datars3dt0001gs3)[0])

# counter = 0
# for a in range(np.shape(h3distdatarspgs3)[0]):
#     h3distdatarspgs3[counter] = h3dist(rot2hyp(datarh3dt0001gs3[a][0:3]),rot2hyp(datarh3dt0001gs3[a][3:6]))
#     s3distdatarspgs3[counter] = r4dist(rot2r4(datars3dt0001gs3[a][0:3]),rot2r4(datars3dt0001gs3[a][3:6]))
#     counter += 1

# ax.plot(t_arr0001,h3distdatarspgs3,'r',label = r'$\mathbf{H}^3$', linewidth=2)
# ax.plot(t_arr0001,s3distdatarspgs3,'k',label = r'$\mathbf{S}^3$', linewidth=2)
# ax.legend(fontsize="20")
# # ax.set_title('Vertex Separation vs. Time')

# plt.title(r'\textbf{Elastic Rod Body Dynamics}', fontsize=20)
# plt.ylabel(r'\textbf{Vertex Separation}', fontsize=20)
# plt.xlabel(r'\textbf{Time}', fontsize=20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)

# # ax.set_xlabel('Time')
# # ax.set_ylabel('Vertex Separation')
# plt.show()


######
# Spring Potential Rod Separation Error as function of distance from origin plot for H3 and S3
######

# ### Gauss s1

# h3disterrdatarspgs1dt1raw     = np.zeros(np.shape(datarh3dt1gs1)[0])
# h3disterrdatarspgs1dt01raw    = np.zeros(np.shape(datarh3dt01gs1)[0])
# h3disterrdatarspgs1dt001raw   = np.zeros(np.shape(datarh3dt001gs1)[0])
# h3disterrdatarspgs1dt0001raw  = np.zeros(np.shape(datarh3dt0001gs1)[0])
# h3disterrdatarspgs1dt00001raw = np.zeros(np.shape(datarh3dt00001gs1)[0])

# h3distgs1dt1     = np.zeros(np.shape(datarh3dt1gs1)[0])
# h3distgs1dt01    = np.zeros(np.shape(datarh3dt01gs1)[0])
# h3distgs1dt001   = np.zeros(np.shape(datarh3dt001gs1)[0])
# h3distgs1dt0001  = np.zeros(np.shape(datarh3dt0001gs1)[0])
# h3distgs1dt00001 = np.zeros(np.shape(datarh3dt00001gs1)[0])

# s3disterrdatarspgs1dt1raw     = np.zeros(np.shape(datars3dt1gs1)[0])
# s3disterrdatarspgs1dt01raw    = np.zeros(np.shape(datars3dt01gs1)[0])
# s3disterrdatarspgs1dt001raw   = np.zeros(np.shape(datars3dt001gs1)[0])
# s3disterrdatarspgs1dt0001raw  = np.zeros(np.shape(datars3dt0001gs1)[0])
# s3disterrdatarspgs1dt00001raw = np.zeros(np.shape(datars3dt00001gs1)[0])

# s3distgs1dt1     = np.zeros(np.shape(datars3dt1gs1)[0])
# s3distgs1dt01    = np.zeros(np.shape(datars3dt01gs1)[0])
# s3distgs1dt001   = np.zeros(np.shape(datars3dt001gs1)[0])
# s3distgs1dt0001  = np.zeros(np.shape(datars3dt0001gs1)[0])
# s3distgs1dt00001 = np.zeros(np.shape(datars3dt00001gs1)[0])

# ### Gauss s2

# h3disterrdatarspgs2dt1raw     = np.zeros(np.shape(datarh3dt1gs2)[0])
# h3disterrdatarspgs2dt01raw    = np.zeros(np.shape(datarh3dt01gs2)[0])
# h3disterrdatarspgs2dt001raw   = np.zeros(np.shape(datarh3dt001gs2)[0])
# h3disterrdatarspgs2dt0001raw  = np.zeros(np.shape(datarh3dt0001gs2)[0])
# h3disterrdatarspgs2dt00001raw = np.zeros(np.shape(datarh3dt00001gs2)[0])

# h3distgs2dt1     = np.zeros(np.shape(datarh3dt1gs2)[0])
# h3distgs2dt01    = np.zeros(np.shape(datarh3dt01gs2)[0])
# h3distgs2dt001   = np.zeros(np.shape(datarh3dt001gs2)[0])
# h3distgs2dt0001  = np.zeros(np.shape(datarh3dt0001gs2)[0])
# h3distgs2dt00001 = np.zeros(np.shape(datarh3dt00001gs2)[0])

# s3disterrdatarspgs2dt1raw     = np.zeros(np.shape(datars3dt1gs2)[0])
# s3disterrdatarspgs2dt01raw    = np.zeros(np.shape(datars3dt01gs2)[0])
# s3disterrdatarspgs2dt001raw   = np.zeros(np.shape(datars3dt001gs2)[0])
# s3disterrdatarspgs2dt0001raw  = np.zeros(np.shape(datars3dt0001gs2)[0])
# s3disterrdatarspgs2dt00001raw = np.zeros(np.shape(datars3dt00001gs2)[0])

# s3distgs2dt1     = np.zeros(np.shape(datars3dt1gs2)[0])
# s3distgs2dt01    = np.zeros(np.shape(datars3dt01gs2)[0])
# s3distgs2dt001   = np.zeros(np.shape(datars3dt001gs2)[0])
# s3distgs2dt0001  = np.zeros(np.shape(datars3dt0001gs2)[0])
# s3distgs2dt00001 = np.zeros(np.shape(datars3dt00001gs2)[0])

# ### Gauss s3

# h3disterrdatarspgs3dt1raw     = np.zeros(np.shape(datarh3dt1gs3)[0])
# h3disterrdatarspgs3dt01raw    = np.zeros(np.shape(datarh3dt01gs3)[0])
# h3disterrdatarspgs3dt001raw   = np.zeros(np.shape(datarh3dt001gs3)[0])
# h3disterrdatarspgs3dt0001raw  = np.zeros(np.shape(datarh3dt0001gs3)[0])
# h3disterrdatarspgs3dt00001raw = np.zeros(np.shape(datarh3dt00001gs3)[0])

# h3distgs3dt1     = np.zeros(np.shape(datarh3dt1gs3)[0])
# h3distgs3dt01    = np.zeros(np.shape(datarh3dt01gs3)[0])
# h3distgs3dt001   = np.zeros(np.shape(datarh3dt001gs3)[0])
# h3distgs3dt0001  = np.zeros(np.shape(datarh3dt0001gs3)[0])
# h3distgs3dt00001 = np.zeros(np.shape(datarh3dt00001gs3)[0])

# s3disterrdatarspgs3dt1raw     = np.zeros(np.shape(datars3dt1gs3)[0])
# s3disterrdatarspgs3dt01raw    = np.zeros(np.shape(datars3dt01gs3)[0])
# s3disterrdatarspgs3dt001raw   = np.zeros(np.shape(datars3dt001gs3)[0])
# s3disterrdatarspgs3dt0001raw  = np.zeros(np.shape(datars3dt0001gs3)[0])
# s3disterrdatarspgs3dt00001raw = np.zeros(np.shape(datars3dt00001gs3)[0])

# s3distgs3dt1     = np.zeros(np.shape(datars3dt1gs3)[0])
# s3distgs3dt01    = np.zeros(np.shape(datars3dt01gs3)[0])
# s3distgs3dt001   = np.zeros(np.shape(datars3dt001gs3)[0])
# s3distgs3dt0001  = np.zeros(np.shape(datars3dt0001gs3)[0])
# s3distgs3dt00001 = np.zeros(np.shape(datars3dt00001gs3)[0])


# counter = 0
# print("Collating dt = .1 Data")
# for a in range(np.shape(datarh3dt1gs1)[0]):
#     h3disterrdatarspgs1dt1raw[counter] = h3dist(rot2hyp(datarh3dt1gs1[a][0:3]),rot2hyp(datarh3dt1gs1[a][3:6]))
#     h3disterrdatarspgs2dt1raw[counter] = h3dist(rot2hyp(datarh3dt1gs2[a][0:3]),rot2hyp(datarh3dt1gs2[a][3:6]))
#     h3disterrdatarspgs3dt1raw[counter] = h3dist(rot2hyp(datarh3dt1gs3[a][0:3]),rot2hyp(datarh3dt1gs3[a][3:6]))

#     h3distgs1dt1[counter] = arccosh(cosh(datarh3dt1gs1[a][0])/cosh(.5*h3disterrdatarspgs1dt1raw[counter]))
#     h3distgs2dt1[counter] = arccosh(cosh(datarh3dt1gs2[a][0])/cosh(.5*h3disterrdatarspgs2dt1raw[counter]))
#     h3distgs3dt1[counter] = arccosh(cosh(datarh3dt1gs3[a][0])/cosh(.5*h3disterrdatarspgs3dt1raw[counter]))

#     s3disterrdatarspgs1dt1raw[counter] = r4dist(rot2r4(datars3dt1gs1[a][0:3]),rot2r4(datars3dt1gs1[a][3:6]))
#     s3disterrdatarspgs2dt1raw[counter] = r4dist(rot2r4(datars3dt1gs2[a][0:3]),rot2r4(datars3dt1gs2[a][3:6]))
#     s3disterrdatarspgs3dt1raw[counter] = r4dist(rot2r4(datars3dt1gs3[a][0:3]),rot2r4(datars3dt1gs3[a][3:6]))

#     s3distgs1dt1[counter] = arccos(cos(datars3dt1gs1[a][0])/cos(.5*s3disterrdatarspgs1dt1raw[counter]))
#     s3distgs2dt1[counter] = arccos(cos(datars3dt1gs2[a][0])/cos(.5*s3disterrdatarspgs2dt1raw[counter]))
#     s3distgs3dt1[counter] = arccos(cos(datars3dt1gs3[a][0])/cos(.5*s3disterrdatarspgs3dt1raw[counter]))

#     counter += 1

# h3disterrdatarspgs1dt1 = (2.*datarh3gtdt00001gs3[::10000,0] - h3disterrdatarspgs1dt1raw)/(2.*datarh3gtdt00001gs3[::10000,0])
# h3disterrdatarspgs2dt1 = (2.*datarh3gtdt00001gs3[::10000,0] - h3disterrdatarspgs2dt1raw)/(2.*datarh3gtdt00001gs3[::10000,0])
# h3disterrdatarspgs3dt1 = (2.*datarh3gtdt00001gs3[::10000,0] - h3disterrdatarspgs3dt1raw)/(2.*datarh3gtdt00001gs3[::10000,0])

# s3disterrdatarspgs1dt1 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]) - s3disterrdatarspgs1dt1raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]))
# s3disterrdatarspgs2dt1 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]) - s3disterrdatarspgs2dt1raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]))
# s3disterrdatarspgs3dt1 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]) - s3disterrdatarspgs3dt1raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]))

# counter = 0
# print("Collating dt = .01 Data")
# for a in range(np.shape(datarh3dt01gs1)[0]):
#     h3disterrdatarspgs1dt01raw[counter] = h3dist(rot2hyp(datarh3dt01gs1[a][0:3]),rot2hyp(datarh3dt01gs1[a][3:6]))
#     h3disterrdatarspgs2dt01raw[counter] = h3dist(rot2hyp(datarh3dt01gs2[a][0:3]),rot2hyp(datarh3dt01gs2[a][3:6]))
#     h3disterrdatarspgs3dt01raw[counter] = h3dist(rot2hyp(datarh3dt01gs3[a][0:3]),rot2hyp(datarh3dt01gs3[a][3:6]))

#     h3distgs1dt01[counter] = arccosh(cosh(datarh3dt01gs1[a][0])/cosh(.5*h3disterrdatarspgs1dt01raw[counter]))
#     h3distgs2dt01[counter] = arccosh(cosh(datarh3dt01gs2[a][0])/cosh(.5*h3disterrdatarspgs2dt01raw[counter]))
#     h3distgs3dt01[counter] = arccosh(cosh(datarh3dt01gs3[a][0])/cosh(.5*h3disterrdatarspgs3dt01raw[counter]))

#     s3disterrdatarspgs1dt01raw[counter] = r4dist(rot2r4(datars3dt01gs1[a][0:3]),rot2r4(datars3dt01gs1[a][3:6]))
#     s3disterrdatarspgs2dt01raw[counter] = r4dist(rot2r4(datars3dt01gs2[a][0:3]),rot2r4(datars3dt01gs2[a][3:6]))
#     s3disterrdatarspgs3dt01raw[counter] = r4dist(rot2r4(datars3dt01gs3[a][0:3]),rot2r4(datars3dt01gs3[a][3:6]))

#     s3distgs1dt01[counter] = arccos(cos(datars3dt01gs1[a][0])/cos(.5*s3disterrdatarspgs1dt01raw[counter]))
#     s3distgs2dt01[counter] = arccos(cos(datars3dt01gs2[a][0])/cos(.5*s3disterrdatarspgs2dt01raw[counter]))
#     s3distgs3dt01[counter] = arccos(cos(datars3dt01gs3[a][0])/cos(.5*s3disterrdatarspgs3dt01raw[counter]))

#     counter += 1

# h3disterrdatarspgs1dt01 = (2.*datarh3gtdt00001gs3[::1000,0] - h3disterrdatarspgs1dt01raw)/(2.*datarh3gtdt00001gs3[::1000,0])
# h3disterrdatarspgs2dt01 = (2.*datarh3gtdt00001gs3[::1000,0] - h3disterrdatarspgs2dt01raw)/(2.*datarh3gtdt00001gs3[::1000,0])
# h3disterrdatarspgs3dt01 = (2.*datarh3gtdt00001gs3[::1000,0] - h3disterrdatarspgs3dt01raw)/(2.*datarh3gtdt00001gs3[::1000,0])

# s3disterrdatarspgs1dt01 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]) - s3disterrdatarspgs1dt01raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]))
# s3disterrdatarspgs2dt01 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]) - s3disterrdatarspgs2dt01raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]))
# s3disterrdatarspgs3dt01 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]) - s3disterrdatarspgs3dt01raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]))

# counter = 0
# print("Collating dt = .001 Data")
# for a in range(np.shape(datarh3dt001gs1)[0]):
#     h3disterrdatarspgs1dt001raw[counter] = h3dist(rot2hyp(datarh3dt001gs1[a][0:3]),rot2hyp(datarh3dt001gs1[a][3:6]))
#     h3disterrdatarspgs2dt001raw[counter] = h3dist(rot2hyp(datarh3dt001gs2[a][0:3]),rot2hyp(datarh3dt001gs2[a][3:6]))
#     h3disterrdatarspgs3dt001raw[counter] = h3dist(rot2hyp(datarh3dt001gs3[a][0:3]),rot2hyp(datarh3dt001gs3[a][3:6]))

#     h3distgs1dt001[counter] = arccosh(cosh(datarh3dt001gs1[a][0])/cosh(.5*h3disterrdatarspgs1dt001raw[counter]))
#     h3distgs2dt001[counter] = arccosh(cosh(datarh3dt001gs2[a][0])/cosh(.5*h3disterrdatarspgs2dt001raw[counter]))
#     h3distgs3dt001[counter] = arccosh(cosh(datarh3dt001gs3[a][0])/cosh(.5*h3disterrdatarspgs3dt001raw[counter]))

#     s3disterrdatarspgs1dt001raw[counter] = r4dist(rot2r4(datars3dt001gs1[a][0:3]),rot2r4(datars3dt001gs1[a][3:6]))
#     s3disterrdatarspgs2dt001raw[counter] = r4dist(rot2r4(datars3dt001gs2[a][0:3]),rot2r4(datars3dt001gs2[a][3:6]))
#     s3disterrdatarspgs3dt001raw[counter] = r4dist(rot2r4(datars3dt001gs3[a][0:3]),rot2r4(datars3dt001gs3[a][3:6]))

#     s3distgs1dt001[counter] = arccos(cos(datars3dt001gs1[a][0])/cos(.5*s3disterrdatarspgs1dt001raw[counter]))
#     s3distgs2dt001[counter] = arccos(cos(datars3dt001gs2[a][0])/cos(.5*s3disterrdatarspgs2dt001raw[counter]))
#     s3distgs3dt001[counter] = arccos(cos(datars3dt001gs3[a][0])/cos(.5*s3disterrdatarspgs3dt001raw[counter]))

#     counter += 1

# h3disterrdatarspgs1dt001 = (2.*datarh3gtdt00001gs3[::100,0] - h3disterrdatarspgs1dt001raw)/(2.*datarh3gtdt00001gs3[::100,0])
# h3disterrdatarspgs2dt001 = (2.*datarh3gtdt00001gs3[::100,0] - h3disterrdatarspgs2dt001raw)/(2.*datarh3gtdt00001gs3[::100,0])
# h3disterrdatarspgs3dt001 = (2.*datarh3gtdt00001gs3[::100,0] - h3disterrdatarspgs3dt001raw)/(2.*datarh3gtdt00001gs3[::100,0])

# s3disterrdatarspgs1dt001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]) - s3disterrdatarspgs1dt001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]))
# s3disterrdatarspgs2dt001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]) - s3disterrdatarspgs2dt001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]))
# s3disterrdatarspgs3dt001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]) - s3disterrdatarspgs3dt001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]))

# counter = 0
# print("Collating dt = .0001 Data")
# for a in range(np.shape(datarh3dt0001gs1)[0]):
#     h3disterrdatarspgs1dt0001raw[counter] = h3dist(rot2hyp(datarh3dt0001gs1[a][0:3]),rot2hyp(datarh3dt0001gs1[a][3:6]))
#     h3disterrdatarspgs2dt0001raw[counter] = h3dist(rot2hyp(datarh3dt0001gs2[a][0:3]),rot2hyp(datarh3dt0001gs2[a][3:6]))
#     h3disterrdatarspgs3dt0001raw[counter] = h3dist(rot2hyp(datarh3dt0001gs3[a][0:3]),rot2hyp(datarh3dt0001gs3[a][3:6]))

#     h3distgs1dt0001[counter] = arccosh(cosh(datarh3dt0001gs1[a][0])/cosh(.5*h3disterrdatarspgs1dt0001raw[counter]))
#     h3distgs2dt0001[counter] = arccosh(cosh(datarh3dt0001gs2[a][0])/cosh(.5*h3disterrdatarspgs2dt0001raw[counter]))
#     h3distgs3dt0001[counter] = arccosh(cosh(datarh3dt0001gs3[a][0])/cosh(.5*h3disterrdatarspgs3dt0001raw[counter]))

#     s3disterrdatarspgs1dt0001raw[counter] = r4dist(rot2r4(datars3dt0001gs1[a][0:3]),rot2r4(datars3dt0001gs1[a][3:6]))
#     s3disterrdatarspgs2dt0001raw[counter] = r4dist(rot2r4(datars3dt0001gs2[a][0:3]),rot2r4(datars3dt0001gs2[a][3:6]))
#     s3disterrdatarspgs3dt0001raw[counter] = r4dist(rot2r4(datars3dt0001gs3[a][0:3]),rot2r4(datars3dt0001gs3[a][3:6]))

#     s3distgs1dt0001[counter] = arccos(cos(datars3dt0001gs1[a][0])/cos(.5*s3disterrdatarspgs1dt0001raw[counter]))
#     s3distgs2dt0001[counter] = arccos(cos(datars3dt0001gs2[a][0])/cos(.5*s3disterrdatarspgs2dt0001raw[counter]))
#     s3distgs3dt0001[counter] = arccos(cos(datars3dt0001gs3[a][0])/cos(.5*s3disterrdatarspgs3dt0001raw[counter]))

#     counter += 1

# h3disterrdatarspgs1dt0001 = (2.*datarh3gtdt00001gs3[::10,0] - h3disterrdatarspgs1dt0001raw)/(2.*datarh3gtdt00001gs3[::10,0])
# h3disterrdatarspgs2dt0001 = (2.*datarh3gtdt00001gs3[::10,0] - h3disterrdatarspgs2dt0001raw)/(2.*datarh3gtdt00001gs3[::10,0])
# h3disterrdatarspgs3dt0001 = (2.*datarh3gtdt00001gs3[::10,0] - h3disterrdatarspgs3dt0001raw)/(2.*datarh3gtdt00001gs3[::10,0])

# s3disterrdatarspgs1dt0001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]) - s3disterrdatarspgs1dt0001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]))
# s3disterrdatarspgs2dt0001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]) - s3disterrdatarspgs2dt0001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]))
# s3disterrdatarspgs3dt0001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]) - s3disterrdatarspgs3dt0001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]))

# counter = 0
# print("Collating dt = .00001 Data")
# for a in range(np.shape(datarh3dt00001gs1)[0]):
#     h3disterrdatarspgs1dt00001raw[counter] = h3dist(rot2hyp(datarh3dt00001gs1[a][0:3]),rot2hyp(datarh3dt00001gs1[a][3:6]))
#     h3disterrdatarspgs2dt00001raw[counter] = h3dist(rot2hyp(datarh3dt00001gs2[a][0:3]),rot2hyp(datarh3dt00001gs2[a][3:6]))
#     h3disterrdatarspgs3dt00001raw[counter] = h3dist(rot2hyp(datarh3dt00001gs3[a][0:3]),rot2hyp(datarh3dt00001gs3[a][3:6]))

#     h3distgs1dt00001[counter] = arccosh(cosh(datarh3dt00001gs1[a][0])/cosh(.5*h3disterrdatarspgs1dt00001raw[counter]))
#     h3distgs2dt00001[counter] = arccosh(cosh(datarh3dt00001gs2[a][0])/cosh(.5*h3disterrdatarspgs2dt00001raw[counter]))
#     h3distgs3dt00001[counter] = arccosh(cosh(datarh3dt00001gs3[a][0])/cosh(.5*h3disterrdatarspgs3dt00001raw[counter]))

#     s3disterrdatarspgs1dt00001raw[counter] = r4dist(rot2r4(datars3dt00001gs1[a][0:3]),rot2r4(datars3dt00001gs1[a][3:6]))
#     s3disterrdatarspgs2dt00001raw[counter] = r4dist(rot2r4(datars3dt00001gs2[a][0:3]),rot2r4(datars3dt00001gs2[a][3:6]))
#     s3disterrdatarspgs3dt00001raw[counter] = r4dist(rot2r4(datars3dt00001gs3[a][0:3]),rot2r4(datars3dt00001gs3[a][3:6]))

#     s3distgs1dt00001[counter] = arccos(cos(datars3dt00001gs1[a][0])/cos(.5*s3disterrdatarspgs1dt00001raw[counter]))
#     s3distgs2dt00001[counter] = arccos(cos(datars3dt00001gs2[a][0])/cos(.5*s3disterrdatarspgs2dt00001raw[counter]))
#     s3distgs3dt00001[counter] = arccos(cos(datars3dt00001gs3[a][0])/cos(.5*s3disterrdatarspgs3dt00001raw[counter]))

#     counter += 1

# h3disterrdatarspgs1dt00001 = (2.*datarh3gtdt00001gs3[::1,0] - h3disterrdatarspgs1dt00001raw)/(2.*datarh3gtdt00001gs3[::1,0])
# h3disterrdatarspgs2dt00001 = (2.*datarh3gtdt00001gs3[::1,0] - h3disterrdatarspgs2dt00001raw)/(2.*datarh3gtdt00001gs3[::1,0])
# h3disterrdatarspgs3dt00001 = (2.*datarh3gtdt00001gs3[::1,0] - h3disterrdatarspgs3dt00001raw)/(2.*datarh3gtdt00001gs3[::1,0])

# s3disterrdatarspgs1dt00001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]) - s3disterrdatarspgs1dt00001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]))
# s3disterrdatarspgs2dt00001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]) - s3disterrdatarspgs2dt00001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]))
# s3disterrdatarspgs3dt00001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]) - s3disterrdatarspgs3dt00001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]))


# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# ax.plot(t_arr00001,h3disterrdatarspgs1dt00001, color='r',label = r'$\mathbf{H}^3$ gs1')
# ax.plot(t_arr00001,h3disterrdatarspgs2dt00001, color='b',label = r'$\mathbf{H}^3$ gs2')
# ax.plot(t_arr00001,h3disterrdatarspgs3dt00001, color='k',label = r'$\mathbf{H}^3$ gs3')

# # ax.plot(t_arr00001,s3disterrdatarspgs1dt00001, color='r',label = r'$\mathbf{S}^3$ gs1')
# # ax.plot(t_arr00001,s3disterrdatarspgs2dt00001, color='b',label = r'$\mathbf{S}^3$ gs2')
# # ax.plot(t_arr00001,s3disterrdatarspgs3dt00001, color='k',label = r'$\mathbf{S}^3$ gs3')

# # h3 
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-2}'
# # ylabels[1] = r'\textbf{-2.0}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{2.0}'
# # ylabels[4] = r'\textbf{4.0}'
# # ylabels[5] = r'\textbf{6}'
# # ylabels[6] = r'\textbf{8}'
# # ylabels[7] = r'\textbf{10}'

# # h3 zoom
# ylabels = [item.get_text() for item in ax.get_yticklabels()]
# ylabels[0] = r'\textbf{-1.0}'
# ylabels[1] = r'\textbf{-0.5}'
# ylabels[2] = r'\textbf{0.0}'
# ylabels[3] = r'\textbf{0.5}'
# ylabels[4] = r'\textbf{1.0}'
# # ylabels[5] = r'\textbf{6}'
# # ylabels[6] = r'\textbf{8}'
# # ylabels[7] = r'\textbf{10}'
# xlabels = [item.get_text() for item in ax.get_xticklabels()]
# xlabels[0] = r'\textbf{0.0}'
# xlabels[1] = r'\textbf{2.0}'
# xlabels[2] = r'\textbf{4.0}'
# xlabels[3] = r'\textbf{6.0}'
# xlabels[4] = r'\textbf{8.0}'
# xlabels[5] = r'\textbf{5.0}'
# xlabels[6] = r'\textbf{6.0}'
# xlabels[7] = r'\textbf{7.0}'

# # s3 
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-1.0}'
# # ylabels[1] = r'\textbf{-1.0}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{1.0}'
# # ylabels[4] = r'\textbf{1.0}'
# # ylabels[5] = r'\textbf{6}'
# # ylabels[6] = r'\textbf{8}'
# # ylabels[7] = r'\textbf{10}'

# # s3 
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-1.0}'
# # ylabels[1] = r'\textbf{-0.5}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{0.5}'
# # ylabels[4] = r'\textbf{1.0}'
# # ylabels[5] = r'\textbf{6}'
# # ylabels[6] = r'\textbf{8}'
# # ylabels[7] = r'\textbf{10}'

# ax.set_xticklabels(xlabels)
# ax.set_yticklabels(ylabels)
# ax.legend(fontsize="20",loc="upper left")
# plt.title(r'\textbf{Elastic Rod Body in $\mathbf{H}^3$}', fontsize=20)
# plt.ylabel(r'\textbf{Relative Error of \\ Vertex Separation ($\times 10^{-11}$)}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)

# # For zoom plots
# ax.set_xlim(0.0,7.0)
# ax.set_ylim(0-1e-11,0+1e-11) # h3
# # ax.set_ylim(0-1e-12,0+1e-12) # s3

# plt.tight_layout()

# plt.show()


##### old stuff

# ax.plot(h3distgs1dt00001,h3disterrdatarspgs1dt00001, color='r',label = "H3 Gauss s1 dt00001")
# ax.plot(h3distgs2dt00001,h3disterrdatarspgs2dt00001, color='b',label = "H3 Gauss s2 dt00001")
# ax.plot(h3distgs3dt00001,h3disterrdatarspgs3dt00001, color='k',label = "H3 Gauss s3 dt00001")

# ax.plot(s3distgs1dt1,s3disterrdatarspgs1dt1, color='r',label = "S3 Gauss s1 dt1")
# ax.plot(s3distgs2dt1,s3disterrdatarspgs2dt1, color='b',label = "S3 Gauss s2 dt1")
# ax.plot(s3distgs3dt1,s3disterrdatarspgs3dt1, color='k',label = "S3 Gauss s3 dt1")

# ax.plot(t_arr0001,h3disterrdatarspgs1dt0001, color='r',label = "H3 Gauss s1 dt0001")
# ax.plot(t_arr0001,h3disterrdatarspgs2dt0001, color='b',label = "H3 Gauss s2 dt0001")
# ax.plot(t_arr0001,h3disterrdatarspgs3dt0001, color='k',label = "H3 Gauss s3 dt0001")

# ax.plot(t_arr00001,s3disterrdatarspgs1dt00001, color='r',label = "S3 Gauss s1 dt00001")
# ax.plot(t_arr00001,s3disterrdatarspgs2dt00001, color='b',label = "S3 Gauss s2 dt00001")
# ax.plot(t_arr00001,s3disterrdatarspgs3dt00001, color='k',label = "S3 Gauss s3 dt00001")

# ax.scatter(timevals,s3errvalss1,marker = 'o',color='r',label = "S3 Gauss s1")
# ax.scatter(timevals,s3errvalss2,marker = 'o',color='b',label = "S3 Gauss s2")
# ax.scatter(timevals,s3errvalss3,marker = 'o',color='k',label = "S3 Gauss s3")

# ax.legend()
# # ax.set_title('Vertex Separation vs. Time')
# ax.set_title('Vertex Separation vs. Midpoint Distance')
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# ax.set_xlim(0,5)
# ax.set_ylim(0-1e-11,0+1e-11)
# # ax.set_xlabel('Time')
# ax.set_xlabel('Midpoint Distance')
# ax.set_ylabel('Vertex Separation')
# plt.show()





######
# Spring Potential Rod Separation Error Scaling plot for H3 and S3
######

# ### Gauss s1

# h3disterrdatarspgs1dt1raw     = np.zeros(np.shape(datarh3dt1gs1)[0])
# h3disterrdatarspgs1dt01raw    = np.zeros(np.shape(datarh3dt01gs1)[0])
# h3disterrdatarspgs1dt001raw   = np.zeros(np.shape(datarh3dt001gs1)[0])
# h3disterrdatarspgs1dt0001raw  = np.zeros(np.shape(datarh3dt0001gs1)[0])
# h3disterrdatarspgs1dt00001raw = np.zeros(np.shape(datarh3dt00001gs1)[0])

# s3disterrdatarspgs1dt1raw     = np.zeros(np.shape(datars3dt1gs1)[0])
# s3disterrdatarspgs1dt01raw    = np.zeros(np.shape(datars3dt01gs1)[0])
# s3disterrdatarspgs1dt001raw   = np.zeros(np.shape(datars3dt001gs1)[0])
# s3disterrdatarspgs1dt0001raw  = np.zeros(np.shape(datars3dt0001gs1)[0])
# s3disterrdatarspgs1dt00001raw = np.zeros(np.shape(datars3dt00001gs1)[0])

# ### Gauss s2

# h3disterrdatarspgs2dt1raw     = np.zeros(np.shape(datarh3dt1gs2)[0])
# h3disterrdatarspgs2dt01raw    = np.zeros(np.shape(datarh3dt01gs2)[0])
# h3disterrdatarspgs2dt001raw   = np.zeros(np.shape(datarh3dt001gs2)[0])
# h3disterrdatarspgs2dt0001raw  = np.zeros(np.shape(datarh3dt0001gs2)[0])
# h3disterrdatarspgs2dt00001raw = np.zeros(np.shape(datarh3dt00001gs2)[0])

# s3disterrdatarspgs2dt1raw     = np.zeros(np.shape(datars3dt1gs2)[0])
# s3disterrdatarspgs2dt01raw    = np.zeros(np.shape(datars3dt01gs2)[0])
# s3disterrdatarspgs2dt001raw   = np.zeros(np.shape(datars3dt001gs2)[0])
# s3disterrdatarspgs2dt0001raw  = np.zeros(np.shape(datars3dt0001gs2)[0])
# s3disterrdatarspgs2dt00001raw = np.zeros(np.shape(datars3dt00001gs2)[0])

# ### Gauss s3

# h3disterrdatarspgs3dt1raw     = np.zeros(np.shape(datarh3dt1gs3)[0])
# h3disterrdatarspgs3dt01raw    = np.zeros(np.shape(datarh3dt01gs3)[0])
# h3disterrdatarspgs3dt001raw   = np.zeros(np.shape(datarh3dt001gs3)[0])
# h3disterrdatarspgs3dt0001raw  = np.zeros(np.shape(datarh3dt0001gs3)[0])
# h3disterrdatarspgs3dt00001raw = np.zeros(np.shape(datarh3dt00001gs3)[0])

# s3disterrdatarspgs3dt1raw     = np.zeros(np.shape(datars3dt1gs3)[0])
# s3disterrdatarspgs3dt01raw    = np.zeros(np.shape(datars3dt01gs3)[0])
# s3disterrdatarspgs3dt001raw   = np.zeros(np.shape(datars3dt001gs3)[0])
# s3disterrdatarspgs3dt0001raw  = np.zeros(np.shape(datars3dt0001gs3)[0])
# s3disterrdatarspgs3dt00001raw = np.zeros(np.shape(datars3dt00001gs3)[0])


# counter = 0
# print("Collating dt = .1 Data")
# for a in range(np.shape(datarh3dt1gs1)[0]):
#     h3disterrdatarspgs1dt1raw[counter] = h3dist(rot2hyp(datarh3dt1gs1[a][0:3]),rot2hyp(datarh3dt1gs1[a][3:6]))
#     h3disterrdatarspgs2dt1raw[counter] = h3dist(rot2hyp(datarh3dt1gs2[a][0:3]),rot2hyp(datarh3dt1gs2[a][3:6]))
#     h3disterrdatarspgs3dt1raw[counter] = h3dist(rot2hyp(datarh3dt1gs3[a][0:3]),rot2hyp(datarh3dt1gs3[a][3:6]))

#     s3disterrdatarspgs1dt1raw[counter] = r4dist(rot2r4(datars3dt1gs1[a][0:3]),rot2r4(datars3dt1gs1[a][3:6]))
#     s3disterrdatarspgs2dt1raw[counter] = r4dist(rot2r4(datars3dt1gs2[a][0:3]),rot2r4(datars3dt1gs2[a][3:6]))
#     s3disterrdatarspgs3dt1raw[counter] = r4dist(rot2r4(datars3dt1gs3[a][0:3]),rot2r4(datars3dt1gs3[a][3:6]))

#     counter += 1

# h3disterrdatarspgs1dt1 = (2.*datarh3gtdt00001gs3[::10000,0] - h3disterrdatarspgs1dt1raw)/(2.*datarh3gtdt00001gs3[::10000,0])
# h3disterrdatarspgs2dt1 = (2.*datarh3gtdt00001gs3[::10000,0] - h3disterrdatarspgs2dt1raw)/(2.*datarh3gtdt00001gs3[::10000,0])
# h3disterrdatarspgs3dt1 = (2.*datarh3gtdt00001gs3[::10000,0] - h3disterrdatarspgs3dt1raw)/(2.*datarh3gtdt00001gs3[::10000,0])

# s3disterrdatarspgs1dt1 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]) - s3disterrdatarspgs1dt1raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]))
# s3disterrdatarspgs2dt1 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]) - s3disterrdatarspgs2dt1raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]))
# s3disterrdatarspgs3dt1 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]) - s3disterrdatarspgs3dt1raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10000,0]))

# counter = 0
# print("Collating dt = .01 Data")
# for a in range(np.shape(datarh3dt01gs1)[0]):
#     h3disterrdatarspgs1dt01raw[counter] = h3dist(rot2hyp(datarh3dt01gs1[a][0:3]),rot2hyp(datarh3dt01gs1[a][3:6]))
#     h3disterrdatarspgs2dt01raw[counter] = h3dist(rot2hyp(datarh3dt01gs2[a][0:3]),rot2hyp(datarh3dt01gs2[a][3:6]))
#     h3disterrdatarspgs3dt01raw[counter] = h3dist(rot2hyp(datarh3dt01gs3[a][0:3]),rot2hyp(datarh3dt01gs3[a][3:6]))

#     s3disterrdatarspgs1dt01raw[counter] = r4dist(rot2r4(datars3dt01gs1[a][0:3]),rot2r4(datars3dt01gs1[a][3:6]))
#     s3disterrdatarspgs2dt01raw[counter] = r4dist(rot2r4(datars3dt01gs2[a][0:3]),rot2r4(datars3dt01gs2[a][3:6]))
#     s3disterrdatarspgs3dt01raw[counter] = r4dist(rot2r4(datars3dt01gs3[a][0:3]),rot2r4(datars3dt01gs3[a][3:6]))

#     counter += 1

# h3disterrdatarspgs1dt01 = (2.*datarh3gtdt00001gs3[::1000,0] - h3disterrdatarspgs1dt01raw)/(2.*datarh3gtdt00001gs3[::1000,0])
# h3disterrdatarspgs2dt01 = (2.*datarh3gtdt00001gs3[::1000,0] - h3disterrdatarspgs2dt01raw)/(2.*datarh3gtdt00001gs3[::1000,0])
# h3disterrdatarspgs3dt01 = (2.*datarh3gtdt00001gs3[::1000,0] - h3disterrdatarspgs3dt01raw)/(2.*datarh3gtdt00001gs3[::1000,0])

# s3disterrdatarspgs1dt01 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]) - s3disterrdatarspgs1dt01raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]))
# s3disterrdatarspgs2dt01 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]) - s3disterrdatarspgs2dt01raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]))
# s3disterrdatarspgs3dt01 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]) - s3disterrdatarspgs3dt01raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1000,0]))

# counter = 0
# print("Collating dt = .001 Data")
# for a in range(np.shape(datarh3dt001gs1)[0]):
#     h3disterrdatarspgs1dt001raw[counter] = h3dist(rot2hyp(datarh3dt001gs1[a][0:3]),rot2hyp(datarh3dt001gs1[a][3:6]))
#     h3disterrdatarspgs2dt001raw[counter] = h3dist(rot2hyp(datarh3dt001gs2[a][0:3]),rot2hyp(datarh3dt001gs2[a][3:6]))
#     h3disterrdatarspgs3dt001raw[counter] = h3dist(rot2hyp(datarh3dt001gs3[a][0:3]),rot2hyp(datarh3dt001gs3[a][3:6]))

#     s3disterrdatarspgs1dt001raw[counter] = r4dist(rot2r4(datars3dt001gs1[a][0:3]),rot2r4(datars3dt001gs1[a][3:6]))
#     s3disterrdatarspgs2dt001raw[counter] = r4dist(rot2r4(datars3dt001gs2[a][0:3]),rot2r4(datars3dt001gs2[a][3:6]))
#     s3disterrdatarspgs3dt001raw[counter] = r4dist(rot2r4(datars3dt001gs3[a][0:3]),rot2r4(datars3dt001gs3[a][3:6]))

#     counter += 1

# h3disterrdatarspgs1dt001 = (2.*datarh3gtdt00001gs3[::100,0] - h3disterrdatarspgs1dt001raw)/(2.*datarh3gtdt00001gs3[::100,0])
# h3disterrdatarspgs2dt001 = (2.*datarh3gtdt00001gs3[::100,0] - h3disterrdatarspgs2dt001raw)/(2.*datarh3gtdt00001gs3[::100,0])
# h3disterrdatarspgs3dt001 = (2.*datarh3gtdt00001gs3[::100,0] - h3disterrdatarspgs3dt001raw)/(2.*datarh3gtdt00001gs3[::100,0])

# s3disterrdatarspgs1dt001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]) - s3disterrdatarspgs1dt001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]))
# s3disterrdatarspgs2dt001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]) - s3disterrdatarspgs2dt001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]))
# s3disterrdatarspgs3dt001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]) - s3disterrdatarspgs3dt001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::100,0]))

# counter = 0
# print("Collating dt = .0001 Data")
# for a in range(np.shape(datarh3dt0001gs1)[0]):
#     h3disterrdatarspgs1dt0001raw[counter] = h3dist(rot2hyp(datarh3dt0001gs1[a][0:3]),rot2hyp(datarh3dt0001gs1[a][3:6]))
#     h3disterrdatarspgs2dt0001raw[counter] = h3dist(rot2hyp(datarh3dt0001gs2[a][0:3]),rot2hyp(datarh3dt0001gs2[a][3:6]))
#     h3disterrdatarspgs3dt0001raw[counter] = h3dist(rot2hyp(datarh3dt0001gs3[a][0:3]),rot2hyp(datarh3dt0001gs3[a][3:6]))

#     s3disterrdatarspgs1dt0001raw[counter] = r4dist(rot2r4(datars3dt0001gs1[a][0:3]),rot2r4(datars3dt0001gs1[a][3:6]))
#     s3disterrdatarspgs2dt0001raw[counter] = r4dist(rot2r4(datars3dt0001gs2[a][0:3]),rot2r4(datars3dt0001gs2[a][3:6]))
#     s3disterrdatarspgs3dt0001raw[counter] = r4dist(rot2r4(datars3dt0001gs3[a][0:3]),rot2r4(datars3dt0001gs3[a][3:6]))

#     counter += 1

# h3disterrdatarspgs1dt0001 = (2.*datarh3gtdt00001gs3[::10,0] - h3disterrdatarspgs1dt0001raw)/(2.*datarh3gtdt00001gs3[::10,0])
# h3disterrdatarspgs2dt0001 = (2.*datarh3gtdt00001gs3[::10,0] - h3disterrdatarspgs2dt0001raw)/(2.*datarh3gtdt00001gs3[::10,0])
# h3disterrdatarspgs3dt0001 = (2.*datarh3gtdt00001gs3[::10,0] - h3disterrdatarspgs3dt0001raw)/(2.*datarh3gtdt00001gs3[::10,0])

# s3disterrdatarspgs1dt0001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]) - s3disterrdatarspgs1dt0001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]))
# s3disterrdatarspgs2dt0001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]) - s3disterrdatarspgs2dt0001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]))
# s3disterrdatarspgs3dt0001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]) - s3disterrdatarspgs3dt0001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::10,0]))

# counter = 0
# print("Collating dt = .00001 Data")
# for a in range(np.shape(datarh3dt00001gs1)[0]):
#     h3disterrdatarspgs1dt00001raw[counter] = h3dist(rot2hyp(datarh3dt00001gs1[a][0:3]),rot2hyp(datarh3dt00001gs1[a][3:6]))
#     h3disterrdatarspgs2dt00001raw[counter] = h3dist(rot2hyp(datarh3dt00001gs2[a][0:3]),rot2hyp(datarh3dt00001gs2[a][3:6]))
#     h3disterrdatarspgs3dt00001raw[counter] = h3dist(rot2hyp(datarh3dt00001gs3[a][0:3]),rot2hyp(datarh3dt00001gs3[a][3:6]))

#     s3disterrdatarspgs1dt00001raw[counter] = r4dist(rot2r4(datars3dt00001gs1[a][0:3]),rot2r4(datars3dt00001gs1[a][3:6]))
#     s3disterrdatarspgs2dt00001raw[counter] = r4dist(rot2r4(datars3dt00001gs2[a][0:3]),rot2r4(datars3dt00001gs2[a][3:6]))
#     s3disterrdatarspgs3dt00001raw[counter] = r4dist(rot2r4(datars3dt00001gs3[a][0:3]),rot2r4(datars3dt00001gs3[a][3:6]))

#     counter += 1

# h3disterrdatarspgs1dt00001 = (2.*datarh3gtdt00001gs3[::1,0] - h3disterrdatarspgs1dt00001raw)/(2.*datarh3gtdt00001gs3[::1,0])
# h3disterrdatarspgs2dt00001 = (2.*datarh3gtdt00001gs3[::1,0] - h3disterrdatarspgs2dt00001raw)/(2.*datarh3gtdt00001gs3[::1,0])
# h3disterrdatarspgs3dt00001 = (2.*datarh3gtdt00001gs3[::1,0] - h3disterrdatarspgs3dt00001raw)/(2.*datarh3gtdt00001gs3[::1,0])

# s3disterrdatarspgs1dt00001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]) - s3disterrdatarspgs1dt00001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]))
# s3disterrdatarspgs2dt00001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]) - s3disterrdatarspgs2dt00001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]))
# s3disterrdatarspgs3dt00001 = (2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]) - s3disterrdatarspgs3dt00001raw)/(2.*(np.pi/2. - datars3gtdt00001gs3[::1,0]))

# timevals = [.1,.01,.001,.0001,.00001]

# timevals2 = [.1**2,.01**2,.001**2,.0001**2,.00001**2]
# timevals4 = [.1**4,.01**4,.001**4,.0001**4,.00001**4]
# timevals6 = [.1**6,.01**6,.001**6,.0001**6,.00001**6]

# h3errvalss1 = [max(abs(h3disterrdatarspgs1dt1)),max(abs(h3disterrdatarspgs1dt01)),max(abs(h3disterrdatarspgs1dt001)),max(abs(h3disterrdatarspgs1dt0001)),max(abs(h3disterrdatarspgs1dt00001))]
# h3errvalss2 = [max(abs(h3disterrdatarspgs2dt1)),max(abs(h3disterrdatarspgs2dt01)),max(abs(h3disterrdatarspgs2dt001)),max(abs(h3disterrdatarspgs2dt0001)),max(abs(h3disterrdatarspgs2dt00001))]
# h3errvalss3 = [max(abs(h3disterrdatarspgs3dt1)),max(abs(h3disterrdatarspgs3dt01)),max(abs(h3disterrdatarspgs3dt001)),max(abs(h3disterrdatarspgs3dt0001)),max(abs(h3disterrdatarspgs3dt00001))]

# s3errvalss1 = [max(abs(s3disterrdatarspgs1dt1)),max(abs(s3disterrdatarspgs1dt01)),max(abs(s3disterrdatarspgs1dt001)),max(abs(s3disterrdatarspgs1dt0001)),max(abs(s3disterrdatarspgs1dt00001))]
# s3errvalss2 = [max(abs(s3disterrdatarspgs2dt1)),max(abs(s3disterrdatarspgs2dt01)),max(abs(s3disterrdatarspgs2dt001)),max(abs(s3disterrdatarspgs2dt0001)),max(abs(s3disterrdatarspgs2dt00001))]
# s3errvalss3 = [max(abs(s3disterrdatarspgs3dt1)),max(abs(s3disterrdatarspgs3dt01)),max(abs(s3disterrdatarspgs3dt001)),max(abs(s3disterrdatarspgs3dt0001)),max(abs(s3disterrdatarspgs3dt00001))]


# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# # ax.scatter(timevals,h3errvalss1,marker = 's', color='r',label = r"$\mathbf{H}^3$ gs1")
# # ax.scatter(timevals,h3errvalss2,marker = 's', color='b',label = r"$\mathbf{H}^3$ gs2")
# # ax.scatter(timevals,h3errvalss3,marker = 's', color='k',label = r"$\mathbf{H}^3$ gs3")
# ax.scatter(timevals,s3errvalss1,marker = 's',color='r',label = r"$\mathbf{S}^3$ gs1")
# ax.scatter(timevals,s3errvalss2,marker = 's',color='b',label = r"$\mathbf{S}^3$ gs2")
# ax.scatter(timevals,s3errvalss3,marker = 's',color='k',label = r"$\mathbf{S}^3$ gs3")
# ax.plot(timevals,timevals2,color='r',label = r"$dt^2$")
# ax.plot(timevals,timevals4,color='b',label = r"$dt^4$")
# ax.plot(timevals,timevals6,color='k',label = r"$dt^6$")

# # h3 
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-2}'
# # ylabels[1] = r'\textbf{-2}'
# # ylabels[2] = r'\textbf{0}'
# # ylabels[3] = r'\textbf{2}'
# # ylabels[4] = r'\textbf{4}'
# # ylabels[5] = r'\textbf{6}'
# # ylabels[6] = r'\textbf{8}'
# # ylabels[7] = r'\textbf{10}'


# # ax.set_xticklabels(xlabels)
# # ax.set_yticklabels(ylabels)
# ax.legend(fontsize="15",loc="upper left",ncol=2)
# plt.title(r'\textbf{Method Error Scaling for Rod Body in $\mathbf{S}^3$}', fontsize=20)
# # plt.title(r'\textbf{Method Error Scaling for Rod Body in $\mathbf{H}^3$}', fontsize=20)
# plt.ylabel(r'\textbf{Maximum Relative \\ \\ Error of Vertex Separation}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time Step Size dt}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)

# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim(1e-14,1e-0)

# plt.tight_layout()

# plt.show()


##### old stuff

# fig,ax=plt.subplots(1,1)

# # ax.scatter(timevals,h3errvalss1,marker = 's', color='r',label = "H3 Gauss s1")
# # ax.scatter(timevals,h3errvalss2,marker = 's', color='b',label = "H3 Gauss s2")
# # ax.scatter(timevals,h3errvalss3,marker = 's', color='k',label = "H3 Gauss s3")
# ax.scatter(timevals,s3errvalss1,marker = 'o',color='r',label = "S3 Gauss s1")
# ax.scatter(timevals,s3errvalss2,marker = 'o',color='b',label = "S3 Gauss s2")
# ax.scatter(timevals,s3errvalss3,marker = 'o',color='k',label = "S3 Gauss s3")
# ax.plot(timevals,timevals2,color='r',label = "h^2")
# ax.plot(timevals,timevals4,color='b',label = "h^4")
# ax.plot(timevals,timevals6,color='k',label = "h^6")

# ax.legend()
# ax.set_title('Vertex Separation Percent Error vs. Time')
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim(1e-14,1e-1)
# ax.set_xlabel('Time')
# ax.set_ylabel('Vertex Separation Percent Error')
# plt.show()



########################################################

######### Spring Triangle plots for H3 and S3 ##########

########################################################




# ####
# Spring Potential Triangle Separation plot for H3 and S3
# ####

# fig,ax=plt.subplots(1,1)

s3distdatatsp12gs3 = np.zeros(np.shape(datats3dt00001gs3)[0])
h3distdatatsp12gs3 = np.zeros(np.shape(datath3dt00001gs3)[0])

h3distdatatsp13gs3 = np.zeros(np.shape(datath3dt00001gs3)[0])
s3distdatatsp13gs3 = np.zeros(np.shape(datats3dt00001gs3)[0])

h3distdatatsp23gs3 = np.zeros(np.shape(datath3dt00001gs3)[0])
s3distdatatsp23gs3 = np.zeros(np.shape(datats3dt00001gs3)[0])

counter = 0
for a in range(np.shape(datath3dt00001gs3)[0]):
    h3distdatatsp12gs3[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][0:3]),rot2hyp(datath3dt00001gs3[a][3:6]))
    s3distdatatsp12gs3[counter] = r4dist(rot2r4(datats3dt00001gs3[a][0:3]),rot2r4(datats3dt00001gs3[a][3:6]))

    h3distdatatsp13gs3[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][0:3]),rot2hyp(datath3dt00001gs3[a][6:9]))
    s3distdatatsp13gs3[counter] = r4dist(rot2r4(datats3dt00001gs3[a][0:3]),rot2r4(datats3dt00001gs3[a][6:9]))

    h3distdatatsp23gs3[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][3:6]),rot2hyp(datath3dt00001gs3[a][6:9]))
    s3distdatatsp23gs3[counter] = r4dist(rot2r4(datats3dt00001gs3[a][3:6]),rot2r4(datats3dt00001gs3[a][6:9]))
    counter += 1

fig,ax=plt.subplots(1,1)
fig.canvas.draw()

ax.plot(t_arr00001,h3distdatatsp12gs3,'r',label = r"$\mathbf{H}^3$ s(t)")
ax.plot(t_arr00001,s3distdatatsp12gs3,'k',label = r"$\mathbf{S}^3$ s(t)")
# ax.plot(t_arr00001,h3distdatatsp13gs3,'r',label = r"$\mathbf{H}^3$ s13")
# ax.plot(t_arr00001,s3distdatatsp13gs3,'k',label = r"$\mathbf{S}^3$ s13")
ax.plot(t_arr00001,h3distdatatsp23gs3,'b',label = r"$\mathbf{H}^3$ 2l(t)")
ax.plot(t_arr00001,s3distdatatsp23gs3,'g',label = r"$\mathbf{S}^3$ 2l(t)")

# h3 
# ylabels = [item.get_text() for item in ax.get_yticklabels()]
# ylabels[0] = r'\textbf{-2}'
# ylabels[1] = r'\textbf{-2}'
# ylabels[2] = r'\textbf{0}'
# ylabels[3] = r'\textbf{2}'
# ylabels[4] = r'\textbf{4}'
# ylabels[5] = r'\textbf{6}'
# ylabels[6] = r'\textbf{8}'
# ylabels[7] = r'\textbf{10}'


# ax.set_xticklabels(xlabels)
# ax.set_yticklabels(ylabels)
ax.legend(fontsize="15",loc="upper right")
plt.title(r'\textbf{Elastic Triangle Body Dynamics}', fontsize=20)
plt.ylabel(r'\textbf{Vertex Separation}', fontsize=20, labelpad = 20)
plt.xlabel(r'\textbf{Time}', fontsize=20, labelpad = 20)

ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)


plt.tight_layout()

plt.show()

############## old stuff

# h3distdatatsp12gs3 = np.zeros(np.shape(datath3dt0001gs3)[0])
# s3distdatatsp12gs3 = np.zeros(np.shape(datats3dt0001gs3)[0])

# h3distdatatsp13gs3 = np.zeros(np.shape(datath3dt0001gs3)[0])
# s3distdatatsp13gs3 = np.zeros(np.shape(datats3dt0001gs3)[0])

# h3distdatatsp23gs3 = np.zeros(np.shape(datath3dt0001gs3)[0])
# s3distdatatsp23gs3 = np.zeros(np.shape(datats3dt0001gs3)[0])

# counter = 0
# for a in range(np.shape(datath3dt0001gs3)[0]):
#     h3distdatatsp12gs3[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][0:3]),rot2hyp(datath3dt0001gs3[a][3:6]))
#     s3distdatatsp12gs3[counter] = r4dist(rot2r4(datats3dt0001gs3[a][0:3]),rot2r4(datats3dt0001gs3[a][3:6]))

#     h3distdatatsp13gs3[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][0:3]),rot2hyp(datath3dt0001gs3[a][6:9]))
#     s3distdatatsp13gs3[counter] = r4dist(rot2r4(datats3dt0001gs3[a][0:3]),rot2r4(datats3dt0001gs3[a][6:9]))

#     h3distdatatsp23gs3[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][3:6]),rot2hyp(datath3dt0001gs3[a][6:9]))
#     s3distdatatsp23gs3[counter] = r4dist(rot2r4(datats3dt0001gs3[a][3:6]),rot2r4(datats3dt0001gs3[a][6:9]))
#     counter += 1

# ax.plot(t_arr0001,h3distdatatsp12gs3,'r',label = "Hyperbolic side 12")
# ax.plot(t_arr0001,s3distdatatsp12gs3,'k',label = "Spherical side 12")
# ax.plot(t_arr0001,h3distdatatsp13gs3,'r',label = "Hyperbolic side 13")
# ax.plot(t_arr0001,s3distdatatsp13gs3,'k',label = "Spherical side 13")
# ax.plot(t_arr0001,h3distdatatsp23gs3,'r',label = "Hyperbolic side 23")
# ax.plot(t_arr0001,s3distdatatsp23gs3,'k',label = "Spherical side 23")
# ax.legend()
# ax.set_title('Vertex Separation vs. Time')
# ax.set_xlabel('Time')
# ax.set_ylabel('Vertex Separation')
# plt.show()


######
# Spring Potential Triangle Separation Error versus midpooint distance for H3 and S3
######

# ### Gauss s1

# h3disterrdatatsp12gs1dt1raw     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3disterrdatatsp12gs1dt01raw    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3disterrdatatsp12gs1dt001raw   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3disterrdatatsp12gs1dt0001raw  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3disterrdatatsp12gs1dt00001raw = np.zeros(np.shape(datath3dt00001gs1)[0])

# h3distsp1gs1dt1     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3distsp1gs1dt01    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3distsp1gs1dt001   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3distsp1gs1dt0001  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3distsp1gs1dt00001 = np.zeros(np.shape(datath3dt00001gs1)[0])

# h3disterrdatatsp13gs1dt1raw     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3disterrdatatsp13gs1dt01raw    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3disterrdatatsp13gs1dt001raw   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3disterrdatatsp13gs1dt0001raw  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3disterrdatatsp13gs1dt00001raw = np.zeros(np.shape(datath3dt00001gs1)[0])

# h3distsp2gs1dt1     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3distsp2gs1dt01    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3distsp2gs1dt001   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3distsp2gs1dt0001  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3distsp2gs1dt00001 = np.zeros(np.shape(datath3dt00001gs1)[0])

# h3disterrdatatsp23gs1dt1raw     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3disterrdatatsp23gs1dt01raw    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3disterrdatatsp23gs1dt001raw   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3disterrdatatsp23gs1dt0001raw  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3disterrdatatsp23gs1dt00001raw = np.zeros(np.shape(datath3dt00001gs1)[0])

# h3distsp3gs1dt1     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3distsp3gs1dt01    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3distsp3gs1dt001   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3distsp3gs1dt0001  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3distsp3gs1dt00001 = np.zeros(np.shape(datath3dt00001gs1)[0])

# s3disterrdatatsp12gs1dt1raw     = np.zeros(np.shape(datats3dt1gs1)[0])
# s3disterrdatatsp12gs1dt01raw    = np.zeros(np.shape(datats3dt01gs1)[0])
# s3disterrdatatsp12gs1dt001raw   = np.zeros(np.shape(datats3dt001gs1)[0])
# s3disterrdatatsp12gs1dt0001raw  = np.zeros(np.shape(datats3dt0001gs1)[0])
# s3disterrdatatsp12gs1dt00001raw = np.zeros(np.shape(datats3dt00001gs1)[0])

# s3disterrdatatsp13gs1dt1raw     = np.zeros(np.shape(datats3dt1gs1)[0])
# s3disterrdatatsp13gs1dt01raw    = np.zeros(np.shape(datats3dt01gs1)[0])
# s3disterrdatatsp13gs1dt001raw   = np.zeros(np.shape(datats3dt001gs1)[0])
# s3disterrdatatsp13gs1dt0001raw  = np.zeros(np.shape(datats3dt0001gs1)[0])
# s3disterrdatatsp13gs1dt00001raw = np.zeros(np.shape(datats3dt00001gs1)[0])

# s3disterrdatatsp23gs1dt1raw     = np.zeros(np.shape(datats3dt1gs1)[0])
# s3disterrdatatsp23gs1dt01raw    = np.zeros(np.shape(datats3dt01gs1)[0])
# s3disterrdatatsp23gs1dt001raw   = np.zeros(np.shape(datats3dt001gs1)[0])
# s3disterrdatatsp23gs1dt0001raw  = np.zeros(np.shape(datats3dt0001gs1)[0])
# s3disterrdatatsp23gs1dt00001raw = np.zeros(np.shape(datats3dt00001gs1)[0])

# ### Gauss s2

# h3disterrdatatsp12gs2dt1raw     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3disterrdatatsp12gs2dt01raw    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3disterrdatatsp12gs2dt001raw   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3disterrdatatsp12gs2dt0001raw  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3disterrdatatsp12gs2dt00001raw = np.zeros(np.shape(datath3dt00001gs2)[0])

# h3distsp1gs2dt1     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3distsp1gs2dt01    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3distsp1gs2dt001   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3distsp1gs2dt0001  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3distsp1gs2dt00001 = np.zeros(np.shape(datath3dt00001gs2)[0])

# h3disterrdatatsp13gs2dt1raw     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3disterrdatatsp13gs2dt01raw    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3disterrdatatsp13gs2dt001raw   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3disterrdatatsp13gs2dt0001raw  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3disterrdatatsp13gs2dt00001raw = np.zeros(np.shape(datath3dt00001gs2)[0])

# h3distsp2gs2dt1     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3distsp2gs2dt01    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3distsp2gs2dt001   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3distsp2gs2dt0001  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3distsp2gs2dt00001 = np.zeros(np.shape(datath3dt00001gs2)[0])

# h3disterrdatatsp23gs2dt1raw     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3disterrdatatsp23gs2dt01raw    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3disterrdatatsp23gs2dt001raw   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3disterrdatatsp23gs2dt0001raw  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3disterrdatatsp23gs2dt00001raw = np.zeros(np.shape(datath3dt00001gs2)[0])

# h3distsp3gs2dt1     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3distsp3gs2dt01    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3distsp3gs2dt001   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3distsp3gs2dt0001  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3distsp3gs2dt00001 = np.zeros(np.shape(datath3dt00001gs2)[0])

# s3disterrdatatsp12gs2dt1raw     = np.zeros(np.shape(datats3dt1gs2)[0])
# s3disterrdatatsp12gs2dt01raw    = np.zeros(np.shape(datats3dt01gs2)[0])
# s3disterrdatatsp12gs2dt001raw   = np.zeros(np.shape(datats3dt001gs2)[0])
# s3disterrdatatsp12gs2dt0001raw  = np.zeros(np.shape(datats3dt0001gs2)[0])
# s3disterrdatatsp12gs2dt00001raw = np.zeros(np.shape(datats3dt00001gs2)[0])

# s3disterrdatatsp13gs2dt1raw     = np.zeros(np.shape(datats3dt1gs2)[0])
# s3disterrdatatsp13gs2dt01raw    = np.zeros(np.shape(datats3dt01gs2)[0])
# s3disterrdatatsp13gs2dt001raw   = np.zeros(np.shape(datats3dt001gs2)[0])
# s3disterrdatatsp13gs2dt0001raw  = np.zeros(np.shape(datats3dt0001gs2)[0])
# s3disterrdatatsp13gs2dt00001raw = np.zeros(np.shape(datats3dt00001gs2)[0])

# s3disterrdatatsp23gs2dt1raw     = np.zeros(np.shape(datats3dt1gs2)[0])
# s3disterrdatatsp23gs2dt01raw    = np.zeros(np.shape(datats3dt01gs2)[0])
# s3disterrdatatsp23gs2dt001raw   = np.zeros(np.shape(datats3dt001gs2)[0])
# s3disterrdatatsp23gs2dt0001raw  = np.zeros(np.shape(datats3dt0001gs2)[0])
# s3disterrdatatsp23gs2dt00001raw = np.zeros(np.shape(datats3dt00001gs2)[0])

# ### Gauss s3

# h3disterrdatatsp12gs3dt1raw     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3disterrdatatsp12gs3dt01raw    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3disterrdatatsp12gs3dt001raw   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3disterrdatatsp12gs3dt0001raw  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3disterrdatatsp12gs3dt00001raw = np.zeros(np.shape(datath3dt00001gs3)[0])

# h3distsp1gs3dt1     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3distsp1gs3dt01    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3distsp1gs3dt001   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3distsp1gs3dt0001  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3distsp1gs3dt00001 = np.zeros(np.shape(datath3dt00001gs3)[0])

# h3disterrdatatsp13gs3dt1raw     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3disterrdatatsp13gs3dt01raw    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3disterrdatatsp13gs3dt001raw   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3disterrdatatsp13gs3dt0001raw  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3disterrdatatsp13gs3dt00001raw = np.zeros(np.shape(datath3dt00001gs3)[0])

# h3distsp2gs3dt1     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3distsp2gs3dt01    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3distsp2gs3dt001   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3distsp2gs3dt0001  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3distsp2gs3dt00001 = np.zeros(np.shape(datath3dt00001gs3)[0])

# h3disterrdatatsp23gs3dt1raw     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3disterrdatatsp23gs3dt01raw    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3disterrdatatsp23gs3dt001raw   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3disterrdatatsp23gs3dt0001raw  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3disterrdatatsp23gs3dt00001raw = np.zeros(np.shape(datath3dt00001gs3)[0])

# h3distsp3gs3dt1     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3distsp3gs3dt01    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3distsp3gs3dt001   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3distsp3gs3dt0001  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3distsp3gs3dt00001 = np.zeros(np.shape(datath3dt00001gs3)[0])

# s3disterrdatatsp12gs3dt1raw     = np.zeros(np.shape(datats3dt1gs3)[0])
# s3disterrdatatsp12gs3dt01raw    = np.zeros(np.shape(datats3dt01gs3)[0])
# s3disterrdatatsp12gs3dt001raw   = np.zeros(np.shape(datats3dt001gs3)[0])
# s3disterrdatatsp12gs3dt0001raw  = np.zeros(np.shape(datats3dt0001gs3)[0])
# s3disterrdatatsp12gs3dt00001raw = np.zeros(np.shape(datats3dt00001gs3)[0])

# s3disterrdatatsp13gs3dt1raw     = np.zeros(np.shape(datats3dt1gs3)[0])
# s3disterrdatatsp13gs3dt01raw    = np.zeros(np.shape(datats3dt01gs3)[0])
# s3disterrdatatsp13gs3dt001raw   = np.zeros(np.shape(datats3dt001gs3)[0])
# s3disterrdatatsp13gs3dt0001raw  = np.zeros(np.shape(datats3dt0001gs3)[0])
# s3disterrdatatsp13gs3dt00001raw = np.zeros(np.shape(datats3dt00001gs3)[0])

# s3disterrdatatsp23gs3dt1raw     = np.zeros(np.shape(datats3dt1gs3)[0])
# s3disterrdatatsp23gs3dt01raw    = np.zeros(np.shape(datats3dt01gs3)[0])
# s3disterrdatatsp23gs3dt001raw   = np.zeros(np.shape(datats3dt001gs3)[0])
# s3disterrdatatsp23gs3dt0001raw  = np.zeros(np.shape(datats3dt0001gs3)[0])
# s3disterrdatatsp23gs3dt00001raw = np.zeros(np.shape(datats3dt00001gs3)[0])


# counter = 0
# print("Collating dt = .1 Data")
# for a in range(np.shape(datath3dt1gs1)[0]):
#     h3disterrdatatsp12gs1dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs1[a][0:3]),rot2hyp(datath3dt1gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs2[a][0:3]),rot2hyp(datath3dt1gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs3[a][0:3]),rot2hyp(datath3dt1gs3[a][3:6]))

#     h3distsp1gs1dt1[counter] = datath3dt1gs1[a][0]
#     h3distsp1gs2dt1[counter] = datath3dt1gs2[a][0]
#     h3distsp1gs3dt1[counter] = datath3dt1gs3[a][0]

#     h3disterrdatatsp13gs1dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs1[a][0:3]),rot2hyp(datath3dt1gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs2[a][0:3]),rot2hyp(datath3dt1gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs3[a][0:3]),rot2hyp(datath3dt1gs3[a][6:9]))

#     h3distsp2gs1dt1[counter] = arccosh(cosh(datath3dt1gs1[a][3])/cosh(.5*h3disterrdatatsp23gs1dt1raw[counter]))
#     h3distsp2gs2dt1[counter] = arccosh(cosh(datath3dt1gs2[a][3])/cosh(.5*h3disterrdatatsp23gs2dt1raw[counter]))
#     h3distsp2gs3dt1[counter] = arccosh(cosh(datath3dt1gs3[a][3])/cosh(.5*h3disterrdatatsp23gs3dt1raw[counter]))

#     h3disterrdatatsp23gs1dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs1[a][3:6]),rot2hyp(datath3dt1gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs2[a][3:6]),rot2hyp(datath3dt1gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs3[a][3:6]),rot2hyp(datath3dt1gs3[a][6:9]))

#     h3distsp3gs1dt1[counter] = arccosh(cosh(datath3dt1gs1[a][6])/cosh(.5*h3disterrdatatsp23gs1dt1raw[counter]))
#     h3distsp3gs2dt1[counter] = arccosh(cosh(datath3dt1gs2[a][6])/cosh(.5*h3disterrdatatsp23gs2dt1raw[counter]))
#     h3distsp3gs3dt1[counter] = arccosh(cosh(datath3dt1gs3[a][6])/cosh(.5*h3disterrdatatsp23gs3dt1raw[counter]))

#     s3disterrdatatsp12gs1dt1raw[counter] = r4dist(rot2r4(datats3dt1gs1[a][0:3]),rot2r4(datats3dt1gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt1raw[counter] = r4dist(rot2r4(datats3dt1gs2[a][0:3]),rot2r4(datats3dt1gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt1raw[counter] = r4dist(rot2r4(datats3dt1gs3[a][0:3]),rot2r4(datats3dt1gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt1raw[counter] = r4dist(rot2r4(datats3dt1gs1[a][0:3]),rot2r4(datats3dt1gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt1raw[counter] = r4dist(rot2r4(datats3dt1gs2[a][0:3]),rot2r4(datats3dt1gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt1raw[counter] = r4dist(rot2r4(datats3dt1gs3[a][0:3]),rot2r4(datats3dt1gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt1raw[counter] = r4dist(rot2r4(datats3dt1gs1[a][3:6]),rot2r4(datats3dt1gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt1raw[counter] = r4dist(rot2r4(datats3dt1gs2[a][3:6]),rot2r4(datats3dt1gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt1raw[counter] = r4dist(rot2r4(datats3dt1gs3[a][3:6]),rot2r4(datats3dt1gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp12gs1dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp12gs2dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp12gs2dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp12gs3dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp12gs3dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))

# h3disterrdatatsp13gs1dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp13gs1dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp13gs2dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp13gs2dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp13gs3dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp13gs3dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))

# h3disterrdatatsp23gs1dt1 = (2.*datath3gtdt00001gs3[::10000,0] - h3disterrdatatsp23gs1dt1raw)/(2.*datath3gtdt00001gs3[::10000,0])
# h3disterrdatatsp23gs2dt1 = (2.*datath3gtdt00001gs3[::10000,0] - h3disterrdatatsp23gs2dt1raw)/(2.*datath3gtdt00001gs3[::10000,0])
# h3disterrdatatsp23gs3dt1 = (2.*datath3gtdt00001gs3[::10000,0] - h3disterrdatatsp23gs3dt1raw)/(2.*datath3gtdt00001gs3[::10000,0])

# s3disterrdatatsp12gs1dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp12gs1dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp12gs2dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp12gs2dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp12gs3dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp12gs3dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))

# s3disterrdatatsp13gs1dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp13gs1dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp13gs2dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp13gs2dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp13gs3dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp13gs3dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))

# s3disterrdatatsp23gs1dt1 = (2.*datats3gtdt00001gs3[::10000,0] - s3disterrdatatsp23gs1dt1raw)/(2.*datats3gtdt00001gs3[::10000,0])
# s3disterrdatatsp23gs2dt1 = (2.*datats3gtdt00001gs3[::10000,0] - s3disterrdatatsp23gs2dt1raw)/(2.*datats3gtdt00001gs3[::10000,0])
# s3disterrdatatsp23gs3dt1 = (2.*datats3gtdt00001gs3[::10000,0] - s3disterrdatatsp23gs3dt1raw)/(2.*datats3gtdt00001gs3[::10000,0])

# counter = 0
# print("Collating dt = .01 Data")
# for a in range(np.shape(datath3dt01gs1)[0]):
#     h3disterrdatatsp12gs1dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs1[a][0:3]),rot2hyp(datath3dt01gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs2[a][0:3]),rot2hyp(datath3dt01gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs3[a][0:3]),rot2hyp(datath3dt01gs3[a][3:6]))

#     h3distsp1gs1dt01[counter] = datath3dt01gs1[a][0]
#     h3distsp1gs2dt01[counter] = datath3dt01gs2[a][0]
#     h3distsp1gs3dt01[counter] = datath3dt01gs3[a][0]

#     h3disterrdatatsp13gs1dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs1[a][0:3]),rot2hyp(datath3dt01gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs2[a][0:3]),rot2hyp(datath3dt01gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs3[a][0:3]),rot2hyp(datath3dt01gs3[a][6:9]))

#     h3distsp2gs1dt01[counter] = arccosh(cosh(datath3dt01gs1[a][3])/cosh(.5*h3disterrdatatsp23gs1dt01raw[counter]))
#     h3distsp2gs2dt01[counter] = arccosh(cosh(datath3dt01gs2[a][3])/cosh(.5*h3disterrdatatsp23gs2dt01raw[counter]))
#     h3distsp2gs3dt01[counter] = arccosh(cosh(datath3dt01gs3[a][3])/cosh(.5*h3disterrdatatsp23gs3dt01raw[counter]))

#     h3disterrdatatsp23gs1dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs1[a][3:6]),rot2hyp(datath3dt01gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs2[a][3:6]),rot2hyp(datath3dt01gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs3[a][3:6]),rot2hyp(datath3dt01gs3[a][6:9]))

#     h3distsp2gs1dt01[counter] = arccosh(cosh(datath3dt01gs1[a][6])/cosh(.5*h3disterrdatatsp23gs1dt01raw[counter]))
#     h3distsp2gs2dt01[counter] = arccosh(cosh(datath3dt01gs2[a][6])/cosh(.5*h3disterrdatatsp23gs2dt01raw[counter]))
#     h3distsp2gs3dt01[counter] = arccosh(cosh(datath3dt01gs3[a][6])/cosh(.5*h3disterrdatatsp23gs3dt01raw[counter]))

#     s3disterrdatatsp12gs1dt01raw[counter] = r4dist(rot2r4(datats3dt01gs1[a][0:3]),rot2r4(datats3dt01gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt01raw[counter] = r4dist(rot2r4(datats3dt01gs2[a][0:3]),rot2r4(datats3dt01gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt01raw[counter] = r4dist(rot2r4(datats3dt01gs3[a][0:3]),rot2r4(datats3dt01gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt01raw[counter] = r4dist(rot2r4(datats3dt01gs1[a][0:3]),rot2r4(datats3dt01gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt01raw[counter] = r4dist(rot2r4(datats3dt01gs2[a][0:3]),rot2r4(datats3dt01gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt01raw[counter] = r4dist(rot2r4(datats3dt01gs3[a][0:3]),rot2r4(datats3dt01gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt01raw[counter] = r4dist(rot2r4(datats3dt01gs1[a][3:6]),rot2r4(datats3dt01gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt01raw[counter] = r4dist(rot2r4(datats3dt01gs2[a][3:6]),rot2r4(datats3dt01gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt01raw[counter] = r4dist(rot2r4(datats3dt01gs3[a][3:6]),rot2r4(datats3dt01gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp12gs1dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp12gs2dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp12gs2dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp12gs3dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp12gs3dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))

# h3disterrdatatsp13gs1dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp13gs1dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp13gs2dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp13gs2dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp13gs3dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp13gs3dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))

# h3disterrdatatsp23gs1dt01 = (2.*datath3gtdt00001gs3[::1000,0] - h3disterrdatatsp23gs1dt01raw)/(2.*datath3gtdt00001gs3[::1000,0])
# h3disterrdatatsp23gs2dt01 = (2.*datath3gtdt00001gs3[::1000,0] - h3disterrdatatsp23gs2dt01raw)/(2.*datath3gtdt00001gs3[::1000,0])
# h3disterrdatatsp23gs3dt01 = (2.*datath3gtdt00001gs3[::1000,0] - h3disterrdatatsp23gs3dt01raw)/(2.*datath3gtdt00001gs3[::1000,0])

# s3disterrdatatsp12gs1dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp12gs1dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp12gs2dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp12gs2dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp12gs3dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp12gs3dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))

# s3disterrdatatsp13gs1dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp13gs1dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp13gs2dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp13gs2dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp13gs3dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp13gs3dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))

# s3disterrdatatsp23gs1dt01 = (2.*datats3gtdt00001gs3[::1000,0] - s3disterrdatatsp23gs1dt01raw)/(2.*datats3gtdt00001gs3[::1000,0])
# s3disterrdatatsp23gs2dt01 = (2.*datats3gtdt00001gs3[::1000,0] - s3disterrdatatsp23gs2dt01raw)/(2.*datats3gtdt00001gs3[::1000,0])
# s3disterrdatatsp23gs3dt01 = (2.*datats3gtdt00001gs3[::1000,0] - s3disterrdatatsp23gs3dt01raw)/(2.*datats3gtdt00001gs3[::1000,0])

# counter = 0
# print("Collating dt = .001 Data")
# for a in range(np.shape(datath3dt001gs1)[0]):
#     h3disterrdatatsp12gs1dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs1[a][0:3]),rot2hyp(datath3dt001gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs2[a][0:3]),rot2hyp(datath3dt001gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs3[a][0:3]),rot2hyp(datath3dt001gs3[a][3:6]))

#     h3distsp1gs1dt001[counter] = datath3dt001gs1[a][0]
#     h3distsp1gs2dt001[counter] = datath3dt001gs2[a][0]
#     h3distsp1gs3dt001[counter] = datath3dt001gs3[a][0]

#     h3disterrdatatsp13gs1dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs1[a][0:3]),rot2hyp(datath3dt001gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs2[a][0:3]),rot2hyp(datath3dt001gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs3[a][0:3]),rot2hyp(datath3dt001gs3[a][6:9]))

#     h3distsp2gs1dt001[counter] = arccosh(cosh(datath3dt001gs1[a][3])/cosh(.5*h3disterrdatatsp23gs1dt001raw[counter]))
#     h3distsp2gs2dt001[counter] = arccosh(cosh(datath3dt001gs2[a][3])/cosh(.5*h3disterrdatatsp23gs2dt001raw[counter]))
#     h3distsp2gs3dt001[counter] = arccosh(cosh(datath3dt001gs3[a][3])/cosh(.5*h3disterrdatatsp23gs3dt001raw[counter]))

#     h3disterrdatatsp23gs1dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs1[a][3:6]),rot2hyp(datath3dt001gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs2[a][3:6]),rot2hyp(datath3dt001gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs3[a][3:6]),rot2hyp(datath3dt001gs3[a][6:9]))

#     h3distsp2gs1dt001[counter] = arccosh(cosh(datath3dt001gs1[a][6])/cosh(.5*h3disterrdatatsp23gs1dt001raw[counter]))
#     h3distsp2gs2dt001[counter] = arccosh(cosh(datath3dt001gs2[a][6])/cosh(.5*h3disterrdatatsp23gs2dt001raw[counter]))
#     h3distsp2gs3dt001[counter] = arccosh(cosh(datath3dt001gs3[a][6])/cosh(.5*h3disterrdatatsp23gs3dt001raw[counter]))

#     s3disterrdatatsp12gs1dt001raw[counter] = r4dist(rot2r4(datats3dt001gs1[a][0:3]),rot2r4(datats3dt001gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt001raw[counter] = r4dist(rot2r4(datats3dt001gs2[a][0:3]),rot2r4(datats3dt001gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt001raw[counter] = r4dist(rot2r4(datats3dt001gs3[a][0:3]),rot2r4(datats3dt001gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt001raw[counter] = r4dist(rot2r4(datats3dt001gs1[a][0:3]),rot2r4(datats3dt001gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt001raw[counter] = r4dist(rot2r4(datats3dt001gs2[a][0:3]),rot2r4(datats3dt001gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt001raw[counter] = r4dist(rot2r4(datats3dt001gs3[a][0:3]),rot2r4(datats3dt001gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt001raw[counter] = r4dist(rot2r4(datats3dt001gs1[a][3:6]),rot2r4(datats3dt001gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt001raw[counter] = r4dist(rot2r4(datats3dt001gs2[a][3:6]),rot2r4(datats3dt001gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt001raw[counter] = r4dist(rot2r4(datats3dt001gs3[a][3:6]),rot2r4(datats3dt001gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp12gs1dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp12gs2dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp12gs2dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp12gs3dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp12gs3dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))

# h3disterrdatatsp13gs1dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp13gs1dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp13gs2dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp13gs2dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp13gs3dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp13gs3dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))

# h3disterrdatatsp23gs1dt001 = (2.*datath3gtdt00001gs3[::100,0] - h3disterrdatatsp23gs1dt001raw)/(2.*datath3gtdt00001gs3[::100,0])
# h3disterrdatatsp23gs2dt001 = (2.*datath3gtdt00001gs3[::100,0] - h3disterrdatatsp23gs2dt001raw)/(2.*datath3gtdt00001gs3[::100,0])
# h3disterrdatatsp23gs3dt001 = (2.*datath3gtdt00001gs3[::100,0] - h3disterrdatatsp23gs3dt001raw)/(2.*datath3gtdt00001gs3[::100,0])

# s3disterrdatatsp12gs1dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp12gs1dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp12gs2dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp12gs2dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp12gs3dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp12gs3dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))

# s3disterrdatatsp13gs1dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp13gs1dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp13gs2dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp13gs2dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp13gs3dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp13gs3dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))

# s3disterrdatatsp23gs1dt001 = (2.*datats3gtdt00001gs3[::100,0] - s3disterrdatatsp23gs1dt001raw)/(2.*datats3gtdt00001gs3[::100,0])
# s3disterrdatatsp23gs2dt001 = (2.*datats3gtdt00001gs3[::100,0] - s3disterrdatatsp23gs2dt001raw)/(2.*datats3gtdt00001gs3[::100,0])
# s3disterrdatatsp23gs3dt001 = (2.*datats3gtdt00001gs3[::100,0] - s3disterrdatatsp23gs3dt001raw)/(2.*datats3gtdt00001gs3[::100,0])

# counter = 0
# print("Collating dt = .0001 Data")
# for a in range(np.shape(datath3dt0001gs1)[0]):
#     h3disterrdatatsp12gs1dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs1[a][0:3]),rot2hyp(datath3dt0001gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs2[a][0:3]),rot2hyp(datath3dt0001gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][0:3]),rot2hyp(datath3dt0001gs3[a][3:6]))

#     h3distsp1gs1dt0001[counter] = datath3dt0001gs1[a][0]
#     h3distsp1gs2dt0001[counter] = datath3dt0001gs2[a][0]
#     h3distsp1gs3dt0001[counter] = datath3dt0001gs3[a][0]

#     h3disterrdatatsp13gs1dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs1[a][0:3]),rot2hyp(datath3dt0001gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs2[a][0:3]),rot2hyp(datath3dt0001gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][0:3]),rot2hyp(datath3dt0001gs3[a][6:9]))

#     h3distsp2gs1dt0001[counter] = arccosh(cosh(datath3dt0001gs1[a][3])/cosh(.5*h3disterrdatatsp23gs1dt0001raw[counter]))
#     h3distsp2gs2dt0001[counter] = arccosh(cosh(datath3dt0001gs2[a][3])/cosh(.5*h3disterrdatatsp23gs2dt0001raw[counter]))
#     h3distsp2gs3dt0001[counter] = arccosh(cosh(datath3dt0001gs3[a][3])/cosh(.5*h3disterrdatatsp23gs3dt0001raw[counter]))

#     h3disterrdatatsp23gs1dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs1[a][3:6]),rot2hyp(datath3dt0001gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs2[a][3:6]),rot2hyp(datath3dt0001gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][3:6]),rot2hyp(datath3dt0001gs3[a][6:9]))

#     h3distsp2gs1dt0001[counter] = arccosh(cosh(datath3dt0001gs1[a][6])/cosh(.5*h3disterrdatatsp23gs1dt0001raw[counter]))
#     h3distsp2gs2dt0001[counter] = arccosh(cosh(datath3dt0001gs2[a][6])/cosh(.5*h3disterrdatatsp23gs2dt0001raw[counter]))
#     h3distsp2gs3dt0001[counter] = arccosh(cosh(datath3dt0001gs3[a][6])/cosh(.5*h3disterrdatatsp23gs3dt0001raw[counter]))

#     s3disterrdatatsp12gs1dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs1[a][0:3]),rot2r4(datats3dt0001gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs2[a][0:3]),rot2r4(datats3dt0001gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs3[a][0:3]),rot2r4(datats3dt0001gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs1[a][0:3]),rot2r4(datats3dt0001gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs2[a][0:3]),rot2r4(datats3dt0001gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs3[a][0:3]),rot2r4(datats3dt0001gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs1[a][3:6]),rot2r4(datats3dt0001gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs2[a][3:6]),rot2r4(datats3dt0001gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs3[a][3:6]),rot2r4(datats3dt0001gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp12gs1dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp12gs2dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp12gs2dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp12gs3dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp12gs3dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))

# h3disterrdatatsp13gs1dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp13gs1dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp13gs2dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp13gs2dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp13gs3dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp13gs3dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))

# h3disterrdatatsp23gs1dt0001 = (2.*datath3gtdt00001gs3[::10,0] - h3disterrdatatsp23gs1dt0001raw)/(2.*datath3gtdt00001gs3[::10,0])
# h3disterrdatatsp23gs2dt0001 = (2.*datath3gtdt00001gs3[::10,0] - h3disterrdatatsp23gs2dt0001raw)/(2.*datath3gtdt00001gs3[::10,0])
# h3disterrdatatsp23gs3dt0001 = (2.*datath3gtdt00001gs3[::10,0] - h3disterrdatatsp23gs3dt0001raw)/(2.*datath3gtdt00001gs3[::10,0])

# s3disterrdatatsp12gs1dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp12gs1dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp12gs2dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp12gs2dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp12gs3dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp12gs3dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))

# s3disterrdatatsp13gs1dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp13gs1dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp13gs2dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp13gs2dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp13gs3dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp13gs3dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))

# s3disterrdatatsp23gs1dt0001 = (2.*datats3gtdt00001gs3[::10,0] - s3disterrdatatsp23gs1dt0001raw)/(2.*datats3gtdt00001gs3[::10,0])
# s3disterrdatatsp23gs2dt0001 = (2.*datats3gtdt00001gs3[::10,0] - s3disterrdatatsp23gs2dt0001raw)/(2.*datats3gtdt00001gs3[::10,0])
# s3disterrdatatsp23gs3dt0001 = (2.*datats3gtdt00001gs3[::10,0] - s3disterrdatatsp23gs3dt0001raw)/(2.*datats3gtdt00001gs3[::10,0])

# counter = 0
# print("Collating dt = .00001 Data")
# for a in range(np.shape(datath3dt00001gs1)[0]):
#     h3disterrdatatsp12gs1dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs1[a][0:3]),rot2hyp(datath3dt00001gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs2[a][0:3]),rot2hyp(datath3dt00001gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][0:3]),rot2hyp(datath3dt00001gs3[a][3:6]))

#     h3distsp1gs1dt00001[counter] = datath3dt00001gs1[a][0]
#     h3distsp1gs2dt00001[counter] = datath3dt00001gs2[a][0]
#     h3distsp1gs3dt00001[counter] = datath3dt00001gs3[a][0]

#     h3disterrdatatsp13gs1dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs1[a][0:3]),rot2hyp(datath3dt00001gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs2[a][0:3]),rot2hyp(datath3dt00001gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][0:3]),rot2hyp(datath3dt00001gs3[a][6:9]))

#     h3distsp2gs1dt00001[counter] = arccosh(cosh(datath3dt00001gs1[a][3])/cosh(.5*h3disterrdatatsp23gs1dt00001raw[counter]))
#     h3distsp2gs2dt00001[counter] = arccosh(cosh(datath3dt00001gs2[a][3])/cosh(.5*h3disterrdatatsp23gs2dt00001raw[counter]))
#     h3distsp2gs3dt00001[counter] = arccosh(cosh(datath3dt00001gs3[a][3])/cosh(.5*h3disterrdatatsp23gs3dt00001raw[counter]))

#     h3disterrdatatsp23gs1dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs1[a][3:6]),rot2hyp(datath3dt00001gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs2[a][3:6]),rot2hyp(datath3dt00001gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][3:6]),rot2hyp(datath3dt00001gs3[a][6:9]))

#     h3distsp2gs1dt00001[counter] = arccosh(cosh(datath3dt00001gs1[a][6])/cosh(.5*h3disterrdatatsp23gs1dt00001raw[counter]))
#     h3distsp2gs2dt00001[counter] = arccosh(cosh(datath3dt00001gs2[a][6])/cosh(.5*h3disterrdatatsp23gs2dt00001raw[counter]))
#     h3distsp2gs3dt00001[counter] = arccosh(cosh(datath3dt00001gs3[a][6])/cosh(.5*h3disterrdatatsp23gs3dt00001raw[counter]))

#     s3disterrdatatsp12gs1dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs1[a][0:3]),rot2r4(datats3dt00001gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs2[a][0:3]),rot2r4(datats3dt00001gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs3[a][0:3]),rot2r4(datats3dt00001gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs1[a][0:3]),rot2r4(datats3dt00001gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs2[a][0:3]),rot2r4(datats3dt00001gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs3[a][0:3]),rot2r4(datats3dt00001gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs1[a][3:6]),rot2r4(datats3dt00001gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs2[a][3:6]),rot2r4(datats3dt00001gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs3[a][3:6]),rot2r4(datats3dt00001gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp12gs1dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp12gs2dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp12gs2dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp12gs3dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp12gs3dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))

# h3disterrdatatsp13gs1dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp13gs1dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp13gs2dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp13gs2dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp13gs3dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp13gs3dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))

# h3disterrdatatsp23gs1dt00001 = (2.*datath3gtdt00001gs3[::1,0] - h3disterrdatatsp23gs1dt00001raw)/(2.*datath3gtdt00001gs3[::1,0])
# h3disterrdatatsp23gs2dt00001 = (2.*datath3gtdt00001gs3[::1,0] - h3disterrdatatsp23gs2dt00001raw)/(2.*datath3gtdt00001gs3[::1,0])
# h3disterrdatatsp23gs3dt00001 = (2.*datath3gtdt00001gs3[::1,0] - h3disterrdatatsp23gs3dt00001raw)/(2.*datath3gtdt00001gs3[::1,0])

# s3disterrdatatsp12gs1dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp12gs1dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp12gs2dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp12gs2dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp12gs3dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp12gs3dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))

# s3disterrdatatsp13gs1dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp13gs1dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp13gs2dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp13gs2dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp13gs3dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp13gs3dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))

# s3disterrdatatsp23gs1dt00001 = (2.*datats3gtdt00001gs3[::1,0] - s3disterrdatatsp23gs1dt00001raw)/(2.*datats3gtdt00001gs3[::1,0])
# s3disterrdatatsp23gs2dt00001 = (2.*datats3gtdt00001gs3[::1,0] - s3disterrdatatsp23gs2dt00001raw)/(2.*datats3gtdt00001gs3[::1,0])
# s3disterrdatatsp23gs3dt00001 = (2.*datats3gtdt00001gs3[::1,0] - s3disterrdatatsp23gs3dt00001raw)/(2.*datats3gtdt00001gs3[::1,0])


# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# # ax.plot(t_arr00001,h3disterrdatatsp12gs1dt00001, color='r',label = r"$\mathbf{H}^3$ gs1")
# # ax.plot(t_arr00001,h3disterrdatatsp12gs2dt00001, color='b',label = r"$\mathbf{H}^3$ gs2")
# # ax.plot(t_arr00001,h3disterrdatatsp12gs3dt00001, color='k',label = r"$\mathbf{H}^3$ gs3")

# # ax.plot(t_arr00001,h3disterrdatatsp13gs1dt00001, color='r',label = r"$\mathbf{H}^3$ gs1")
# # ax.plot(t_arr00001,h3disterrdatatsp13gs2dt00001, color='b',label = r"$\mathbf{H}^3$ gs2")
# # ax.plot(t_arr00001,h3disterrdatatsp13gs3dt00001, color='k',label = r"$\mathbf{H}^3$ gs3")

# # ax.plot(t_arr00001,h3disterrdatatsp23gs1dt00001, color='r',label = r"$\mathbf{H}^3$ gs1")
# # ax.plot(t_arr00001,h3disterrdatatsp23gs2dt00001, color='b',label = r"$\mathbf{H}^3$ gs2")
# # ax.plot(t_arr00001,h3disterrdatatsp23gs3dt00001, color='k',label = r"$\mathbf{H}^3$ gs3")

# # ax.plot(t_arr00001,s3disterrdatatsp12gs1dt00001, color='r',label = r"$\mathbf{S}^3$ gs1")
# # ax.plot(t_arr00001,s3disterrdatatsp12gs2dt00001, color='b',label = r"$\mathbf{S}^3$ gs2")
# # ax.plot(t_arr00001,s3disterrdatatsp12gs3dt00001, color='k',label = r"$\mathbf{S}^3$ gs3")

# # ax.plot(t_arr00001,s3disterrdatatsp13gs1dt00001, color='r',label = r"$\mathbf{S}^3$ gs1")
# # ax.plot(t_arr00001,s3disterrdatatsp13gs2dt00001, color='b',label = r"$\mathbf{S}^3$ gs2")
# # ax.plot(t_arr00001,s3disterrdatatsp13gs3dt00001, color='k',label = r"$\mathbf{S}^3$ gs3")

# ax.plot(t_arr00001,s3disterrdatatsp23gs1dt00001, color='r',label = r"$\mathbf{S}^3$ gs1")
# ax.plot(t_arr00001,s3disterrdatatsp23gs2dt00001, color='b',label = r"$\mathbf{S}^3$ gs2")
# ax.plot(t_arr00001,s3disterrdatatsp23gs3dt00001, color='k',label = r"$\mathbf{S}^3$ gs3")


# # # for h data s12
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-1.5}'
# # ylabels[1] = r'\textbf{-1.0}'
# # ylabels[2] = r'\textbf{-0.5}'
# # ylabels[3] = r'\textbf{0.0}'
# # ylabels[4] = r'\textbf{0.5}'
# # ylabels[5] = r'\textbf{1.0}'

# # # for h data s13
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-1.0}'
# # ylabels[1] = r'\textbf{-.5}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{0.5}'
# # ylabels[4] = r'\textbf{1.0}'
# # ylabels[5] = r'\textbf{1.5}'

# # # for h data s23
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-4.0}'
# # ylabels[1] = r'\textbf{-2.0}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{2.0}'
# # ylabels[4] = r'\textbf{4.0}'
# # ylabels[5] = r'\textbf{6.0}'

# # for s data s12
# ylabels = [item.get_text() for item in ax.get_yticklabels()]
# ylabels[0] = r'\textbf{-5.0}'
# ylabels[1] = r'\textbf{-5.0}'
# ylabels[2] = r'\textbf{0.0}'
# ylabels[3] = r'\textbf{5.0}'
# ylabels[4] = r'\textbf{5.0}'
# ylabels[5] = r'\textbf{7.5}'

# # ax.set_xticklabels(xlabels)
# ax.set_yticklabels(ylabels)
# ax.legend(fontsize="15",loc="upper left")
# plt.title(r'\textbf{Triangle Body in $\mathbf{S}^3$}', fontsize=20)
# # plt.title(r'\textbf{Triangle Body in $\mathbf{H}^3$}', fontsize=20)
# plt.ylabel(r'\textbf{Relative Error of \\ Vertex Separation ($\times 10^{-11}$)}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)


# plt.tight_layout()

# plt.show()


######### old stuff

# fig,ax=plt.subplots(1,1)

# ax.plot(h3distsp2gs1dt00001,h3disterrdatatsp23gs1dt00001, color='r',label = "H3 Gauss s1 dt00001")
# ax.plot(h3distsp2gs2dt00001,h3disterrdatatsp23gs2dt00001, color='b',label = "H3 Gauss s2 dt00001")
# ax.plot(h3distsp2gs3dt00001,h3disterrdatatsp23gs3dt00001, color='k',label = "H3 Gauss s3 dt00001")

# # ax.plot(t_arr00001,s3disterrdatatsp12gs1dt00001, color='r',label = "S3 Gauss s1 dt00001")
# # ax.plot(t_arr00001,s3disterrdatatsp12gs2dt00001, color='b',label = "S3 Gauss s2 dt00001")
# # ax.plot(t_arr00001,s3disterrdatatsp12gs3dt00001, color='k',label = "S3 Gauss s3 dt00001")

# # ax.scatter(timevals,h3errvals12s1,marker = 's', color='r',label = "H3 Gauss 12 s1")
# # ax.scatter(timevals,h3errvals13s2,marker = 's', color='b',label = "H3 Gauss 13 s2")
# # ax.scatter(timevals,h3errvals23s3,marker = 's', color='k',label = "H3 Gauss 23 s3")
# # ax.scatter(timevals,s3errvals12s1,marker = 'o',color='r',label = "S3 Gauss 12 s1")
# # ax.scatter(timevals,s3errvals13s2,marker = 'o',color='b',label = "S3 Gauss 13 s2")
# # ax.scatter(timevals,s3errvals23s3,marker = 'o',color='k',label = "S3 Gauss 23 s3")
# ax.legend()
# ax.set_title('Percent Error of Vertex Separation vs. Time')
# # ax.set_title('Vertex Separation vs. Midpoint Distance')
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# # ax.set_xlim(0,5)
# # ax.set_ylim(0-1e-11,0+1e-11)
# ax.set_xlabel('Time')
# # ax.set_xlabel('Midpoint Distance')
# ax.set_ylabel('Percent Error of Vertex Separation')
# plt.show()




######
# Spring Potential Triangle Separation Error Scaling plot for H3 and S3
######

# ### Gauss s1

# h3disterrdatatsp12gs1dt1raw     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3disterrdatatsp12gs1dt01raw    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3disterrdatatsp12gs1dt001raw   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3disterrdatatsp12gs1dt0001raw  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3disterrdatatsp12gs1dt00001raw = np.zeros(np.shape(datath3dt00001gs1)[0])

# h3disterrdatatsp13gs1dt1raw     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3disterrdatatsp13gs1dt01raw    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3disterrdatatsp13gs1dt001raw   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3disterrdatatsp13gs1dt0001raw  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3disterrdatatsp13gs1dt00001raw = np.zeros(np.shape(datath3dt00001gs1)[0])

# h3disterrdatatsp23gs1dt1raw     = np.zeros(np.shape(datath3dt1gs1)[0])
# h3disterrdatatsp23gs1dt01raw    = np.zeros(np.shape(datath3dt01gs1)[0])
# h3disterrdatatsp23gs1dt001raw   = np.zeros(np.shape(datath3dt001gs1)[0])
# h3disterrdatatsp23gs1dt0001raw  = np.zeros(np.shape(datath3dt0001gs1)[0])
# h3disterrdatatsp23gs1dt00001raw = np.zeros(np.shape(datath3dt00001gs1)[0])

# s3disterrdatatsp12gs1dt1raw     = np.zeros(np.shape(datats3dt1gs1)[0])
# s3disterrdatatsp12gs1dt01raw    = np.zeros(np.shape(datats3dt01gs1)[0])
# s3disterrdatatsp12gs1dt001raw   = np.zeros(np.shape(datats3dt001gs1)[0])
# s3disterrdatatsp12gs1dt0001raw  = np.zeros(np.shape(datats3dt0001gs1)[0])
# s3disterrdatatsp12gs1dt00001raw = np.zeros(np.shape(datats3dt00001gs1)[0])

# s3disterrdatatsp13gs1dt1raw     = np.zeros(np.shape(datats3dt1gs1)[0])
# s3disterrdatatsp13gs1dt01raw    = np.zeros(np.shape(datats3dt01gs1)[0])
# s3disterrdatatsp13gs1dt001raw   = np.zeros(np.shape(datats3dt001gs1)[0])
# s3disterrdatatsp13gs1dt0001raw  = np.zeros(np.shape(datats3dt0001gs1)[0])
# s3disterrdatatsp13gs1dt00001raw = np.zeros(np.shape(datats3dt00001gs1)[0])

# s3disterrdatatsp23gs1dt1raw     = np.zeros(np.shape(datats3dt1gs1)[0])
# s3disterrdatatsp23gs1dt01raw    = np.zeros(np.shape(datats3dt01gs1)[0])
# s3disterrdatatsp23gs1dt001raw   = np.zeros(np.shape(datats3dt001gs1)[0])
# s3disterrdatatsp23gs1dt0001raw  = np.zeros(np.shape(datats3dt0001gs1)[0])
# s3disterrdatatsp23gs1dt00001raw = np.zeros(np.shape(datats3dt00001gs1)[0])

# ### Gauss s2

# h3disterrdatatsp12gs2dt1raw     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3disterrdatatsp12gs2dt01raw    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3disterrdatatsp12gs2dt001raw   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3disterrdatatsp12gs2dt0001raw  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3disterrdatatsp12gs2dt00001raw = np.zeros(np.shape(datath3dt00001gs2)[0])

# h3disterrdatatsp13gs2dt1raw     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3disterrdatatsp13gs2dt01raw    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3disterrdatatsp13gs2dt001raw   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3disterrdatatsp13gs2dt0001raw  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3disterrdatatsp13gs2dt00001raw = np.zeros(np.shape(datath3dt00001gs2)[0])

# h3disterrdatatsp23gs2dt1raw     = np.zeros(np.shape(datath3dt1gs2)[0])
# h3disterrdatatsp23gs2dt01raw    = np.zeros(np.shape(datath3dt01gs2)[0])
# h3disterrdatatsp23gs2dt001raw   = np.zeros(np.shape(datath3dt001gs2)[0])
# h3disterrdatatsp23gs2dt0001raw  = np.zeros(np.shape(datath3dt0001gs2)[0])
# h3disterrdatatsp23gs2dt00001raw = np.zeros(np.shape(datath3dt00001gs2)[0])

# s3disterrdatatsp12gs2dt1raw     = np.zeros(np.shape(datats3dt1gs2)[0])
# s3disterrdatatsp12gs2dt01raw    = np.zeros(np.shape(datats3dt01gs2)[0])
# s3disterrdatatsp12gs2dt001raw   = np.zeros(np.shape(datats3dt001gs2)[0])
# s3disterrdatatsp12gs2dt0001raw  = np.zeros(np.shape(datats3dt0001gs2)[0])
# s3disterrdatatsp12gs2dt00001raw = np.zeros(np.shape(datats3dt00001gs2)[0])

# s3disterrdatatsp13gs2dt1raw     = np.zeros(np.shape(datats3dt1gs2)[0])
# s3disterrdatatsp13gs2dt01raw    = np.zeros(np.shape(datats3dt01gs2)[0])
# s3disterrdatatsp13gs2dt001raw   = np.zeros(np.shape(datats3dt001gs2)[0])
# s3disterrdatatsp13gs2dt0001raw  = np.zeros(np.shape(datats3dt0001gs2)[0])
# s3disterrdatatsp13gs2dt00001raw = np.zeros(np.shape(datats3dt00001gs2)[0])

# s3disterrdatatsp23gs2dt1raw     = np.zeros(np.shape(datats3dt1gs2)[0])
# s3disterrdatatsp23gs2dt01raw    = np.zeros(np.shape(datats3dt01gs2)[0])
# s3disterrdatatsp23gs2dt001raw   = np.zeros(np.shape(datats3dt001gs2)[0])
# s3disterrdatatsp23gs2dt0001raw  = np.zeros(np.shape(datats3dt0001gs2)[0])
# s3disterrdatatsp23gs2dt00001raw = np.zeros(np.shape(datats3dt00001gs2)[0])

# ### Gauss s3

# h3disterrdatatsp12gs3dt1raw     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3disterrdatatsp12gs3dt01raw    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3disterrdatatsp12gs3dt001raw   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3disterrdatatsp12gs3dt0001raw  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3disterrdatatsp12gs3dt00001raw = np.zeros(np.shape(datath3dt00001gs3)[0])

# h3disterrdatatsp13gs3dt1raw     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3disterrdatatsp13gs3dt01raw    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3disterrdatatsp13gs3dt001raw   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3disterrdatatsp13gs3dt0001raw  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3disterrdatatsp13gs3dt00001raw = np.zeros(np.shape(datath3dt00001gs3)[0])

# h3disterrdatatsp23gs3dt1raw     = np.zeros(np.shape(datath3dt1gs3)[0])
# h3disterrdatatsp23gs3dt01raw    = np.zeros(np.shape(datath3dt01gs3)[0])
# h3disterrdatatsp23gs3dt001raw   = np.zeros(np.shape(datath3dt001gs3)[0])
# h3disterrdatatsp23gs3dt0001raw  = np.zeros(np.shape(datath3dt0001gs3)[0])
# h3disterrdatatsp23gs3dt00001raw = np.zeros(np.shape(datath3dt00001gs3)[0])

# s3disterrdatatsp12gs3dt1raw     = np.zeros(np.shape(datats3dt1gs3)[0])
# s3disterrdatatsp12gs3dt01raw    = np.zeros(np.shape(datats3dt01gs3)[0])
# s3disterrdatatsp12gs3dt001raw   = np.zeros(np.shape(datats3dt001gs3)[0])
# s3disterrdatatsp12gs3dt0001raw  = np.zeros(np.shape(datats3dt0001gs3)[0])
# s3disterrdatatsp12gs3dt00001raw = np.zeros(np.shape(datats3dt00001gs3)[0])

# s3disterrdatatsp13gs3dt1raw     = np.zeros(np.shape(datats3dt1gs3)[0])
# s3disterrdatatsp13gs3dt01raw    = np.zeros(np.shape(datats3dt01gs3)[0])
# s3disterrdatatsp13gs3dt001raw   = np.zeros(np.shape(datats3dt001gs3)[0])
# s3disterrdatatsp13gs3dt0001raw  = np.zeros(np.shape(datats3dt0001gs3)[0])
# s3disterrdatatsp13gs3dt00001raw = np.zeros(np.shape(datats3dt00001gs3)[0])

# s3disterrdatatsp23gs3dt1raw     = np.zeros(np.shape(datats3dt1gs3)[0])
# s3disterrdatatsp23gs3dt01raw    = np.zeros(np.shape(datats3dt01gs3)[0])
# s3disterrdatatsp23gs3dt001raw   = np.zeros(np.shape(datats3dt001gs3)[0])
# s3disterrdatatsp23gs3dt0001raw  = np.zeros(np.shape(datats3dt0001gs3)[0])
# s3disterrdatatsp23gs3dt00001raw = np.zeros(np.shape(datats3dt00001gs3)[0])


# counter = 0
# print("Collating dt = .1 Data")
# for a in range(np.shape(datath3dt1gs1)[0]):
#     h3disterrdatatsp12gs1dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs1[a][0:3]),rot2hyp(datath3dt1gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs2[a][0:3]),rot2hyp(datath3dt1gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs3[a][0:3]),rot2hyp(datath3dt1gs3[a][3:6]))

#     h3disterrdatatsp13gs1dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs1[a][0:3]),rot2hyp(datath3dt1gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs2[a][0:3]),rot2hyp(datath3dt1gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs3[a][0:3]),rot2hyp(datath3dt1gs3[a][6:9]))

#     h3disterrdatatsp23gs1dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs1[a][3:6]),rot2hyp(datath3dt1gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs2[a][3:6]),rot2hyp(datath3dt1gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt1raw[counter] = h3dist(rot2hyp(datath3dt1gs3[a][3:6]),rot2hyp(datath3dt1gs3[a][6:9]))

#     s3disterrdatatsp12gs1dt1raw[counter] = r4dist(rot2r4(datats3dt1gs1[a][0:3]),rot2r4(datats3dt1gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt1raw[counter] = r4dist(rot2r4(datats3dt1gs2[a][0:3]),rot2r4(datats3dt1gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt1raw[counter] = r4dist(rot2r4(datats3dt1gs3[a][0:3]),rot2r4(datats3dt1gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt1raw[counter] = r4dist(rot2r4(datats3dt1gs1[a][0:3]),rot2r4(datats3dt1gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt1raw[counter] = r4dist(rot2r4(datats3dt1gs2[a][0:3]),rot2r4(datats3dt1gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt1raw[counter] = r4dist(rot2r4(datats3dt1gs3[a][0:3]),rot2r4(datats3dt1gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt1raw[counter] = r4dist(rot2r4(datats3dt1gs1[a][3:6]),rot2r4(datats3dt1gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt1raw[counter] = r4dist(rot2r4(datats3dt1gs2[a][3:6]),rot2r4(datats3dt1gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt1raw[counter] = r4dist(rot2r4(datats3dt1gs3[a][3:6]),rot2r4(datats3dt1gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp12gs1dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp12gs2dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp12gs2dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp12gs3dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp12gs3dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))

# h3disterrdatatsp13gs1dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp13gs1dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp13gs2dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp13gs2dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))
# h3disterrdatatsp13gs3dt1 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])) - h3disterrdatatsp13gs3dt1raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10000,0])*np.cosh(datath3gtdt00001gs3[::10000,1])))

# h3disterrdatatsp23gs1dt1 = (2.*datath3gtdt00001gs3[::10000,0] - h3disterrdatatsp23gs1dt1raw)/(2.*datath3gtdt00001gs3[::10000,0])
# h3disterrdatatsp23gs2dt1 = (2.*datath3gtdt00001gs3[::10000,0] - h3disterrdatatsp23gs2dt1raw)/(2.*datath3gtdt00001gs3[::10000,0])
# h3disterrdatatsp23gs3dt1 = (2.*datath3gtdt00001gs3[::10000,0] - h3disterrdatatsp23gs3dt1raw)/(2.*datath3gtdt00001gs3[::10000,0])

# s3disterrdatatsp12gs1dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp12gs1dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp12gs2dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp12gs2dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp12gs3dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp12gs3dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))

# s3disterrdatatsp13gs1dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp13gs1dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp13gs2dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp13gs2dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))
# s3disterrdatatsp13gs3dt1 = (np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])) - s3disterrdatatsp13gs3dt1raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10000,0])*np.cos(datats3gtdt00001gs3[::10000,1])))

# s3disterrdatatsp23gs1dt1 = (2.*datats3gtdt00001gs3[::10000,0] - s3disterrdatatsp23gs1dt1raw)/(2.*datats3gtdt00001gs3[::10000,0])
# s3disterrdatatsp23gs2dt1 = (2.*datats3gtdt00001gs3[::10000,0] - s3disterrdatatsp23gs2dt1raw)/(2.*datats3gtdt00001gs3[::10000,0])
# s3disterrdatatsp23gs3dt1 = (2.*datats3gtdt00001gs3[::10000,0] - s3disterrdatatsp23gs3dt1raw)/(2.*datats3gtdt00001gs3[::10000,0])

# counter = 0
# print("Collating dt = .01 Data")
# for a in range(np.shape(datath3dt01gs1)[0]):
#     h3disterrdatatsp12gs1dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs1[a][0:3]),rot2hyp(datath3dt01gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs2[a][0:3]),rot2hyp(datath3dt01gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs3[a][0:3]),rot2hyp(datath3dt01gs3[a][3:6]))

#     h3disterrdatatsp13gs1dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs1[a][0:3]),rot2hyp(datath3dt01gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs2[a][0:3]),rot2hyp(datath3dt01gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs3[a][0:3]),rot2hyp(datath3dt01gs3[a][6:9]))

#     h3disterrdatatsp23gs1dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs1[a][3:6]),rot2hyp(datath3dt01gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs2[a][3:6]),rot2hyp(datath3dt01gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt01raw[counter] = h3dist(rot2hyp(datath3dt01gs3[a][3:6]),rot2hyp(datath3dt01gs3[a][6:9]))

#     s3disterrdatatsp12gs1dt01raw[counter] = r4dist(rot2r4(datats3dt01gs1[a][0:3]),rot2r4(datats3dt01gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt01raw[counter] = r4dist(rot2r4(datats3dt01gs2[a][0:3]),rot2r4(datats3dt01gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt01raw[counter] = r4dist(rot2r4(datats3dt01gs3[a][0:3]),rot2r4(datats3dt01gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt01raw[counter] = r4dist(rot2r4(datats3dt01gs1[a][0:3]),rot2r4(datats3dt01gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt01raw[counter] = r4dist(rot2r4(datats3dt01gs2[a][0:3]),rot2r4(datats3dt01gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt01raw[counter] = r4dist(rot2r4(datats3dt01gs3[a][0:3]),rot2r4(datats3dt01gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt01raw[counter] = r4dist(rot2r4(datats3dt01gs1[a][3:6]),rot2r4(datats3dt01gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt01raw[counter] = r4dist(rot2r4(datats3dt01gs2[a][3:6]),rot2r4(datats3dt01gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt01raw[counter] = r4dist(rot2r4(datats3dt01gs3[a][3:6]),rot2r4(datats3dt01gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp12gs1dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp12gs2dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp12gs2dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp12gs3dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp12gs3dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))

# h3disterrdatatsp13gs1dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp13gs1dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp13gs2dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp13gs2dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))
# h3disterrdatatsp13gs3dt01 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])) - h3disterrdatatsp13gs3dt01raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1000,0])*np.cosh(datath3gtdt00001gs3[::1000,1])))

# h3disterrdatatsp23gs1dt01 = (2.*datath3gtdt00001gs3[::1000,0] - h3disterrdatatsp23gs1dt01raw)/(2.*datath3gtdt00001gs3[::1000,0])
# h3disterrdatatsp23gs2dt01 = (2.*datath3gtdt00001gs3[::1000,0] - h3disterrdatatsp23gs2dt01raw)/(2.*datath3gtdt00001gs3[::1000,0])
# h3disterrdatatsp23gs3dt01 = (2.*datath3gtdt00001gs3[::1000,0] - h3disterrdatatsp23gs3dt01raw)/(2.*datath3gtdt00001gs3[::1000,0])

# s3disterrdatatsp12gs1dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp12gs1dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp12gs2dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp12gs2dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp12gs3dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp12gs3dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))

# s3disterrdatatsp13gs1dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp13gs1dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp13gs2dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp13gs2dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))
# s3disterrdatatsp13gs3dt01 = (np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])) - s3disterrdatatsp13gs3dt01raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1000,0])*np.cos(datats3gtdt00001gs3[::1000,1])))

# s3disterrdatatsp23gs1dt01 = (2.*datats3gtdt00001gs3[::1000,0] - s3disterrdatatsp23gs1dt01raw)/(2.*datats3gtdt00001gs3[::1000,0])
# s3disterrdatatsp23gs2dt01 = (2.*datats3gtdt00001gs3[::1000,0] - s3disterrdatatsp23gs2dt01raw)/(2.*datats3gtdt00001gs3[::1000,0])
# s3disterrdatatsp23gs3dt01 = (2.*datats3gtdt00001gs3[::1000,0] - s3disterrdatatsp23gs3dt01raw)/(2.*datats3gtdt00001gs3[::1000,0])

# counter = 0
# print("Collating dt = .001 Data")
# for a in range(np.shape(datath3dt001gs1)[0]):
#     h3disterrdatatsp12gs1dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs1[a][0:3]),rot2hyp(datath3dt001gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs2[a][0:3]),rot2hyp(datath3dt001gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs3[a][0:3]),rot2hyp(datath3dt001gs3[a][3:6]))

#     h3disterrdatatsp13gs1dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs1[a][0:3]),rot2hyp(datath3dt001gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs2[a][0:3]),rot2hyp(datath3dt001gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs3[a][0:3]),rot2hyp(datath3dt001gs3[a][6:9]))

#     h3disterrdatatsp23gs1dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs1[a][3:6]),rot2hyp(datath3dt001gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs2[a][3:6]),rot2hyp(datath3dt001gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt001raw[counter] = h3dist(rot2hyp(datath3dt001gs3[a][3:6]),rot2hyp(datath3dt001gs3[a][6:9]))

#     s3disterrdatatsp12gs1dt001raw[counter] = r4dist(rot2r4(datats3dt001gs1[a][0:3]),rot2r4(datats3dt001gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt001raw[counter] = r4dist(rot2r4(datats3dt001gs2[a][0:3]),rot2r4(datats3dt001gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt001raw[counter] = r4dist(rot2r4(datats3dt001gs3[a][0:3]),rot2r4(datats3dt001gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt001raw[counter] = r4dist(rot2r4(datats3dt001gs1[a][0:3]),rot2r4(datats3dt001gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt001raw[counter] = r4dist(rot2r4(datats3dt001gs2[a][0:3]),rot2r4(datats3dt001gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt001raw[counter] = r4dist(rot2r4(datats3dt001gs3[a][0:3]),rot2r4(datats3dt001gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt001raw[counter] = r4dist(rot2r4(datats3dt001gs1[a][3:6]),rot2r4(datats3dt001gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt001raw[counter] = r4dist(rot2r4(datats3dt001gs2[a][3:6]),rot2r4(datats3dt001gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt001raw[counter] = r4dist(rot2r4(datats3dt001gs3[a][3:6]),rot2r4(datats3dt001gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp12gs1dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp12gs2dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp12gs2dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp12gs3dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp12gs3dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))

# h3disterrdatatsp13gs1dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp13gs1dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp13gs2dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp13gs2dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))
# h3disterrdatatsp13gs3dt001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])) - h3disterrdatatsp13gs3dt001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::100,0])*np.cosh(datath3gtdt00001gs3[::100,1])))

# h3disterrdatatsp23gs1dt001 = (2.*datath3gtdt00001gs3[::100,0] - h3disterrdatatsp23gs1dt001raw)/(2.*datath3gtdt00001gs3[::100,0])
# h3disterrdatatsp23gs2dt001 = (2.*datath3gtdt00001gs3[::100,0] - h3disterrdatatsp23gs2dt001raw)/(2.*datath3gtdt00001gs3[::100,0])
# h3disterrdatatsp23gs3dt001 = (2.*datath3gtdt00001gs3[::100,0] - h3disterrdatatsp23gs3dt001raw)/(2.*datath3gtdt00001gs3[::100,0])

# s3disterrdatatsp12gs1dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp12gs1dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp12gs2dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp12gs2dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp12gs3dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp12gs3dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))

# s3disterrdatatsp13gs1dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp13gs1dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp13gs2dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp13gs2dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))
# s3disterrdatatsp13gs3dt001 = (np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])) - s3disterrdatatsp13gs3dt001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::100,0])*np.cos(datats3gtdt00001gs3[::100,1])))

# s3disterrdatatsp23gs1dt001 = (2.*datats3gtdt00001gs3[::100,0] - s3disterrdatatsp23gs1dt001raw)/(2.*datats3gtdt00001gs3[::100,0])
# s3disterrdatatsp23gs2dt001 = (2.*datats3gtdt00001gs3[::100,0] - s3disterrdatatsp23gs2dt001raw)/(2.*datats3gtdt00001gs3[::100,0])
# s3disterrdatatsp23gs3dt001 = (2.*datats3gtdt00001gs3[::100,0] - s3disterrdatatsp23gs3dt001raw)/(2.*datats3gtdt00001gs3[::100,0])

# counter = 0
# print("Collating dt = .0001 Data")
# for a in range(np.shape(datath3dt0001gs1)[0]):
#     h3disterrdatatsp12gs1dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs1[a][0:3]),rot2hyp(datath3dt0001gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs2[a][0:3]),rot2hyp(datath3dt0001gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][0:3]),rot2hyp(datath3dt0001gs3[a][3:6]))

#     h3disterrdatatsp13gs1dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs1[a][0:3]),rot2hyp(datath3dt0001gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs2[a][0:3]),rot2hyp(datath3dt0001gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][0:3]),rot2hyp(datath3dt0001gs3[a][6:9]))

#     h3disterrdatatsp23gs1dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs1[a][3:6]),rot2hyp(datath3dt0001gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs2[a][3:6]),rot2hyp(datath3dt0001gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt0001raw[counter] = h3dist(rot2hyp(datath3dt0001gs3[a][3:6]),rot2hyp(datath3dt0001gs3[a][6:9]))

#     s3disterrdatatsp12gs1dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs1[a][0:3]),rot2r4(datats3dt0001gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs2[a][0:3]),rot2r4(datats3dt0001gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs3[a][0:3]),rot2r4(datats3dt0001gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs1[a][0:3]),rot2r4(datats3dt0001gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs2[a][0:3]),rot2r4(datats3dt0001gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs3[a][0:3]),rot2r4(datats3dt0001gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs1[a][3:6]),rot2r4(datats3dt0001gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs2[a][3:6]),rot2r4(datats3dt0001gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt0001raw[counter] = r4dist(rot2r4(datats3dt0001gs3[a][3:6]),rot2r4(datats3dt0001gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp12gs1dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp12gs2dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp12gs2dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp12gs3dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp12gs3dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))

# h3disterrdatatsp13gs1dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp13gs1dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp13gs2dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp13gs2dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))
# h3disterrdatatsp13gs3dt0001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])) - h3disterrdatatsp13gs3dt0001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::10,0])*np.cosh(datath3gtdt00001gs3[::10,1])))

# h3disterrdatatsp23gs1dt0001 = (2.*datath3gtdt00001gs3[::10,0] - h3disterrdatatsp23gs1dt0001raw)/(2.*datath3gtdt00001gs3[::10,0])
# h3disterrdatatsp23gs2dt0001 = (2.*datath3gtdt00001gs3[::10,0] - h3disterrdatatsp23gs2dt0001raw)/(2.*datath3gtdt00001gs3[::10,0])
# h3disterrdatatsp23gs3dt0001 = (2.*datath3gtdt00001gs3[::10,0] - h3disterrdatatsp23gs3dt0001raw)/(2.*datath3gtdt00001gs3[::10,0])

# s3disterrdatatsp12gs1dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp12gs1dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp12gs2dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp12gs2dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp12gs3dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp12gs3dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))

# s3disterrdatatsp13gs1dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp13gs1dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp13gs2dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp13gs2dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))
# s3disterrdatatsp13gs3dt0001 = (np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])) - s3disterrdatatsp13gs3dt0001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::10,0])*np.cos(datats3gtdt00001gs3[::10,1])))

# s3disterrdatatsp23gs1dt0001 = (2.*datats3gtdt00001gs3[::10,0] - s3disterrdatatsp23gs1dt0001raw)/(2.*datats3gtdt00001gs3[::10,0])
# s3disterrdatatsp23gs2dt0001 = (2.*datats3gtdt00001gs3[::10,0] - s3disterrdatatsp23gs2dt0001raw)/(2.*datats3gtdt00001gs3[::10,0])
# s3disterrdatatsp23gs3dt0001 = (2.*datats3gtdt00001gs3[::10,0] - s3disterrdatatsp23gs3dt0001raw)/(2.*datats3gtdt00001gs3[::10,0])

# counter = 0
# print("Collating dt = .00001 Data")
# for a in range(np.shape(datath3dt00001gs1)[0]):
#     h3disterrdatatsp12gs1dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs1[a][0:3]),rot2hyp(datath3dt00001gs1[a][3:6]))
#     h3disterrdatatsp12gs2dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs2[a][0:3]),rot2hyp(datath3dt00001gs2[a][3:6]))
#     h3disterrdatatsp12gs3dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][0:3]),rot2hyp(datath3dt00001gs3[a][3:6]))

#     h3disterrdatatsp13gs1dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs1[a][0:3]),rot2hyp(datath3dt00001gs1[a][6:9]))
#     h3disterrdatatsp13gs2dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs2[a][0:3]),rot2hyp(datath3dt00001gs2[a][6:9]))
#     h3disterrdatatsp13gs3dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][0:3]),rot2hyp(datath3dt00001gs3[a][6:9]))

#     h3disterrdatatsp23gs1dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs1[a][3:6]),rot2hyp(datath3dt00001gs1[a][6:9]))
#     h3disterrdatatsp23gs2dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs2[a][3:6]),rot2hyp(datath3dt00001gs2[a][6:9]))
#     h3disterrdatatsp23gs3dt00001raw[counter] = h3dist(rot2hyp(datath3dt00001gs3[a][3:6]),rot2hyp(datath3dt00001gs3[a][6:9]))

#     s3disterrdatatsp12gs1dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs1[a][0:3]),rot2r4(datats3dt00001gs1[a][3:6]))
#     s3disterrdatatsp12gs2dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs2[a][0:3]),rot2r4(datats3dt00001gs2[a][3:6]))
#     s3disterrdatatsp12gs3dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs3[a][0:3]),rot2r4(datats3dt00001gs3[a][3:6]))

#     s3disterrdatatsp13gs1dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs1[a][0:3]),rot2r4(datats3dt00001gs1[a][6:9]))
#     s3disterrdatatsp13gs2dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs2[a][0:3]),rot2r4(datats3dt00001gs2[a][6:9]))
#     s3disterrdatatsp13gs3dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs3[a][0:3]),rot2r4(datats3dt00001gs3[a][6:9]))

#     s3disterrdatatsp23gs1dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs1[a][3:6]),rot2r4(datats3dt00001gs1[a][6:9]))
#     s3disterrdatatsp23gs2dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs2[a][3:6]),rot2r4(datats3dt00001gs2[a][6:9]))
#     s3disterrdatatsp23gs3dt00001raw[counter] = r4dist(rot2r4(datats3dt00001gs3[a][3:6]),rot2r4(datats3dt00001gs3[a][6:9]))

#     counter += 1

# h3disterrdatatsp12gs1dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp12gs1dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp12gs2dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp12gs2dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp12gs3dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp12gs3dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))

# h3disterrdatatsp13gs1dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp13gs1dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp13gs2dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp13gs2dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))
# h3disterrdatatsp13gs3dt00001 = (np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])) - h3disterrdatatsp13gs3dt00001raw)/(np.arccosh(np.cosh(datath3gtdt00001gs3[::1,0])*np.cosh(datath3gtdt00001gs3[::1,1])))

# h3disterrdatatsp23gs1dt00001 = (2.*datath3gtdt00001gs3[::1,0] - h3disterrdatatsp23gs1dt00001raw)/(2.*datath3gtdt00001gs3[::1,0])
# h3disterrdatatsp23gs2dt00001 = (2.*datath3gtdt00001gs3[::1,0] - h3disterrdatatsp23gs2dt00001raw)/(2.*datath3gtdt00001gs3[::1,0])
# h3disterrdatatsp23gs3dt00001 = (2.*datath3gtdt00001gs3[::1,0] - h3disterrdatatsp23gs3dt00001raw)/(2.*datath3gtdt00001gs3[::1,0])

# s3disterrdatatsp12gs1dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp12gs1dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp12gs2dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp12gs2dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp12gs3dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp12gs3dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))

# s3disterrdatatsp13gs1dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp13gs1dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp13gs2dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp13gs2dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))
# s3disterrdatatsp13gs3dt00001 = (np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])) - s3disterrdatatsp13gs3dt00001raw)/(np.arccos(np.cos(datats3gtdt00001gs3[::1,0])*np.cos(datats3gtdt00001gs3[::1,1])))

# s3disterrdatatsp23gs1dt00001 = (2.*datats3gtdt00001gs3[::1,0] - s3disterrdatatsp23gs1dt00001raw)/(2.*datats3gtdt00001gs3[::1,0])
# s3disterrdatatsp23gs2dt00001 = (2.*datats3gtdt00001gs3[::1,0] - s3disterrdatatsp23gs2dt00001raw)/(2.*datats3gtdt00001gs3[::1,0])
# s3disterrdatatsp23gs3dt00001 = (2.*datats3gtdt00001gs3[::1,0] - s3disterrdatatsp23gs3dt00001raw)/(2.*datats3gtdt00001gs3[::1,0])

# timevals = [.1,.01,.001,.0001,.00001]

# timevals2 = [.1**2,.01**2,.001**2,.0001**2,.00001**2]
# timevals4 = [.1**4,.01**4,.001**4,.0001**4,.00001**4]
# timevals6 = [.1**6,.01**6,.001**6,.0001**6,.00001**6]

# h3errvals12s1 = [max(h3disterrdatatsp12gs1dt1),max(h3disterrdatatsp12gs1dt01),max(h3disterrdatatsp12gs1dt001),max(h3disterrdatatsp12gs1dt0001),max(h3disterrdatatsp12gs1dt00001)]
# h3errvals12s2 = [max(h3disterrdatatsp12gs2dt1),max(h3disterrdatatsp12gs2dt01),max(h3disterrdatatsp12gs2dt001),max(h3disterrdatatsp12gs2dt0001),max(h3disterrdatatsp12gs2dt00001)]
# h3errvals12s3 = [max(h3disterrdatatsp12gs3dt1),max(h3disterrdatatsp12gs3dt01),max(h3disterrdatatsp12gs3dt001),max(h3disterrdatatsp12gs3dt0001),max(h3disterrdatatsp12gs3dt00001)]

# h3errvals13s1 = [max(h3disterrdatatsp13gs1dt1),max(h3disterrdatatsp13gs1dt01),max(h3disterrdatatsp13gs1dt001),max(h3disterrdatatsp13gs1dt0001),max(h3disterrdatatsp13gs1dt00001)]
# h3errvals13s2 = [max(h3disterrdatatsp13gs2dt1),max(h3disterrdatatsp13gs2dt01),max(h3disterrdatatsp13gs2dt001),max(h3disterrdatatsp13gs2dt0001),max(h3disterrdatatsp13gs2dt00001)]
# h3errvals13s3 = [max(h3disterrdatatsp13gs3dt1),max(h3disterrdatatsp13gs3dt01),max(h3disterrdatatsp13gs3dt001),max(h3disterrdatatsp13gs3dt0001),max(h3disterrdatatsp13gs3dt00001)]

# h3errvals23s1 = [max(h3disterrdatatsp23gs1dt1),max(h3disterrdatatsp23gs1dt01),max(h3disterrdatatsp23gs1dt001),max(h3disterrdatatsp23gs1dt0001),max(h3disterrdatatsp23gs1dt00001)]
# h3errvals23s2 = [max(h3disterrdatatsp23gs2dt1),max(h3disterrdatatsp23gs2dt01),max(h3disterrdatatsp23gs2dt001),max(h3disterrdatatsp23gs2dt0001),max(h3disterrdatatsp23gs2dt00001)]
# h3errvals23s3 = [max(h3disterrdatatsp23gs3dt1),max(h3disterrdatatsp23gs3dt01),max(h3disterrdatatsp23gs3dt001),max(h3disterrdatatsp23gs3dt0001),max(h3disterrdatatsp23gs3dt00001)]

# s3errvals12s1 = [max(s3disterrdatatsp12gs1dt1),max(s3disterrdatatsp12gs1dt01),max(s3disterrdatatsp12gs1dt001),max(s3disterrdatatsp12gs1dt0001),max(s3disterrdatatsp12gs1dt00001)]
# s3errvals12s2 = [max(s3disterrdatatsp12gs2dt1),max(s3disterrdatatsp12gs2dt01),max(s3disterrdatatsp12gs2dt001),max(s3disterrdatatsp12gs2dt0001),max(s3disterrdatatsp12gs2dt00001)]
# s3errvals12s3 = [max(s3disterrdatatsp12gs3dt1),max(s3disterrdatatsp12gs3dt01),max(s3disterrdatatsp12gs3dt001),max(s3disterrdatatsp12gs3dt0001),max(s3disterrdatatsp12gs3dt00001)]

# s3errvals13s1 = [max(s3disterrdatatsp13gs1dt1),max(s3disterrdatatsp13gs1dt01),max(s3disterrdatatsp13gs1dt001),max(s3disterrdatatsp13gs1dt0001),max(s3disterrdatatsp13gs1dt00001)]
# s3errvals13s2 = [max(s3disterrdatatsp13gs2dt1),max(s3disterrdatatsp13gs2dt01),max(s3disterrdatatsp13gs2dt001),max(s3disterrdatatsp13gs2dt0001),max(s3disterrdatatsp13gs2dt00001)]
# s3errvals13s3 = [max(s3disterrdatatsp13gs3dt1),max(s3disterrdatatsp13gs3dt01),max(s3disterrdatatsp13gs3dt001),max(s3disterrdatatsp13gs3dt0001),max(s3disterrdatatsp13gs3dt00001)]

# s3errvals23s1 = [max(s3disterrdatatsp23gs1dt1),max(s3disterrdatatsp23gs1dt01),max(s3disterrdatatsp23gs1dt001),max(s3disterrdatatsp23gs1dt0001),max(s3disterrdatatsp23gs1dt00001)]
# s3errvals23s2 = [max(s3disterrdatatsp23gs2dt1),max(s3disterrdatatsp23gs2dt01),max(s3disterrdatatsp23gs2dt001),max(s3disterrdatatsp23gs2dt0001),max(s3disterrdatatsp23gs2dt00001)]
# s3errvals23s3 = [max(s3disterrdatatsp23gs3dt1),max(s3disterrdatatsp23gs3dt01),max(s3disterrdatatsp23gs3dt001),max(s3disterrdatatsp23gs3dt0001),max(s3disterrdatatsp23gs3dt00001)]

# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# # ax.scatter(timevals,h3errvals12s1,marker = 's', color='r',label = r"$\mathbf{H}^3$ gs1 s12")
# # ax.scatter(timevals,h3errvals12s2,marker = 's', color='b',label = r"$\mathbf{H}^3$ gs2 s12")
# # ax.scatter(timevals,h3errvals12s3,marker = 's', color='k',label = r"$\mathbf{H}^3$ gs3 s12")

# # ax.scatter(timevals,h3errvals13s1,marker = 's', color='r',label = r"$\mathbf{H}^3$ gs1 s13")
# # ax.scatter(timevals,h3errvals13s2,marker = 's', color='b',label = r"$\mathbf{H}^3$ gs2 s13")
# # ax.scatter(timevals,h3errvals13s3,marker = 's', color='k',label = r"$\mathbf{H}^3$ gs3 s13")

# # ax.scatter(timevals,h3errvals23s1,marker = 's', color='r',label = r"$\mathbf{H}^3$ gs1 s23")
# # ax.scatter(timevals,h3errvals23s2,marker = 's', color='b',label = r"$\mathbf{H}^3$ gs2 s23")
# # ax.scatter(timevals,h3errvals23s3,marker = 's', color='k',label = r"$\mathbf{H}^3$ gs3 s23")

# # ax.scatter(timevals,s3errvals12s1,marker = 's', color='r',label = r"$\mathbf{S}^3$ gs1 s12")
# # ax.scatter(timevals,s3errvals12s2,marker = 's', color='b',label = r"$\mathbf{S}^3$ gs2 s12")
# # ax.scatter(timevals,s3errvals12s3,marker = 's', color='k',label = r"$\mathbf{S}^3$ gs3 s12")

# # ax.scatter(timevals,s3errvals13s1,marker = 's', color='r',label = r"$\mathbf{S}^3$ gs1 s13")
# # ax.scatter(timevals,s3errvals13s2,marker = 's', color='b',label = r"$\mathbf{S}^3$ gs2 s13")
# # ax.scatter(timevals,s3errvals13s3,marker = 's', color='k',label = r"$\mathbf{S}^3$ gs3 s13")

# ax.scatter(timevals,s3errvals23s1,marker = 's', color='r',label = r"$\mathbf{S}^3$ gs1 s23")
# ax.scatter(timevals,s3errvals23s2,marker = 's', color='b',label = r"$\mathbf{S}^3$ gs2 s23")
# ax.scatter(timevals,s3errvals23s3,marker = 's', color='k',label = r"$\mathbf{S}^3$ gs3 s23")

# ax.plot(timevals,timevals2,color='r',label = r"$dt^2$")
# ax.plot(timevals,timevals4,color='b',label = r"$dt^4$")
# ax.plot(timevals,timevals6,color='k',label = r"$dt^6$")



# # # for h data s12
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-1.5}'
# # ylabels[1] = r'\textbf{-1.0}'
# # ylabels[2] = r'\textbf{-0.5}'
# # ylabels[3] = r'\textbf{0.0}'
# # ylabels[4] = r'\textbf{0.5}'
# # ylabels[5] = r'\textbf{1.0}'


# # ax.set_xticklabels(xlabels)
# # ax.set_yticklabels(ylabels)
# ax.legend(fontsize="15",loc="upper left",ncol=2)
# plt.title(r'\textbf{ Method Error Scaling for Triangle Body in $\mathbf{S}^3$}', fontsize=20)
# # plt.title(r'\textbf{ Method Error Scaling for Triangle Body in $\mathbf{H}^3$}', fontsize=20)
# plt.ylabel(r'\textbf{Maximum Relative \\ \\ Error of Vertex Separation}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time Step Size dt}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)

# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim(1e-14,1e-0)

# plt.tight_layout()

# plt.show()

###### old stuff


# fig,ax=plt.subplots(1,1)

# # ax.scatter(timevals,h3errvals12s1,marker = 's', color='r',label = "H3 Gauss 12 s1")
# # ax.scatter(timevals,h3errvals13s2,marker = 's', color='b',label = "H3 Gauss 13 s2")
# # ax.scatter(timevals,h3errvals23s3,marker = 's', color='k',label = "H3 Gauss 23 s3")
# ax.scatter(timevals,s3errvals12s1,marker = 'o',color='r',label = "S3 Gauss 12 s1")
# ax.scatter(timevals,s3errvals13s2,marker = 'o',color='b',label = "S3 Gauss 13 s2")
# ax.scatter(timevals,s3errvals23s3,marker = 'o',color='k',label = "S3 Gauss 23 s3")
# ax.plot(timevals,timevals2,color='r',label = "h^2")
# ax.plot(timevals,timevals4,color='b',label = "h^4")
# ax.plot(timevals,timevals6,color='k',label = "h^6")
# ax.legend()
# ax.set_title('Vertex Separation Percent Error vs. Time')
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim(1e-14,1e-1)
# ax.set_xlabel('Time')
# ax.set_ylabel('Vertex Separation Percent Error')
# plt.show()


##################################################

######### Rigid Rod plots for H3 and S3 ##########

##################################################


#####
# Rigid Rod Killing flow difference plot for H3 and S3
#####

# ### Radau s2

# h3diffdatarrads2dt1raw   = np.zeros((np.shape(datarigrh3dt1rs2)[0],2))
# h3diffdatarrads2dt01raw  = np.zeros((np.shape(datarigrh3dt01rs2)[0],2))
# h3diffdatarrads2dt001raw = np.zeros((np.shape(datarigrh3dt001rs2)[0],2))

# s3diffdatarrads2dt1raw   = np.zeros((np.shape(datarigrs3dt1rs2)[0],2))
# s3diffdatarrads2dt01raw  = np.zeros((np.shape(datarigrs3dt01rs2)[0],2))
# s3diffdatarrads2dt001raw = np.zeros((np.shape(datarigrs3dt001rs2)[0],2))

# ### Radau s3

# h3diffdatarrads3dt1raw   = np.zeros((np.shape(datarigrh3dt1rs3)[0],2))
# h3diffdatarrads3dt01raw  = np.zeros((np.shape(datarigrh3dt01rs3)[0],2))
# h3diffdatarrads3dt001raw = np.zeros((np.shape(datarigrh3dt001rs3)[0],2))

# s3diffdatarrads3dt1raw   = np.zeros((np.shape(datarigrs3dt1rs3)[0],2))
# s3diffdatarrads3dt01raw  = np.zeros((np.shape(datarigrs3dt01rs3)[0],2))
# s3diffdatarrads3dt001raw = np.zeros((np.shape(datarigrs3dt001rs3)[0],2))


# counter = 0
# print("Collating dt = .1 Data")
# for a in range(np.shape(datarigrh3dt1rs2)[0]):
#     # Additional angle for trig function conditioning of the data
#     h3diffdatarrads2dt1raw[counter][0] = h3distb2((boostxh3(a*.1) @ rot2hyp(datarigrh3dt1rs2[0][0:3])),rot2hyp(datarigrh3dt1rs2[a][0:3])) - 1.
#     h3diffdatarrads2dt1raw[counter][1] = h3distb2((boostxh3(a*.1) @ rot2hyp(datarigrh3dt1rs2[0][3:6])),rot2hyp(datarigrh3dt1rs2[a][3:6])) - 1.
#     h3diffdatarrads3dt1raw[counter][0] = h3distb2((boostxh3(a*.1) @ rot2hyp(datarigrh3dt1rs3[0][0:3])),rot2hyp(datarigrh3dt1rs3[a][0:3])) - 1.
#     h3diffdatarrads3dt1raw[counter][1] = h3distb2((boostxh3(a*.1) @ rot2hyp(datarigrh3dt1rs3[0][3:6])),rot2hyp(datarigrh3dt1rs3[a][3:6])) - 1.

#     s3diffdatarrads2dt1raw[counter][0] = r4distb2((rotzh3(a*.1) @ rot2r4(datarigrs3dt1rs2[0][0:3])),rot2r4(datarigrs3dt1rs2[a][0:3])) - 1.
#     s3diffdatarrads2dt1raw[counter][1] = r4distb2((rotzh3(a*.1) @ rot2r4(datarigrs3dt1rs2[0][3:6])),rot2r4(datarigrs3dt1rs2[a][3:6])) - 1.
#     s3diffdatarrads3dt1raw[counter][0] = r4distb2((rotzh3(a*.1) @ rot2r4(datarigrs3dt1rs3[0][0:3])),rot2r4(datarigrs3dt1rs3[a][0:3])) - 1.
#     s3diffdatarrads3dt1raw[counter][1] = r4distb2((rotzh3(a*.1) @ rot2r4(datarigrs3dt1rs3[0][3:6])),rot2r4(datarigrs3dt1rs3[a][3:6])) - 1.

#     # s3lamdatarrads2dt1raw[counter] = (s3anal_value - datarigrs3dt1rs2[a][-1])/s3anal_value
#     # s3lamdatarrads3dt1raw[counter] = (s3anal_value - datarigrs3dt1rs3[a][-1])/s3anal_value

#     counter += 1


# counter = 0
# print("Collating dt = .01 Data")
# for a in range(np.shape(datarigrh3dt01rs2)[0]):
#     h3diffdatarrads2dt01raw[counter][0] = h3distb2((boostxh3(a*.01) @ rot2hyp(datarigrh3dt01rs2[0][0:3])),rot2hyp(datarigrh3dt01rs2[a][0:3])) - 1.
#     h3diffdatarrads2dt01raw[counter][1] = h3distb2((boostxh3(a*.01) @ rot2hyp(datarigrh3dt01rs2[0][3:6])),rot2hyp(datarigrh3dt01rs2[a][3:6])) - 1.
#     h3diffdatarrads3dt01raw[counter][0] = h3distb2((boostxh3(a*.01) @ rot2hyp(datarigrh3dt01rs3[0][0:3])),rot2hyp(datarigrh3dt01rs3[a][0:3])) - 1.
#     h3diffdatarrads3dt01raw[counter][1] = h3distb2((boostxh3(a*.01) @ rot2hyp(datarigrh3dt01rs3[0][3:6])),rot2hyp(datarigrh3dt01rs3[a][3:6])) - 1.

#     s3diffdatarrads2dt01raw[counter][0] = r4distb2((rotzh3(a*.01) @ rot2r4(datarigrs3dt01rs2[0][0:3])),rot2r4(datarigrs3dt01rs2[a][0:3])) - 1.
#     s3diffdatarrads2dt01raw[counter][1] = r4distb2((rotzh3(a*.01) @ rot2r4(datarigrs3dt01rs2[0][3:6])),rot2r4(datarigrs3dt01rs2[a][3:6])) - 1.
#     s3diffdatarrads3dt01raw[counter][0] = r4distb2((rotzh3(a*.01) @ rot2r4(datarigrs3dt01rs3[0][0:3])),rot2r4(datarigrs3dt01rs3[a][0:3])) - 1.
#     s3diffdatarrads3dt01raw[counter][1] = r4distb2((rotzh3(a*.01) @ rot2r4(datarigrs3dt01rs3[0][3:6])),rot2r4(datarigrs3dt01rs3[a][3:6])) - 1.

#     counter += 1


# counter = 0
# print("Collating dt = .001 Data")
# for a in range(np.shape(datarigrh3dt001rs2)[0]):
#     h3diffdatarrads2dt001raw[counter][0] = h3distb2((boostxh3(a*.001) @ rot2hyp(datarigrh3dt001rs2[0][0:3])),rot2hyp(datarigrh3dt001rs2[a][0:3])) - 1.
#     h3diffdatarrads2dt001raw[counter][1] = h3distb2((boostxh3(a*.001) @ rot2hyp(datarigrh3dt001rs2[0][3:6])),rot2hyp(datarigrh3dt001rs2[a][3:6])) - 1.
#     h3diffdatarrads3dt001raw[counter][0] = h3distb2((boostxh3(a*.001) @ rot2hyp(datarigrh3dt001rs3[0][0:3])),rot2hyp(datarigrh3dt001rs3[a][0:3])) - 1.
#     h3diffdatarrads3dt001raw[counter][1] = h3distb2((boostxh3(a*.001) @ rot2hyp(datarigrh3dt001rs3[0][3:6])),rot2hyp(datarigrh3dt001rs3[a][3:6])) - 1.

#     s3diffdatarrads2dt001raw[counter][0] = r4distb2((rotzh3(a*.001) @ rot2r4(datarigrs3dt001rs2[0][0:3])),rot2r4(datarigrs3dt001rs2[a][0:3])) - 1.
#     s3diffdatarrads2dt001raw[counter][1] = r4distb2((rotzh3(a*.001) @ rot2r4(datarigrs3dt001rs2[0][3:6])),rot2r4(datarigrs3dt001rs2[a][3:6])) - 1.
#     s3diffdatarrads3dt001raw[counter][0] = r4distb2((rotzh3(a*.001) @ rot2r4(datarigrs3dt001rs3[0][0:3])),rot2r4(datarigrs3dt001rs3[a][0:3])) - 1.
#     s3diffdatarrads3dt001raw[counter][1] = r4distb2((rotzh3(a*.001) @ rot2r4(datarigrs3dt001rs3[0][3:6])),rot2r4(datarigrs3dt001rs3[a][3:6])) - 1.

#     counter += 1

# fig,ax=plt.subplots(1,1)

# # ax.plot(t_arr001[1:],s3lamdatarrads2dt001raw[1:],'r',label = "Spherical rs2")
# # ax.plot(t_arr1[1:],s3lamdatarrads3dt1raw[1:],'g',label = "Spherical rs3")
# # ax.plot(t_arr1[1:],h3diffdatarrads2dt1raw[1:],'k',label = "Hyperbolic rs2")
# # ax.plot(t_arr001[:],h3diffdatarrads2dt001raw[:,0],'b',label = "Hyperbolic rs2")
# # ax.plot(t_arr001[:],h3diffdatarrads2dt001raw[:,1],'r',label = "Hyperbolic rs2")
# # ax.plot(t_arr001[:],h3diffdatarrads3dt001raw[:,0],'g',label = "Hyperbolic rs3")
# # ax.plot(t_arr001[:],h3diffdatarrads3dt001raw[:,1],'k',label = "Hyperbolic rs3")
# # ax.plot(t_arr001[:],s3diffdatarrads2dt001raw[:,0],'b',label = "Spherical rs2")
# # ax.plot(t_arr001[:],s3diffdatarrads2dt001raw[:,1],'r',label = "Spherical rs2")
# # ax.plot(t_arr001[:],s3diffdatarrads3dt001raw[:,0],'g',label = "Spherical rs3")
# # ax.plot(t_arr001[:],s3diffdatarrads3dt001raw[:,1],'k',label = "Spherical rs3")

# # Correct Plotting - Only for 1 particles

# # ax.plot(t_arr001[:],h3diffdatarrads2dt001raw[:,0],'r',label = "H3 rod rs2 dt=.001")
# # ax.plot(t_arr01[:],h3diffdatarrads2dt01raw[:,0],'k',label = "H3 rod rs2 dt=.01")
# # ax.plot(t_arr1[:],h3diffdatarrads2dt1raw[:,0],'b',label = "H3 rod rs2 dt=.1")
# # ax.plot(t_arr001[:],h3diffdatarrads3dt001raw[:,0],'r',label = "H3 rod rs3 dt=.001")
# # ax.plot(t_arr01[:],h3diffdatarrads3dt01raw[:,0],'k',label = "H3 rod rs3 dt=.01")
# # ax.plot(t_arr1[:],h3diffdatarrads3dt1raw[:,0],'b',label = "H3 rod rs3 dt=.1")

# # ax.plot(t_arr001[:],s3diffdatarrads2dt001raw[:,0],'r',label = "S3 rod rs2 dt=.001")
# # ax.plot(t_arr01[:],s3diffdatarrads2dt01raw[:,0],'k',label = "S3 rod rs2 dt=.01")
# # ax.plot(t_arr1[:],s3diffdatarrads2dt1raw[:,0],'b',label = "S3 rod rs2 dt=.1")
# ax.plot(t_arr001[:],s3diffdatarrads3dt001raw[:,0],'r',label = "S3 rod rs3 dt=.001")
# ax.plot(t_arr01[:],s3diffdatarrads3dt01raw[:,0],'k',label = "S3 rod rs3 dt=.01")
# ax.plot(t_arr1[:],s3diffdatarrads3dt1raw[:,0],'b',label = "S3 rod rs3 dt=.1")



# # ax.plot(t_arr1[1:],datarigrs3dt1rs3[1:],'k',label = "Spherical rs3")

# ax.legend()
# ax.set_title('Percent Difference Error vs. Time')
# ax.set_xlabel('Time')
# ax.set_ylabel('Percent Difference Error')
# plt.show()




#####
# Rigid Rod Lagrange Multiplier plot for H3 and S3
#####

# def h3transD12(a1, b1, g1, a2, b2, g2):
#     return -sinh(a1)*sinh(a2) + cosh(a1)*cosh(a2)*(cosh(b1)*cosh(b2)*cosh(g1 - g2) - sinh(b1)*sinh(b2))

# # First Derivatives
# def h3transda1D12(a1, b1, g1, a2, b2, g2):
#     return -cosh(a1)*sinh(a2) + cosh(a2)*sinh(a1)*(cosh(b1)*cosh(b2)*cosh(g1 - g2) - sinh(b1)*sinh(b2))

# def s3D12(a1, b1, g1, a2, b2, g2):
#     return cos(a1)*cos(a2) + sin(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

# # First Derivatives
# def s3da1D12(a1, b1, g1, a2, b2, g2):
#     return -cos(a2)*sin(a1) + cos(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

# ha1,hb1,hg1 = rot2transh3pos(datarigrh3dt1rs2[0][0:3])
# ha2,hb2,hg2 = rot2transh3pos(datarigrh3dt1rs2[0][3:6])
# had1,hbd1,hgd1 = rot2transh3vel(datarigrh3dt1rs2[0][0:3],datarigrh3dt1rs2[0][6:9])
# had2,hbd2,hgd2 = rot2transh3vel(datarigrh3dt1rs2[0][3:6],datarigrh3dt1rs2[0][9:12])
# sa1,sb1,sg1,sa2,sb2,sg2,sad1,sbd1,sgd1,sad2,sbd2,sgd2,_ = datarigrs3dt1rs2[0]
# v,x,m = 1.,1.,1.

# # Distance Function
# hd12 = h3transD12(ha1, hb1, hg1, ha2, hb2, hg2)
# # First derivatives of distance function
# hda1d12 = h3transda1D12(ha1, hb1, hg1, ha2, hb2, hg2)

# # Distance Function
# sd12 = s3D12(sa1, sb1, sg1, sa2, sb2, sg2)
# # First derivatives of distance function
# sda1d12 = s3da1D12(sa1, sb1, sg1, sa2, sb2, sg2)


# ### Radau s2

# h3lamdatarrads2dt1raw   = np.zeros(np.shape(datarigrh3dt1rs2)[0])
# h3lamdatarrads2dt01raw  = np.zeros(np.shape(datarigrh3dt01rs2)[0])
# h3lamdatarrads2dt001raw = np.zeros(np.shape(datarigrh3dt001rs2)[0])

# s3lamdatarrads2dt1raw   = np.zeros(np.shape(datarigrs3dt1rs2)[0])
# s3lamdatarrads2dt01raw  = np.zeros(np.shape(datarigrs3dt01rs2)[0])
# s3lamdatarrads2dt001raw = np.zeros(np.shape(datarigrs3dt001rs2)[0])

# ### Radau s3

# h3lamdatarrads3dt1raw   = np.zeros(np.shape(datarigrh3dt1rs3)[0])
# h3lamdatarrads3dt01raw  = np.zeros(np.shape(datarigrh3dt01rs3)[0])
# h3lamdatarrads3dt001raw = np.zeros(np.shape(datarigrh3dt001rs3)[0])

# s3lamdatarrads3dt1raw   = np.zeros(np.shape(datarigrs3dt1rs3)[0])
# s3lamdatarrads3dt01raw  = np.zeros(np.shape(datarigrs3dt01rs3)[0])
# s3lamdatarrads3dt001raw = np.zeros(np.shape(datarigrs3dt001rs3)[0])

# h3anal_value = 1./2.*sinh(2. * ha1)*hbd1**2./(hda1d12/(m*sqrt(hd12**2. - 1.)))
# s3anal_value = -1./2.*sin(2. * sa1)*sgd1**2./(sda1d12/(m*sqrt(1. - sd12**2.)))

# counter = 0
# print("Collating dt = .1 Data")
# for a in range(np.shape(datarigrh3dt1rs2)[0]):
#     h3lamdatarrads2dt1raw[counter] = (h3anal_value - datarigrh3dt1rs2[a][-1])/h3anal_value
#     h3lamdatarrads3dt1raw[counter] = (h3anal_value - datarigrh3dt1rs3[a][-1])/h3anal_value

#     s3lamdatarrads2dt1raw[counter] = (s3anal_value - datarigrs3dt1rs2[a][-1])/s3anal_value
#     s3lamdatarrads3dt1raw[counter] = (s3anal_value - datarigrs3dt1rs3[a][-1])/s3anal_value

#     counter += 1


# counter = 0
# print("Collating dt = .01 Data")
# for a in range(np.shape(datarigrh3dt01rs2)[0]):
#     h3lamdatarrads2dt01raw[counter] = (h3anal_value - datarigrh3dt01rs2[a][-1])/h3anal_value
#     h3lamdatarrads3dt01raw[counter] = (h3anal_value - datarigrh3dt01rs3[a][-1])/h3anal_value

#     s3lamdatarrads2dt01raw[counter] = (s3anal_value - datarigrs3dt01rs2[a][-1])/s3anal_value
#     s3lamdatarrads3dt01raw[counter] = (s3anal_value - datarigrs3dt01rs3[a][-1])/s3anal_value

#     counter += 1


# counter = 0
# print("Collating dt = .001 Data")
# for a in range(np.shape(datarigrh3dt001rs2)[0]):
#     h3lamdatarrads2dt001raw[counter] = (h3anal_value - datarigrh3dt001rs2[a][-1])/h3anal_value
#     h3lamdatarrads3dt001raw[counter] = (h3anal_value - datarigrh3dt001rs3[a][-1])/h3anal_value

#     s3lamdatarrads2dt001raw[counter] = (s3anal_value - datarigrs3dt001rs2[a][-1])/s3anal_value
#     s3lamdatarrads3dt001raw[counter] = (s3anal_value - datarigrs3dt001rs3[a][-1])/s3anal_value

#     counter += 1

# timevals = [.1,.01,.001]

# h3errvals2 = [max(abs(h3lamdatarrads2dt1raw[1:])),max(abs(h3lamdatarrads2dt01raw[1:])),max(abs(h3lamdatarrads2dt001raw[1:]))]
# h3errvals3 = [max(abs(h3lamdatarrads3dt1raw[1:])),max(abs(h3lamdatarrads3dt01raw[1:])),max(abs(h3lamdatarrads3dt001raw[1:]))]

# s3errvals2 = [max(abs(s3lamdatarrads2dt1raw[1:])),max(abs(s3lamdatarrads2dt01raw[1:])),max(abs(s3lamdatarrads2dt001raw[1:]))]
# s3errvals3 = [max(abs(s3lamdatarrads3dt1raw[1:])),max(abs(s3lamdatarrads3dt01raw[1:])),max(abs(s3lamdatarrads3dt001raw[1:]))]

# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# # ax.plot(t_arr1[:],h3lamdatarrads2dt1raw[:],'r',label = r"$\mathbf{H}^3$ rs2 dt=.1")
# # ax.plot(t_arr01[:],h3lamdatarrads2dt01raw[:],'b',label = r"$\mathbf{H}^3$ rs2 dt=.01")
# # ax.plot(t_arr001[:],h3lamdatarrads2dt001raw[:],'k',label = r"$\mathbf{H}^3$ rs2 dt=.001")

# # ax.plot(t_arr1[:],h3lamdatarrads3dt1raw[:],'r',label = r"$\mathbf{H}^3$ rs3 dt=.1")
# # ax.plot(t_arr01[:],h3lamdatarrads3dt01raw[:],'b',label = r"$\mathbf{H}^3$ rs3 dt=.01")
# # ax.plot(t_arr001[:],h3lamdatarrads3dt001raw[:],'k',label = r"$\mathbf{H}^3$ rs3 dt=.001")

# # ax.plot(t_arr1[:],s3lamdatarrads2dt1raw[:],'r',label = r"$\mathbf{S}^3$ rs2 dt=.1")
# # ax.plot(t_arr01[:],s3lamdatarrads2dt01raw[:],'b',label = r"$\mathbf{S}^3$ rs2 dt=.01")
# # ax.plot(t_arr001[:],s3lamdatarrads2dt001raw[:],'k',label = r"$\mathbf{S}^3$ rs2 dt=.001")

# ax.plot(t_arr1[:],s3lamdatarrads3dt1raw[:],'r',label = r"$\mathbf{S}^3$ rs3 dt=.1")
# ax.plot(t_arr01[:],s3lamdatarrads3dt01raw[:],'b',label = r"$\mathbf{S}^3$ rs3 dt=.01")
# ax.plot(t_arr001[:],s3lamdatarrads3dt001raw[:],'k',label = r"$\mathbf{S}^3$ rs3 dt=.001")

# # # ax.plot(t_arr001[1:],h3lamdatarrads3dt001raw[1:],'r',label = "H3 Rod rs3 dt=.001")
# # # ax.plot(t_arr01[1:],h3lamdatarrads3dt01raw[1:],'k',label = "H3 Rod rs3 dt=.01")
# # # ax.plot(t_arr1[1:],h3lamdatarrads3dt1raw[1:],'b',label = "H3 Rod rs3 dt=.1")

# # ax.scatter(timevals,h3errvals12s1,marker = 's', color='r',label = r"$\mathbf{H}^3$ gs1 s12")
# # ax.scatter(timevals,h3errvals12s2,marker = 's', color='b',label = r"$\mathbf{H}^3$ gs2 s12")
# # ax.scatter(timevals,h3errvals12s3,marker = 's', color='k',label = r"$\mathbf{H}^3$ gs3 s12")


# # for h data rs3
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-2.0}'
# # ylabels[1] = r'\textbf{-1.0}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{1.0}'
# # ylabels[4] = r'\textbf{2.0}'
# # ylabels[5] = r'\textbf{3.0}'

# # # for s data rs2
# ylabels = [item.get_text() for item in ax.get_yticklabels()]
# ylabels[0] = r'\textbf{0.0}'
# ylabels[1] = r'\textbf{0.0}'
# ylabels[2] = r'\textbf{2.0}'
# ylabels[3] = r'\textbf{4.0}'
# ylabels[4] = r'\textbf{6.0}'
# ylabels[5] = r'\textbf{8.0}'


# # ax.set_xticklabels(xlabels)
# # ax.set_yticklabels(ylabels)
# ax.legend(fontsize="20",loc="upper left")
# # plt.title(r'\textbf{ Rigid Rod Body in $\mathbf{H}^3$}', fontsize=20)
# plt.title(r'\textbf{ Rigid Rod Body in $\mathbf{S}^3$}', fontsize=20)
# plt.ylabel(r'\textbf{Relative Error of \\ \\ Lagrange Multiplier}', fontsize=20, labelpad = 20)
# # plt.ylabel(r'\textbf{Relative Error of \\ \\ Lagrange Multiplier ($\times 10^{-7}$)}', fontsize=20, labelpad = 20)
# # plt.ylabel(r'\textbf{Relative Error of \\ \\ Lagrange Multiplier ($\times 10^{-10}$)}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)


# plt.tight_layout()

# plt.show()


###### old stuff


# fig,ax=plt.subplots(1,1)


# ax.plot(t_arr001[1:],s3lamdatarrads2dt001raw[1:],'r',label = "S3 Rod rs2 dt=.001")
# ax.plot(t_arr01[1:],s3lamdatarrads2dt01raw[1:],'k',label = "S3 Rod rs2 dt=.01")
# ax.plot(t_arr1[1:],s3lamdatarrads2dt1raw[1:],'b',label = "S3 Rod rs2 dt=.1")
# # ax.plot(t_arr001[1:],s3lamdatarrads3dt001raw[1:],'r',label = "S3 Rod rs3 dt=.001")
# # ax.plot(t_arr01[1:],s3lamdatarrads3dt01raw[1:],'k',label = "S3 Rod rs3 dt=.01")
# # ax.plot(t_arr1[1:],s3lamdatarrads3dt1raw[1:],'b',label = "S3 Rod rs3 dt=.1")

# # ax.plot(t_arr001[1:],h3lamdatarrads2dt001raw[1:],'r',label = "H3 Rod rs2 dt=.001")
# # ax.plot(t_arr01[1:],h3lamdatarrads2dt01raw[1:],'k',label = "H3 Rod rs2 dt=.01")
# # ax.plot(t_arr1[1:],h3lamdatarrads2dt1raw[1:],'b',label = "H3 Rod rs2 dt=.1")
# # ax.plot(t_arr001[1:],h3lamdatarrads3dt001raw[1:],'r',label = "H3 Rod rs3 dt=.001")
# # ax.plot(t_arr01[1:],h3lamdatarrads3dt01raw[1:],'k',label = "H3 Rod rs3 dt=.01")
# # ax.plot(t_arr1[1:],h3lamdatarrads3dt1raw[1:],'b',label = "H3 Rod rs3 dt=.1")
# ax.legend()
# ax.set_title('Relative Error of Lagrange Multiplier vs. Time')
# ax.set_xlabel('Time')
# ax.set_ylabel('Relative Error of Lagrange Multiplier')
# plt.show()



#######################################################

######### Rigid Triangle plots for H3 and S3 ##########

#######################################################


#####
# Rigid Triangle Lagrange Multiplier plot for H3 and S3
#####

# def h3transD12(a1, b1, g1, a2, b2, g2):
#     return -sinh(a1)*sinh(a2) + cosh(a1)*cosh(a2)*(cosh(b1)*cosh(b2)*cosh(g1 - g2) - sinh(b1)*sinh(b2))

# # First Derivatives
# def h3transda1D12(a1, b1, g1, a2, b2, g2):
#     return -cosh(a1)*sinh(a2) + cosh(a2)*sinh(a1)*(cosh(b1)*cosh(b2)*cosh(g1 - g2) - sinh(b1)*sinh(b2))

# def h3transdb1D12(a1, b1, g1, a2, b2, g2):
#     return cosh(a1)*cosh(a2)*(cosh(b2)*cosh(g1 - g2)*sinh(b1) - cosh(b1)*sinh(b2))



# def s3D12(a1, b1, g1, a2, b2, g2):
#     return cos(a1)*cos(a2) + sin(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

# # First Derivatives
# def s3da1D12(a1, b1, g1, a2, b2, g2):
#     return -cos(a2)*sin(a1) + cos(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

# def s3db1D12(a1, b1, g1, a2, b2, g2):
#         return sin(a1)*sin(a2)*(-cos(b2)*sin(b1) + cos(b1)*cos(g1 - g2)*sin(b2))

# ha1,hb1,hg1 = rot2transh3pos(datarigth3dt1rs2[0][0:3])
# ha2,hb2,hg2 = rot2transh3pos(datarigth3dt1rs2[0][3:6])
# ha3,hb3,hg3 = rot2transh3pos(datarigth3dt1rs2[0][6:9])
# had1,hbd1,hgd1 = rot2transh3vel(datarigth3dt1rs2[0][0:3],datarigth3dt1rs2[0][9:12])
# had2,hbd2,hgd2 = rot2transh3vel(datarigth3dt1rs2[0][3:6],datarigth3dt1rs2[0][12:15])
# had3,hbd3,hgd3 = rot2transh3vel(datarigth3dt1rs2[0][6:9],datarigth3dt1rs2[0][15:18])
# sa1,sb1,sg1,sa2,sb2,sg2,sa3,sb3,sg3,sad1,sbd1,sgd1,sad2,sbd2,sgd2,sad3,sbd3,sgd3,_,_,_ = datarigts3dt1rs2[0]
# v,x,m = 1.,1.,1.

# # Distance Function
# hd12 = h3transD12(ha1, hb1, hg1, ha2, hb2, hg2)
# hd13 = h3transD12(ha1, hb1, hg1, ha3, hb3, hg3)
# hd23 = h3transD12(ha2, hb2, hg2, ha3, hb3, hg3)
# # First derivatives of distance function
# # side 12
# hda1d12 = h3transda1D12(ha1, hb1, hg1, ha2, hb2, hg2)
# hdb1d12 = h3transdb1D12(ha1, hb1, hg1, ha2, hb2, hg2)

# hda2d12 = h3transda1D12(ha2, hb2, hg2, ha1, hb1, hg1)
# hdb2d12 = h3transdb1D12(ha2, hb2, hg2, ha1, hb1, hg1)

# # side 13
# hda1d13 = h3transda1D12(ha1, hb1, hg1, ha3, hb3, hg3)
# hdb1d13 = h3transdb1D12(ha1, hb1, hg1, ha3, hb3, hg3)

# hda3d13 = h3transda1D12(ha3, hb3, hg3, ha1, hb1, hg1)
# hdb3d13 = h3transdb1D12(ha3, hb3, hg3, ha1, hb1, hg1)

# # side 23
# hda2d23 = h3transda1D12(ha2, hb2, hg2, ha3, hb3, hg3)
# hdb2d23 = h3transdb1D12(ha2, hb2, hg2, ha3, hb3, hg3)

# hda3d23 = h3transda1D12(ha3, hb3, hg3, ha2, hb2, hg2)
# hdb3d23 = h3transdb1D12(ha3, hb3, hg3, ha2, hb2, hg2)


# # Distance Function
# sd12 = s3D12(sa1, sb1, sg1, sa2, sb2, sg2)
# sd13 = s3D12(sa1, sb1, sg1, sa3, sb3, sg3)
# sd23 = s3D12(sa2, sb2, sg2, sa3, sb3, sg3)
# # First derivatives of distance function
# # side 12
# sda1d12 = s3da1D12(sa1, sb1, sg1, sa2, sb2, sg2)
# sdb1d12 = s3db1D12(sa1, sb1, sg1, sa2, sb2, sg2)

# sda2d12 = s3da1D12(sa2, sb2, sg2, sa1, sb1, sg1)
# sdb2d12 = s3db1D12(sa2, sb2, sg2, sa1, sb1, sg1)

# # side 13
# sda1d13 = s3da1D12(sa1, sb1, sg1, sa3, sb3, sg3)
# sdb1d13 = s3db1D12(sa1, sb1, sg1, sa3, sb3, sg3)

# sda3d13 = s3da1D12(sa3, sb3, sg3, sa1, sb1, sg1)
# sdb3d13 = s3db1D12(sa3, sb3, sg3, sa1, sb1, sg1)

# # side 23
# sda2d23 = s3da1D12(sa2, sb2, sg2, sa3, sb3, sg3)
# sdb2d23 = s3db1D12(sa2, sb2, sg2, sa3, sb3, sg3)

# sda3d23 = s3da1D12(sa3, sb3, sg3, sa2, sb2, sg2)
# sdb3d23 = s3db1D12(sa3, sb3, sg3, sa2, sb2, sg2)


# ### Radau s2

# h3lam1datatrads2dt1raw   = np.zeros(np.shape(datarigth3dt1rs2)[0])
# h3lam1datatrads2dt01raw  = np.zeros(np.shape(datarigth3dt01rs2)[0])
# h3lam1datatrads2dt001raw = np.zeros(np.shape(datarigth3dt001rs2)[0])
# h3lam2datatrads2dt1raw   = np.zeros(np.shape(datarigth3dt1rs2)[0])
# h3lam2datatrads2dt01raw  = np.zeros(np.shape(datarigth3dt01rs2)[0])
# h3lam2datatrads2dt001raw = np.zeros(np.shape(datarigth3dt001rs2)[0])
# h3lam3datatrads2dt1raw   = np.zeros(np.shape(datarigth3dt1rs2)[0])
# h3lam3datatrads2dt01raw  = np.zeros(np.shape(datarigth3dt01rs2)[0])
# h3lam3datatrads2dt001raw = np.zeros(np.shape(datarigth3dt001rs2)[0])

# s3lam1datatrads2dt1raw   = np.zeros(np.shape(datarigts3dt1rs2)[0])
# s3lam1datatrads2dt01raw  = np.zeros(np.shape(datarigts3dt01rs2)[0])
# s3lam1datatrads2dt001raw = np.zeros(np.shape(datarigts3dt001rs2)[0])
# s3lam2datatrads2dt1raw   = np.zeros(np.shape(datarigts3dt1rs2)[0])
# s3lam2datatrads2dt01raw  = np.zeros(np.shape(datarigts3dt01rs2)[0])
# s3lam2datatrads2dt001raw = np.zeros(np.shape(datarigts3dt001rs2)[0])
# s3lam3datatrads2dt1raw   = np.zeros(np.shape(datarigts3dt1rs2)[0])
# s3lam3datatrads2dt01raw  = np.zeros(np.shape(datarigts3dt01rs2)[0])
# s3lam3datatrads2dt001raw = np.zeros(np.shape(datarigts3dt001rs2)[0])

# ### Radau s3

# h3lam1datatrads3dt1raw   = np.zeros(np.shape(datarigth3dt1rs3)[0])
# h3lam1datatrads3dt01raw  = np.zeros(np.shape(datarigth3dt01rs3)[0])
# h3lam1datatrads3dt001raw = np.zeros(np.shape(datarigth3dt001rs3)[0])
# h3lam2datatrads3dt1raw   = np.zeros(np.shape(datarigth3dt1rs3)[0])
# h3lam2datatrads3dt01raw  = np.zeros(np.shape(datarigth3dt01rs3)[0])
# h3lam2datatrads3dt001raw = np.zeros(np.shape(datarigth3dt001rs3)[0])
# h3lam3datatrads3dt1raw   = np.zeros(np.shape(datarigth3dt1rs3)[0])
# h3lam3datatrads3dt01raw  = np.zeros(np.shape(datarigth3dt01rs3)[0])
# h3lam3datatrads3dt001raw = np.zeros(np.shape(datarigth3dt001rs3)[0])

# s3lam1datatrads3dt1raw   = np.zeros(np.shape(datarigts3dt1rs3)[0])
# s3lam1datatrads3dt01raw  = np.zeros(np.shape(datarigts3dt01rs3)[0])
# s3lam1datatrads3dt001raw = np.zeros(np.shape(datarigts3dt001rs3)[0])
# s3lam2datatrads3dt1raw   = np.zeros(np.shape(datarigts3dt1rs3)[0])
# s3lam2datatrads3dt01raw  = np.zeros(np.shape(datarigts3dt01rs3)[0])
# s3lam2datatrads3dt001raw = np.zeros(np.shape(datarigts3dt001rs3)[0])
# s3lam3datatrads3dt1raw   = np.zeros(np.shape(datarigts3dt1rs3)[0])
# s3lam3datatrads3dt01raw  = np.zeros(np.shape(datarigts3dt01rs3)[0])
# s3lam3datatrads3dt001raw = np.zeros(np.shape(datarigts3dt001rs3)[0])

# hlam1coa2 = (hda2d12/(m*sqrt(hd12**2. - 1.)))
# hlam1cob2 = (hdb2d12/(m*cosh(ha2)**2.*sqrt(hd12**2. - 1.)))

# hlam2coa3 = (hda3d13/(m*sqrt(hd13**2. - 1.)))
# hlam2cob3 = (hdb3d13/(m*cosh(ha3)**2.*sqrt(hd13**2. - 1.)))

# hlam3coa2 = (hda2d23/(m*sqrt(hd23**2. - 1.)))
# hlam3cob2 = (hdb2d23/(m*cosh(ha2)**2.*sqrt(hd23**2. - 1.)))
# hlam3coa3 = (hda3d23/(m*sqrt(hd23**2. - 1.)))
# hlam3cob3 = (hdb3d23/(m*cosh(ha3)**2.*sqrt(hd23**2. - 1.)))

# h3geoa2 = 1./2.*sinh(2. * ha2)*hbd2**2.
# h3geob2 = -2.*had2*hbd2*tanh(ha2)
# h3geoa3 = 1./2.*sinh(2. * ha3)*hbd3**2.
# h3geob3 = -2.*had3*hbd3*tanh(ha3)

# h3anal_values12 = (hlam3coa2*h3geob2 - hlam3cob2*h3geoa2)/(hlam1cob2*hlam3coa2 - hlam1coa2*hlam3cob2)
# h3anal_values13 = (hlam3coa2*h3geob3 - hlam3cob2*h3geoa3)/(hlam2cob3*hlam3coa2 - hlam2coa3*hlam3cob2)
# h3anal_values23 = (hlam1cob2*h3geoa2 - hlam1coa2*h3geob2)/(hlam1cob2*hlam3coa2 - hlam1coa2*hlam3cob2)


# slam1coa2 = -(sda2d12/(m*sqrt(1. - sd12**2.)))
# slam1cob2 = -(sdb2d12/(m*sin(sa2)**2.*sqrt(1. - sd12**2.)))

# slam2coa3 = -(sda3d13/(m*sqrt(1. - sd13**2.)))
# slam2cob3 = -(sdb3d13/(m*sin(sa3)**2.*sqrt(1. - sd13**2.)))

# slam3coa2 = -(sda2d23/(m*sqrt(1. - sd23**2.)))
# slam3cob2 = -(sdb2d23/(m*sin(sa2)**2.*sqrt(1. - sd23**2.)))
# slam3coa3 = -(sda3d23/(m*sqrt(1. - sd23**2.)))
# slam3cob3 = -(sdb3d23/(m*sin(sa3)**2.*sqrt(1. - sd23**2.)))

# s3geoa2 = 1./2.*sin(2. * sa2)*sgd2**2.
# s3geob2 = -2.*sad2*sgd2/tan(sa2)
# s3geoa3 = 1./2.*sin(2. * sa3)*sgd3**2.
# s3geob3 = -2.*sad3*sgd3/tan(sa3)

# s3anal_values12 = (slam3coa2*s3geob2 - slam3cob2*s3geoa2)/(slam1cob2*slam3coa2 - slam1coa2*slam3cob2)
# s3anal_values13 = (slam3coa2*s3geob3 - slam3cob2*s3geoa3)/(slam2cob3*slam3coa2 - slam2coa3*slam3cob2)
# s3anal_values23 = (slam1cob2*s3geoa2 - slam1coa2*s3geob2)/(slam1cob2*slam3coa2 - slam1coa2*slam3cob2)

# counter = 0
# print("Collating dt = .1 Data")
# for a in range(np.shape(datarigth3dt1rs2)[0]):
#     h3lam1datatrads2dt1raw[counter] = (h3anal_values12 - datarigth3dt1rs2[a][-3])#/h3anal_values12
#     h3lam1datatrads3dt1raw[counter] = (h3anal_values12 - datarigth3dt1rs3[a][-3])#/h3anal_values12
#     h3lam2datatrads2dt1raw[counter] = (h3anal_values13 - datarigth3dt1rs2[a][-2])#/h3anal_values13
#     h3lam2datatrads3dt1raw[counter] = (h3anal_values13 - datarigth3dt1rs3[a][-2])#/h3anal_values13
#     h3lam3datatrads2dt1raw[counter] = (h3anal_values23 - datarigth3dt1rs2[a][-1])#/h3anal_values23
#     h3lam3datatrads3dt1raw[counter] = (h3anal_values23 - datarigth3dt1rs3[a][-1])#/h3anal_values23

#     s3lam1datatrads2dt1raw[counter] = (s3anal_values12 - datarigts3dt1rs2[a][-3])#/s3anal_values12
#     s3lam1datatrads3dt1raw[counter] = (s3anal_values12 - datarigts3dt1rs3[a][-3])#/s3anal_values12
#     s3lam2datatrads2dt1raw[counter] = (s3anal_values13 - datarigts3dt1rs2[a][-2])#/s3anal_values13
#     s3lam2datatrads3dt1raw[counter] = (s3anal_values13 - datarigts3dt1rs3[a][-2])#/s3anal_values13
#     s3lam3datatrads2dt1raw[counter] = (s3anal_values23 - datarigts3dt1rs2[a][-1])#/s3anal_values23
#     s3lam3datatrads3dt1raw[counter] = (s3anal_values23 - datarigts3dt1rs3[a][-1])#/s3anal_values23

#     counter += 1


# counter = 0
# print("Collating dt = .01 Data")
# for a in range(np.shape(datarigth3dt01rs2)[0]):
#     h3lam1datatrads2dt01raw[counter] = (h3anal_values12 - datarigth3dt01rs2[a][-3])#/h3anal_values12
#     h3lam1datatrads3dt01raw[counter] = (h3anal_values12 - datarigth3dt01rs3[a][-3])#/h3anal_values12
#     h3lam2datatrads2dt01raw[counter] = (h3anal_values13 - datarigth3dt01rs2[a][-2])#/h3anal_values13
#     h3lam2datatrads3dt01raw[counter] = (h3anal_values13 - datarigth3dt01rs3[a][-2])#/h3anal_values13
#     h3lam3datatrads2dt01raw[counter] = (h3anal_values23 - datarigth3dt01rs2[a][-1])#/h3anal_values23
#     h3lam3datatrads3dt01raw[counter] = (h3anal_values23 - datarigth3dt01rs3[a][-1])#/h3anal_values23

#     s3lam1datatrads2dt01raw[counter] = (s3anal_values12 - datarigts3dt01rs2[a][-3])#/s3anal_values12
#     s3lam1datatrads3dt01raw[counter] = (s3anal_values12 - datarigts3dt01rs3[a][-3])#/s3anal_values12
#     s3lam2datatrads2dt01raw[counter] = (s3anal_values13 - datarigts3dt01rs2[a][-2])#/s3anal_values13
#     s3lam2datatrads3dt01raw[counter] = (s3anal_values13 - datarigts3dt01rs3[a][-2])#/s3anal_values13
#     s3lam3datatrads2dt01raw[counter] = (s3anal_values23 - datarigts3dt01rs2[a][-1])#/s3anal_values23
#     s3lam3datatrads3dt01raw[counter] = (s3anal_values23 - datarigts3dt01rs3[a][-1])#/s3anal_values23

#     counter += 1


# counter = 0
# print("Collating dt = .001 Data")
# for a in range(np.shape(datarigth3dt001rs2)[0]):
#     h3lam1datatrads2dt001raw[counter] = (h3anal_values12 - datarigth3dt001rs2[a][-3])#/h3anal_values12
#     h3lam1datatrads3dt001raw[counter] = (h3anal_values12 - datarigth3dt001rs3[a][-3])#/h3anal_values12
#     h3lam2datatrads2dt001raw[counter] = (h3anal_values13 - datarigth3dt001rs2[a][-2])#/h3anal_values13
#     h3lam2datatrads3dt001raw[counter] = (h3anal_values13 - datarigth3dt001rs3[a][-2])#/h3anal_values13
#     h3lam3datatrads2dt001raw[counter] = (h3anal_values23 - datarigth3dt001rs2[a][-1])#/h3anal_values23
#     h3lam3datatrads3dt001raw[counter] = (h3anal_values23 - datarigth3dt001rs3[a][-1])#/h3anal_values23

#     s3lam1datatrads2dt001raw[counter] = (s3anal_values12 - datarigts3dt001rs2[a][-3])#/s3anal_values12
#     s3lam1datatrads3dt001raw[counter] = (s3anal_values12 - datarigts3dt001rs3[a][-3])#/s3anal_values12
#     s3lam2datatrads2dt001raw[counter] = (s3anal_values13 - datarigts3dt001rs2[a][-2])#/s3anal_values13
#     s3lam2datatrads3dt001raw[counter] = (s3anal_values13 - datarigts3dt001rs3[a][-2])#/s3anal_values13
#     s3lam3datatrads2dt001raw[counter] = (s3anal_values23 - datarigts3dt001rs2[a][-1])#/s3anal_values23
#     s3lam3datatrads3dt001raw[counter] = (s3anal_values23 - datarigts3dt001rs3[a][-1])#/s3anal_values23

#     counter += 1

# h3errvals12s2 = [max(abs(h3lam1datatrads2dt1raw[1:])),max(abs(h3lam1datatrads2dt01raw[1:])),max(abs(h3lam1datatrads2dt001raw[1:]))]
# h3errvals13s2 = [max(abs(h3lam2datatrads2dt1raw[1:])),max(abs(h3lam2datatrads2dt01raw[1:])),max(abs(h3lam2datatrads2dt001raw[1:]))]
# h3errvals23s2 = [max(abs(h3lam3datatrads2dt1raw[1:])),max(abs(h3lam3datatrads2dt01raw[1:])),max(abs(h3lam3datatrads2dt001raw[1:]))]
# h3errvals12s3 = [max(abs(h3lam1datatrads3dt1raw[1:])),max(abs(h3lam1datatrads3dt01raw[1:])),max(abs(h3lam1datatrads3dt001raw[1:]))]
# h3errvals13s3 = [max(abs(h3lam2datatrads3dt1raw[1:])),max(abs(h3lam2datatrads3dt01raw[1:])),max(abs(h3lam2datatrads3dt001raw[1:]))]
# h3errvals23s3 = [max(abs(h3lam3datatrads3dt1raw[1:])),max(abs(h3lam3datatrads3dt01raw[1:])),max(abs(h3lam3datatrads3dt001raw[1:]))]

# s3errvals12s2 = [max(abs(s3lam1datatrads2dt1raw[1:])),max(abs(s3lam1datatrads2dt01raw[1:])),max(abs(s3lam1datatrads2dt001raw[1:]))]
# s3errvals13s2 = [max(abs(s3lam2datatrads2dt1raw[1:])),max(abs(s3lam2datatrads2dt01raw[1:])),max(abs(s3lam2datatrads2dt001raw[1:]))]
# s3errvals23s2 = [max(abs(s3lam3datatrads2dt1raw[1:])),max(abs(s3lam3datatrads2dt01raw[1:])),max(abs(s3lam3datatrads2dt001raw[1:]))]
# s3errvals12s3 = [max(abs(s3lam1datatrads3dt1raw[1:])),max(abs(s3lam1datatrads3dt01raw[1:])),max(abs(s3lam1datatrads3dt001raw[1:]))]
# s3errvals13s3 = [max(abs(s3lam2datatrads3dt1raw[1:])),max(abs(s3lam2datatrads3dt01raw[1:])),max(abs(s3lam2datatrads3dt001raw[1:]))]
# s3errvals23s3 = [max(abs(s3lam3datatrads3dt1raw[1:])),max(abs(s3lam3datatrads3dt01raw[1:])),max(abs(s3lam3datatrads3dt001raw[1:]))]

# # s3errvals2 = [max(abs(s3lamdatarrads2dt1raw[1:])),max(abs(s3lamdatarrads2dt01raw[1:])),max(abs(s3lamdatarrads2dt001raw[1:]))]
# # s3errvals3 = [max(abs(s3lamdatarrads3dt1raw[1:])),max(abs(s3lamdatarrads3dt01raw[1:])),max(abs(s3lamdatarrads3dt001raw[1:]))]

# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# #### H3

# # ax.plot(t_arr1[:],h3lam1datatrads2dt1raw[:],'r',label = r"$\mathbf{H}^3$ s12 rs2 dt=.1")
# # ax.plot(t_arr01[:],h3lam1datatrads2dt01raw[:],'b',label = r"$\mathbf{H}^3$ s12 rs2 dt=.01")
# # ax.plot(t_arr001[:],h3lam1datatrads2dt001raw[:],'k',label = r"$\mathbf{H}^3$ s12 rs2 dt=.001")

# # ax.plot(t_arr1[:],h3lam1datatrads3dt1raw[:],'r',label = r"$\mathbf{H}^3$ s12 rs3 dt=.1")
# # ax.plot(t_arr01[:],h3lam1datatrads3dt01raw[:],'b',label = r"$\mathbf{H}^3$ s12 rs3 dt=.01")
# # ax.plot(t_arr001[:],h3lam1datatrads3dt001raw[:],'k',label = r"$\mathbf{H}^3$ s12 rs3 dt=.001")


# # ax.plot(t_arr1[:],h3lam2datatrads2dt1raw[:],'r',label = r"$\mathbf{H}^3$ s13 rs2 dt=.1")
# # ax.plot(t_arr01[:],h3lam2datatrads2dt01raw[:],'b',label = r"$\mathbf{H}^3$ s13 rs2 dt=.01")
# # ax.plot(t_arr001[:],h3lam2datatrads2dt001raw[:],'k',label = r"$\mathbf{H}^3$ s13 rs2 dt=.001")

# # ax.plot(t_arr1[:],h3lam2datatrads3dt1raw[:],'r',label = r"$\mathbf{H}^3$ s13 rs3 dt=.1")
# # ax.plot(t_arr01[:],h3lam2datatrads3dt01raw[:],'b',label = r"$\mathbf{H}^3$ s13 rs3 dt=.01")
# # ax.plot(t_arr001[:],h3lam2datatrads3dt001raw[:],'k',label = r"$\mathbf{H}^3$ s13 rs3 dt=.001")


# # ax.plot(t_arr1[:],h3lam3datatrads2dt1raw[:],'r',label = r"$\mathbf{H}^3$ s23 rs2 dt=.1")
# # ax.plot(t_arr01[:],h3lam3datatrads2dt01raw[:],'b',label = r"$\mathbf{H}^3$ s23 rs2 dt=.01")
# # ax.plot(t_arr001[:],h3lam3datatrads2dt001raw[:],'k',label = r"$\mathbf{H}^3$ s23 rs2 dt=.001")

# # ax.plot(t_arr1[:],h3lam3datatrads3dt1raw[:],'r',label = r"$\mathbf{H}^3$ s23 rs3 dt=.1")
# # ax.plot(t_arr01[:],h3lam3datatrads3dt01raw[:],'b',label = r"$\mathbf{H}^3$ s23 rs3 dt=.01")
# # ax.plot(t_arr001[:],h3lam3datatrads3dt001raw[:],'k',label = r"$\mathbf{H}^3$ s23 rs3 dt=.001")

# #### S3

# # ax.plot(t_arr1[:],s3lam1datatrads2dt1raw[:],'r',label = r"$\mathbf{S}^3$ s12 rs2 dt=.1")
# # ax.plot(t_arr01[:],s3lam1datatrads2dt01raw[:],'b',label = r"$\mathbf{S}^3$ s12 rs2 dt=.01")
# # ax.plot(t_arr001[:],s3lam1datatrads2dt001raw[:],'k',label = r"$\mathbf{S}^3$ s12 rs2 dt=.001")

# # ax.plot(t_arr1[:],s3lam1datatrads3dt1raw[:],'r',label = r"$\mathbf{S}^3$ s12 rs3 dt=.1")
# # ax.plot(t_arr01[:],s3lam1datatrads3dt01raw[:],'b',label = r"$\mathbf{S}^3$ s12 rs3 dt=.01")
# # ax.plot(t_arr001[:],s3lam1datatrads3dt001raw[:],'k',label = r"$\mathbf{S}^3$ s12 rs3 dt=.001")


# # ax.plot(t_arr1[:],s3lam2datatrads2dt1raw[:],'r',label = r"$\mathbf{S}^3$ s13 rs2 dt=.1")
# # ax.plot(t_arr01[:],s3lam2datatrads2dt01raw[:],'b',label = r"$\mathbf{S}^3$ s13 rs2 dt=.01")
# # ax.plot(t_arr001[:],s3lam2datatrads2dt001raw[:],'k',label = r"$\mathbf{S}^3$ s13 rs2 dt=.001")

# # ax.plot(t_arr1[:],s3lam2datatrads3dt1raw[:],'r',label = r"$\mathbf{S}^3$ s13 rs3 dt=.1")
# # ax.plot(t_arr01[:],s3lam2datatrads3dt01raw[:],'b',label = r"$\mathbf{S}^3$ s13 rs3 dt=.01")
# # ax.plot(t_arr001[:],s3lam2datatrads3dt001raw[:],'k',label = r"$\mathbf{S}^3$ s13 rs3 dt=.001")


# # ax.plot(t_arr1[:],s3lam3datatrads2dt1raw[:],'r',label = r"$\mathbf{S}^3$ s23 rs2 dt=.1")
# # ax.plot(t_arr01[:],s3lam3datatrads2dt01raw[:],'b',label = r"$\mathbf{S}^3$ s23 rs2 dt=.01")
# # ax.plot(t_arr001[:],s3lam3datatrads2dt001raw[:],'k',label = r"$\mathbf{S}^3$ s23 rs2 dt=.001")

# ax.plot(t_arr1[:],s3lam3datatrads3dt1raw[:],'r',label = r"$\mathbf{S}^3$ s23 rs3 dt=.1")
# ax.plot(t_arr01[:],s3lam3datatrads3dt01raw[:],'b',label = r"$\mathbf{S}^3$ s23 rs3 dt=.01")
# ax.plot(t_arr001[:],s3lam3datatrads3dt001raw[:],'k',label = r"$\mathbf{S}^3$ s23 rs3 dt=.001")




# # ax.scatter(timevals,h3errvals12s1,marker = 's', color='r',label = r"$\mathbf{H}^3$ gs1 s12")
# # ax.scatter(timevals,h3errvals12s2,marker = 's', color='b',label = r"$\mathbf{H}^3$ gs2 s12")
# # ax.scatter(timevals,h3errvals12s3,marker = 's', color='k',label = r"$\mathbf{H}^3$ gs3 s12")


# # for h data rs3
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-2.0}'
# # ylabels[1] = r'\textbf{-1.0}'
# # ylabels[2] = r'\textbf{-0.5}'
# # ylabels[3] = r'\textbf{0.0}'
# # ylabels[4] = r'\textbf{0.5}'
# # ylabels[5] = r'\textbf{1.0}'

# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-2.0}'
# # ylabels[1] = r'\textbf{-1.0}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{1.0}'
# # ylabels[4] = r'\textbf{0.5}'
# # ylabels[5] = r'\textbf{1.0}'

# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-2.0}'
# # ylabels[1] = r'\textbf{-0.5}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{0.5}'
# # ylabels[4] = r'\textbf{1.0}'
# # ylabels[5] = r'\textbf{1.0}'

# # # for s data rs2
# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-4.0}'
# # ylabels[1] = r'\textbf{-2.0}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{2.0}'
# # ylabels[4] = r'\textbf{4.0}'
# # ylabels[5] = r'\textbf{8.0}'

# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-4.0}'
# # ylabels[1] = r'\textbf{-1.0}'
# # ylabels[2] = r'\textbf{0.0}'
# # ylabels[3] = r'\textbf{1.0}'
# # ylabels[4] = r'\textbf{4.0}'
# # ylabels[5] = r'\textbf{8.0}'

# # ylabels = [item.get_text() for item in ax.get_yticklabels()]
# # ylabels[0] = r'\textbf{-4.0}'
# # ylabels[1] = r'\textbf{-4.0}'
# # ylabels[2] = r'\textbf{-2.0}'
# # ylabels[3] = r'\textbf{0.0}'
# # ylabels[4] = r'\textbf{2.0}'
# # ylabels[5] = r'\textbf{8.0}'

# ylabels = [item.get_text() for item in ax.get_yticklabels()]
# ylabels[0] = r'\textbf{-1.00}'
# ylabels[1] = r'\textbf{-1.00}'
# ylabels[2] = r'\textbf{-0.75}'
# ylabels[3] = r'\textbf{-0.50}'
# ylabels[4] = r'\textbf{-0.25}'
# ylabels[5] = r'\textbf{0.00}'


# # ax.set_xticklabels(xlabels)
# ax.set_yticklabels(ylabels)
# ax.legend(fontsize="15",loc="lower left")
# # plt.title(r'\textbf{ Rigid Triangle Body in $\mathbf{H}^3$}', fontsize=20)
# plt.title(r'\textbf{ Rigid Triangle Body in $\mathbf{S}^3$}', fontsize=20)
# # plt.ylabel(r'\textbf{Absolute Error of \\ \\ Lagrange Multiplier}', fontsize=20, labelpad = 20)
# # plt.ylabel(r'\textbf{Absolute Error of \\ \\ Lagrange Multiplier ($\times 10^{-14}$)}', fontsize=20, labelpad = 20)
# plt.ylabel(r'\textbf{Absolute Error of \\ \\ Lagrange Multiplier ($\times 10^{-13}$)}', fontsize=20, labelpad = 20)
# # plt.ylabel(r'\textbf{Absolute Error of \\ \\ Lagrange Multiplier ($\times 10^{-9}$)}', fontsize=20, labelpad = 20)
# # plt.ylabel(r'\textbf{Absolute Error of \\ \\ Lagrange Multiplier ($\times 10^{-7}$)}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)


# plt.tight_layout()

# plt.show()


####### old stuff

# fig,ax=plt.subplots(1,1)

# # ax.plot(t_arr001[1:],s3lamdatarrads3dt001raw[1:],'r',label = "S3 Rod rs3 dt=.001")
# # ax.plot(t_arr01[1:],s3lamdatarrads3dt01raw[1:],'k',label = "S3 Rod rs3 dt=.01")
# # ax.plot(t_arr1[1:],s3lamdatarrads3dt1raw[1:],'b',label = "S3 Rod rs3 dt=.1")

# # ax.plot(t_arr001[1:],h3lam1datatrads2dt001raw[1:],'k',label = "H3 Tri lam1 rs2 dt=.001")
# # ax.plot(t_arr001[1:],h3lam2datatrads2dt001raw[1:],'b',label = "H3 Tri lam2 rs2 dt=.001")
# # ax.plot(t_arr001[1:],h3lam3datatrads2dt001raw[1:],'r',label = "H3 Tri lam3 rs2 dt=.001")
# # ax.plot(t_arr001[1:],h3lam1datatrads3dt001raw[1:],'k',label = "H3 Tri lam1 rs3 dt=.001")
# # ax.plot(t_arr001[1:],h3lam2datatrads3dt001raw[1:],'b',label = "H3 Tri lam2 rs3 dt=.001")
# # ax.plot(t_arr001[1:],h3lam3datatrads3dt001raw[1:],'r',label = "H3 Tri lam3 rs3 dt=.001")
# # ax.plot(t_arr001[1:],datarigth3dt001rs3[1:,-1],'b',label = "Hyperbolic rs3")

# # ax.plot(t_arr001[1:],s3lam1datatrads2dt001raw[1:],'k',label = "S3 Tri lam1 rs2 dt=.001")
# # ax.plot(t_arr001[1:],s3lam2datatrads2dt001raw[1:],'b',label = "S3 Tri lam2 rs2 dt=.001")
# # ax.plot(t_arr001[1:],s3lam3datatrads2dt001raw[1:],'r',label = "S3 Tri lam3 rs2 dt=.001")
# ax.plot(t_arr001[1:],s3lam1datatrads3dt001raw[1:],'k',label = "S3 Tri lam1 rs3 dt=.001")
# ax.plot(t_arr001[1:],s3lam2datatrads3dt001raw[1:],'b',label = "S3 Tri lam2 rs3 dt=.001")
# ax.plot(t_arr001[1:],s3lam3datatrads3dt001raw[1:],'r',label = "S3 Tri lam3 rs3 dt=.001")
# # ax.plot(t_arr001[1:],datarigts3dt001rs3[1:,-3],'b',label = "Spherical rs3")
# ax.legend()
# ax.set_title('Absolute Error of Lagrange Multipliers vs. Time')
# ax.set_xlabel('Time')
# ax.set_ylabel('Absolute Error of Lagrange Multipliers')
# plt.show()






##################################################

######### Rigid Triangle plots for H2xE ##########

##################################################


#####
# Rigid Triangle len and angle plot for H2xD
#####


# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# #### sep

# # ax.plot(datarigth2edt01s12[:,0],(datarigth2edt01s12[:,1]-datarigth2edt01s12[0,1])/(datarigth2edt01s12[0,1]),'r',label = r"$\mathbf{H}^2 \times \mathbf{E}$ s12")
# # ax.plot(datarigth2edt01s13[:,0],(datarigth2edt01s13[:,1]-datarigth2edt01s13[0,1])/(datarigth2edt01s13[0,1]),'b',label = r"$\mathbf{H}^2 \times \mathbf{E}$ s13")
# # ax.plot(datarigth2edt01s23[:,0],(datarigth2edt01s23[:,1]-datarigth2edt01s23[0,1])/(datarigth2edt01s23[0,1]),'k',label = r"$\mathbf{H}^2 \times \mathbf{E}$ s23")

# #### ang

# ax.plot(datarigth2edt01ang1[:,0],datarigth2edt01ang1[:,1],'r',label = r"$\mathbf{H}^2 \times \mathbf{E}$ ang1")
# ax.plot(datarigth2edt01ang2[:,0],datarigth2edt01ang2[:,1],'b',label = r"$\mathbf{H}^2 \times \mathbf{E}$ ang2")
# ax.plot(datarigth2edt01ang3[:,0],datarigth2edt01ang3[:,1],'k',label = r"$\mathbf{H}^2 \times \mathbf{E}$ ang3")




# ylabels = [item.get_text() for item in ax.get_yticklabels()]
# ylabels[0] = r'\textbf{-1.00}'
# ylabels[1] = r'\textbf{-1.0}'
# ylabels[2] = r'\textbf{-0.5}'
# ylabels[3] = r'\textbf{0.0}'
# ylabels[4] = r'\textbf{0.5}'
# ylabels[5] = r'\textbf{0.00}'


# # ax.set_xticklabels(xlabels)
# # ax.set_yticklabels(ylabels)
# ax.legend(fontsize="15",loc="upper right")
# plt.title(r'\textbf{ Rigid Triangle Body in $\mathbf{H}^2\times\mathbf{E}$}', fontsize=20)
# plt.ylabel(r'\textbf{Vertex Angle}', fontsize=20, labelpad = 20)
# # plt.ylabel(r'\textbf{Relative Error of \\ \\ Vertex Separation ($\times 10^{-15}$)}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)


# plt.tight_layout()

# plt.show()





#######################################

######### Curvature Detector ##########

#######################################


#####
# Rigid Triangle len and angle plot for H2xD
#####

# def bumpgeodist(state_vec):

#     def ppqfunc(a1, b1, g1, a2, b2, g2):
#         return np.sqrt((a1 + a2)**2. + (b1 + b2)**2. + (g1 + g2)**2.)

#     def pnqfunc(a1, b1, g1, a2, b2, g2):
#         return np.sqrt((a1 - a2)**2. + (b1 - b2)**2. + (g1 - g2)**2.)

#     a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec

#     # PPQ & PNQ Functions
#     ppq = ppqfunc(a1, b1, g1, a2, b2, g2)
#     pnq = pnqfunc(a1, b1, g1, a2, b2, g2)

#     fp = np.arcsinh(1/np.sqrt(2)*(ppq + pnq)/2) + (ppq + pnq)/2*np.sqrt(((ppq + pnq)/2)**2 + 2)/2
#     fn = np.arcsinh(1/np.sqrt(2)*(ppq - pnq)/2) + (ppq - pnq)/2*np.sqrt(((ppq - pnq)/2)**2 + 2)/2

#     return (fp - fn)

# dt = .01    # Number of steps
# t_max = 50      # Total simulation time
# t_arr = np.arange(0.,t_max+dt,dt)

# distdata3 = np.zeros(np.shape(databumpcurvdet)[0])

# counter = 0
# for a in range(np.shape(databumpcurvdet)[0]):
#     distdata3[counter] = bumpgeodist(databumpcurvdet[a])
#     counter += 1


# fig,ax=plt.subplots(1,1)
# fig.canvas.draw()

# #### sep

# ax.plot(t_arr,distdata3,'r',label = r"$\mathbf{H}^2 \times \mathbf{E}$ s12")
# ax.plot(t_arr,np.full(t_arr.shape,distdata3[0]),'b',label = r"$\mathbf{H}^2 \times \mathbf{E}$ s13")




# ylabels = [item.get_text() for item in ax.get_yticklabels()]
# ylabels[0] = r'\textbf{-1.00}'
# ylabels[1] = r'\textbf{-1.0}'
# ylabels[2] = r'\textbf{-0.5}'
# ylabels[3] = r'\textbf{0.0}'
# ylabels[4] = r'\textbf{0.5}'
# ylabels[5] = r'\textbf{0.00}'


# # ax.set_xticklabels(xlabels)
# # ax.set_yticklabels(ylabels)
# # ax.legend(fontsize="15",loc="upper right")
# plt.title(r'\textbf{ Curvature Detector in IH space}', fontsize=20)
# plt.ylabel(r'\textbf{Vertex Separation}', fontsize=20, labelpad = 20)
# # plt.ylabel(r'\textbf{Relative Error of \\ \\ Vertex Separation ($\times 10^{-15}$)}', fontsize=20, labelpad = 20)
# plt.xlabel(r'\textbf{Time}', fontsize=20, labelpad = 20)

# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)


# plt.tight_layout()

# plt.show()














