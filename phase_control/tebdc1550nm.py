# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 16:20:56 2022

@author: Varahm
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct  19 13:09:02 2022

@author: Varham
"""

#setting up the gds for the directional coupler
import matplotlib.cm as cm

from masktopolygon import Masktopolygon
from IPython import get_ipython
get_ipython().magic('reset -sf')
import math
import importlib.util
from scipy import constants
import matplotlib.pyplot as plt
import numpy as np
import sys
import nazca as nd

spec_win = importlib.util.spec_from_file_location('lumapi', 'C:\\Program Files\\Lumerical\\v241\\api\\python\\lumapi.py')

lumapi = importlib.util.module_from_spec(spec_win) #windows
spec_win.loader.exec_module(lumapi)

T1 = []
T2 = []

FDTD_above=0.75e-6;
FDTD_below=300e-9;

#INITIAL PARAMETERS 
gap = 0.1 #0.19
width1 = 0.2        
width2 =  0.3
width3 = 0.1
L_wg_in_between = 1.5 #4.65
L_taper = 1.00
length_wg1 =2#12.4
length_wg2 =2 #12.4

bend_offset = 0.8
sbends_l1 =1
sbends_l2 =1
radius_bend=20
maxvzWAFER = 0.25e-6

with nd.Cell(name='topcell') as topcell:
    xs = nd.add_xsection('xxx')
    nd.add_layer2xsection('xxx')
    ic = nd.interconnects.Interconnect(xs='xxx', width = width1)
        
    sbend_1 = ic.sbend(radius=radius_bend, width=width1, offset=-sbends_l1,length1=1, length2=1).put()
    sbend_1.raise_pins(['a0', 'b0'], ['a0', 'b0'])   
    
    
    
    strt_up_arm1 = ic.strt(length= length_wg1, width= width1).put(sbend_1.pin['b0'])
    strt_up_arm1.raise_pins(['a0', 'b0'], ['a0', 'b0'])   
    
    
    
    taper1 =ic.taper(width1 = width1, width2 = width3, length=L_taper).put(strt_up_arm1.pin['b0'])
    taper1.raise_pins(['a0', 'b0'], ['a0', 'b0'])  
    
    strt_section = ic.strt(width=width3, length=L_wg_in_between).put(taper1.pin['b0'])
    strt_section.raise_pins(['a0', 'b0'], ['a0', 'b0'])  

    taper2 =ic.taper(width1 = width3, width2 = width1, length=L_taper).put(strt_section.pin['b0'])
    taper2.raise_pins(['a0', 'b0'], ['a0', 'b0'])     
    
        
    strt_up_arm2 = ic.strt(length= length_wg2, width= width1).put(taper2.pin['b0'])
    strt_up_arm2.raise_pins(['a0', 'b0'], ['a0', 'b0'])   
        
    sbend_3 = ic.sbend(radius=radius_bend, width=width1, offset=sbends_l1,length1=1, length2=1).put(strt_up_arm2.pin['b0'])
    sbend_3.raise_pins(['a0', 'b0'], ['a0', 'b0'])   
       
    sbend_2 = ic.sbend(radius=radius_bend, width=width1, offset=sbends_l1,length1=1, length2=1).put(0,-width1-gap-2*sbends_l1)
    sbend_2.raise_pins(['a0', 'b0'], ['a0', 'b0'])   

    strt_up_arm3 = ic.strt(length= length_wg1, width= width1).put(sbend_2.pin['b0'])
    strt_up_arm3.raise_pins(['a0', 'b0'], ['a0', 'b0'])   
    
    
    taper3 =ic.taper(width1 = width1, width2 = width2, length=L_taper).put(strt_up_arm3.pin['b0'])
    taper3.raise_pins(['a0', 'b0'], ['a0', 'b0'])  
    
    strt_section2 = ic.strt(width=width2, length=L_wg_in_between).put(taper3.pin['b0'])
    strt_section2.raise_pins(['a0', 'b0'], ['a0', 'b0'])  
    
    
    taper4 =ic.taper(width1 = width2, width2 = width1, length=L_taper).put(strt_section2.pin['b0'])
    taper4.raise_pins(['a0', 'b0'], ['a0', 'b0'])     
        
    strt_up_arm4 = ic.strt(length= length_wg2, width= width1).put(taper4.pin['b0'])
    strt_up_arm4.raise_pins(['a0', 'b0'], ['a0', 'b0'])   
    
    sbend_4 = ic.sbend(radius=radius_bend, width=width1, offset=-sbends_l1,length1=1, length2=1).put(strt_up_arm4.pin['b0'])
    sbend_4.raise_pins(['a0', 'b0'], ['a0', 'b0'])       
            
    x1 =  sbend_1.pin['a0'].x*1e-6
    y1 =  sbend_1.pin['a0'].y*1e-6  
    
    
        
    ouput_monitor1x = sbend_3.pin['b0'].x*1e-6
    ouput_monitor1y = sbend_3.pin['b0'].y*1e-6

    output_monitor2x =sbend_4.pin['b0'].x*1e-6
    output_monitor2y =sbend_4.pin['b0'].y*1e-6
   
        
    x2 =  sbend_4.pin['b0'].x*1e-6  
    y2 =  sbend_4.pin['b0'].y*1e-6      
        
    mesh_ymax = taper2.pin['b0'].y*1e-6 +width1*0.5*1e-6
    mesh_ymin= taper4.pin['b0'].y*1e-6-width1*1e-6
    
    mesh_xmin= taper1.pin['a0'].x*1e-6 -length_wg1 *0.5*1e-6
    mesh_xmax= taper2.pin['b0'].x*1e-6+length_wg2 *0.5*1e-6

        
nd.export_plt(topcell, 'xxx.gds')
nd.export_gds(topcell, 'xxx.gds')

    
fdtd = lumapi.FDTD()
fdtd.switchtolayout;
fdtd.feval('GDS_2_FDTD');


Thickness_Si = maxvzWAFER
fdtd.setglobalsource("wavelength start",630e-9) #wavelength start must be lower than wavelength stop
fdtd.setglobalsource("wavelength stop",650e-9)
fdtd.setglobalmonitor("use wavelength spacing",0)
fdtd.setglobalmonitor("use source limits",1)
fdtd.setglobalmonitor("frequency points",21)


filename = "Broadband_DC_lms"
fdtd.save(filename)



fdtd.addfdtd(  # FDTD simulation volume
x_min=x1,
x_max=ouput_monitor1x,
y_min=y2 -2*width1*1e-6,
y_max=y1 +2*width1*1e-6,
z_min=-300e-9,	 
z_max=maxvzWAFER+FDTD_above,
mesh_accuracy=3,
x_min_bc="PML", 
x_max_bc="PML",
y_min_bc="PML", 
y_max_bc="PML",
z_min_bc="PML", 
z_max_bc="PML")

fdtd.addmesh(
x_min=mesh_xmin, 
x_max=mesh_xmax, 
y_max=mesh_ymax,
y_min=mesh_ymin,
z_min=0,   
z_max=Thickness_Si,                
override_y_mesh=1, 
override_z_mesh=0, 
override_x_mesh=0,
set_equivalent_index=1,	 
equivalent_y_index=5)



#
fdtd.addmode(
name="source",
injection_axis="x-axis", 
direction="forward", 
y=y2, 
y_span=width1*4.25*1e-6, 
x=0.2e-6,
z=Thickness_Si/2, 
z_span = Thickness_Si*2)



#draw monitor input wf
fdtd.addpower(name = "input",
              simulation_type = "All",
              monitor_type = "2D X-normal", 
              y=ouput_monitor1y,
              x=ouput_monitor1x - 0.5*1e-6,
              y_span =width1*2*1e-6,
              z_span = maxvzWAFER+FDTD_above/2,
              z = 0+Thickness_Si/2
              )





#draw monitor output wg
fdtd.addpower(name = "output",
              simulation_type = "All",
              monitor_type = "2D X-normal", 
              y=output_monitor2y,
              x=output_monitor2x - 0.5*1e-6,
              y_span = width1*2*1e-6,
              z_span = maxvzWAFER+FDTD_above/2,
              z = 0+Thickness_Si/2
              )




#draw field
fdtd.addprofile(name = "profile",
              simulation_type = "All",
              monitor_type = "2D Z-normal", 
              x_min=x1,
              x_max=ouput_monitor1x,
              y_min=y2 -2*width1*1e-6,
              y_max=y1 +2*width1*1e-6,
              z = 0+Thickness_Si/2
              )


save_path   = cwd

filename ="DC.lms"
fdtd.save(save_path + "\\" + filename)
fdtd.run()


T_input = fdtd.getresult("input","T")
T_output = fdtd.getresult("output","T")
inputT = T_input["T"]
outputT = T_output["T"]
diff = inputT - outputT
wavelength_span = T_input["lambda"]

plt.plot(wavelength_span*10e5,inputT , label= 'up')
plt.plot(wavelength_span*10e5,outputT , label= 'down')
plt.xlabel("wavelengths [um]")
plt.title('No SiN on top')
plt.ylabel('Transmission')
plt.yticks(np.linspace(0,1,22))  
plt.legend()
#plt.xticks(np.linspace(1.4, 1.6,11))  
plt.xlim([1.45, 1.6])
