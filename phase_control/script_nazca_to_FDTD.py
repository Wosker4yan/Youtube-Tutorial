# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:09:02 2022

@author: Varahm
"""

#setting up the gds for the directional coupler

from masktopolygon import Masktopolygon
from IPython import get_ipython
get_ipython().magic('reset -sf')
import math
import sys, imp, datetime, os, json
from scipy import constants
import matplotlib.pyplot as plt
import numpy as np
import sys
lumapi = imp.load_source('lumapi', "C:/Program Files/Lumerical/v202/api/python/lumapi.py")
sys.path.append('C:/Users/Vahram/OneDrive/Desktop/Work\LiDAR/Designs/DC/servers_results/Broadband Phase control Section')

cwd      = os.getcwd()
save_path   = cwd

angles =[12.6]
T1 = []
T2 = []
difference = []
for angle in angles:
    import nazca as nd
    with nd.Cell(name='topcell') as topcell:
        #INIT
        w1 = 0.4
        w2= 0.4
        gap=  0.3
        offset_down = 8
        Rc = 26
        Ls = 3.9
        L_in = 2
        offset_up = 6
        R2 =Rc - 0.5*(gap- w1)
        R1 = Rc + 0.5*(gap + w1)
        alfa = angle
        beta = 2*alfa
        length_arm_up = 7
        length_arm_down =5
        FDTD_above=300e-9;
        FDTD_below=300e-9;
        maxvzWAFER = 0.22e-6;
        
        
        
        xs = nd.add_xsection('curved_DC')
        nd.add_layer2xsection('curved_DC', 1)
        ic = nd.interconnects.Interconnect(xs='curved_DC', width = w1)
        
        #UPPER ARM
    
        e1  = ic.bend(angle=alfa, radius=2).put()
        e1.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        strt1 = ic.strt(length =length_arm_up, width=w1).put(e1.pin['a0'])
        e2  = ic.bend(angle=-alfa, radius=R1).put(e1.pin['b0'])
        e2.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        
        get_last_xe2 = e2.pin['b0'].x
        get_last_ye2= e2.pin['b0'].y
        
        e3  = ic.bend(angle=-alfa, radius=R1).put(e2.pin['b0'])
        e3.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        e4  = ic.bend(angle=alfa, radius=2).put(e3.pin['b0'])
        e4.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        strt2 = ic.strt(length =length_arm_up, width=w1).put(e4.pin['b0'])
        
        
        get_x_input_1 = strt1.pin['b0'].x
        get_y_input_1 =strt1.pin['a0'].y
        get_y_output_1 =strt2.pin['b0'].y
        get_x_output_1 = strt2.pin['b0'].x
    
        #lOWER ARM
    
        e5 = ic.bend(angle=-beta, radius=R2).put(get_last_xe2,get_last_ye2-w1-gap)
        e5.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        
        
        e6 = ic.bend(angle=beta, radius=R2).put(e5.pin['a0'])
        e6.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        e7  = ic.bend(angle=-2*alfa, radius=2).put(e6.pin['b0'])
        e7.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        e8  = ic.bend(angle=2*alfa, radius=2).put(e5.pin['b0'])
        e8.raise_pins(['a0', 'b0'], ['a0', 'b0'])
        strt3= ic.strt(length =abs(get_x_input_1-e7.pin['b0'].x), width=w1).put(e7.pin['b0'])
        strt4= ic.strt(length =get_x_output_1-e8.pin['b0'].x, width=w1).put(e8.pin['b0'])
    
        get_down_bend_pin = strt1.pin['a0'].y

     
        
        get_x_input_2 = strt3.pin['b0'].x
        get_x_output_2 =strt4.pin['b0'].x
        get_y_input_2 =strt3.pin['b0'].y
        get_y_output_2 =strt4.pin['b0'].y
        
        
        get_x_min_mesh= strt1.pin['a0'].x
        get_x_max_mesh =strt2.pin['a0'].x
        
        mesh_y_max = get_last_ye2
        mesh_y_min = get_down_bend_pin - w1
        
    nd.export_gds(topcell, 'xxx.gds')
    nd.export_plt(topcell, 'xxx.gds')
    
    cwd  = os.getcwd()
    fdtd = lumapi.FDTD()
    fdtd.switchtolayout;
    fdtd.feval('GDS_2_FDTD');
    
    
    Thickness_Si = maxvzWAFER
    fdtd.setglobalsource("wavelength start",1500e-9) #wavelength start must be lower than wavelength stop
    fdtd.setglobalsource("wavelength stop",1600e-9)
    fdtd.setglobalmonitor("use wavelength spacing",0)
    fdtd.setglobalmonitor("use source limits",1)
    fdtd.setglobalmonitor("frequency points",21)
    
    
    filename = str(alfa)  +"_curved_DC_lms"
    fdtd.save(filename)
    
    
    
    fdtd.addfdtd(  # FDTD simulation volume
    x_min=get_x_input_2*1e-6,
    x_max=get_x_output_2*1e-6,
    y_min=get_y_input_2*1e-6 - w1*1e-6 - 0.5e-6,
    y_max=mesh_y_max*1e-6+1e-6,
    z_min=-FDTD_below,	 
    z_max=maxvzWAFER+FDTD_above,
    mesh_accuracy=3,
    x_min_bc="PML", 
    x_max_bc="PML",
    y_min_bc="PML", 
    y_max_bc="PML",
    z_min_bc="PML", 
    z_max_bc="PML")
    
    fdtd.addmesh(
    x_min=get_x_min_mesh*1e-6, 
    x_max=get_x_max_mesh*1e-6, 
    y_max=(mesh_y_max+w1/2)*1e-6,
    y_min=mesh_y_min*1e-6,
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
    y=0, 
    y_span=w1*2*1e-6, 
    x= (get_x_input_1 + 1.1)*1e-6,
    z_min=-FDTD_below, 
    z_max=maxvzWAFER+FDTD_above)
    
    
    
    #draw monitor input wf
    fdtd.addpower(name = "input",
                  simulation_type = "All",
                  monitor_type = "2D X-normal", 
                  y= get_y_output_1*1e-6,
                  x= get_x_output_1*1e-6-1e-6,
                  y_span = w1*2*1e-6,
                  z_span = maxvzWAFER+FDTD_above/2,
                  z = 0+Thickness_Si/2
                  )
    
    
    
    
    
    #draw monitor output wg
    fdtd.addpower(name = "output",
                  simulation_type = "All",
                  monitor_type = "2D X-normal", 
                  y= get_y_output_2*1e-6,
                  x= get_x_output_1*1e-6-1e-6,
                  y_span = w1*2*1e-6,
                  z_span = maxvzWAFER+FDTD_above/2,
                  z = 0+Thickness_Si/2
                  )


    filename = str(angle)  +"_Curved_DC.lms"
    fdtd.save(save_path + "\\" + filename)
    fdtd.run()
    

    T_input = fdtd.getresult("input","T")
    T_output = fdtd.getresult("output","T")
    inputT = T_input["T"]
    outputT = T_output["T"]
    diff = inputT - outputT
    wavelength_span = T_input["lambda"]
    
    
    
    difference.append(diff)    
    T1.append(inputT)
    T2.append(outputT)
    
    
#    
fig2, ax2 = plt.subplots()
ax2.set_title('Different T1')
for i, trans in enumerate(T2):
    ax2.plot(wavelength_span*10e8, trans, label=angles[i])
            

ax2.set_xlabel('WAVELENGTH [nm]')
ax2.legend()
ax2.set_ylabel('T1')

fig2.tight_layout()
fig2.show()
#

plt.plot(wavelength_span*10e8,np.transpose(T1), label = 'T1')
plt.plot(wavelength_span*10e8,np.transpose(T2), label = 'T2')
plt.title('alfa ='+str(angles[0])+ ' degrees')
plt.xlabel('wavelength [nm]')
plt.ylabel('Transmission')
plt.legend()