#!/usr/bin/env python

import time
import os
import pandas as pd
import org_radwavefn


# Function to convert k quantum number to j quantum number
def j_qnum(k_quantum_num):
  if abs(k_quantum_num)==1:
    return('1/2')
  elif abs(k_quantum_num)==2:
    return('3/2')
  elif abs(k_quantum_num)==3:
    return('5/2')
  elif abs(k_quantum_num)==4:
    return('7/2')
  elif abs(k_quantum_num)==5:
    return('9/2')

# Function to create Veusz graphs for either numerically defined orbitals or
#   orbital equations.
def create_veusz_file(mode,atom):

  # Define lists to hold legend labels, data series, and names for each
  #  page in the Veusz plot.
  legend_labels_large=[]
  legend_labels_small=[]
  data_series_rad=[]
  data_series_lar=[]
  data_series_sma=[]
  page_names=[]

  # Plot Grasp2K numerical orbitals as we need them P(r)/r and Q(r)/r
  if mode=='plot_grasp_orbitals':

    # Create the .vsz file.
    with open(atom.element+'_grasp_numerical_orbitals.vsz','w') as g:

      # Loop over the number of orbitals
      for i in range(atom.orbital_info.num_orbitals):

        # There exists a problem in veusz with the page names that have _- in them.
        #   This fix avoids that error by substituting an 'n' for the negative sign.
        orb_text=atom.element+"_"+str(atom.orbital_info.n_qnums[i])+"_"+str(atom.orbital_info.k_qnums[i]).replace('-', 'n')

        # Append to the list of page names.
        page_names.append(orb_text)

        # Append the names for each radial, large, and small data set to the
        #  respective lists.
        data_series_rad.append(orb_text+'_radialcomp')
        data_series_lar.append(orb_text+'_largecomp')
        data_series_sma.append(orb_text+'_smallcomp')

        # Make the legend names for each orbital in this case we're only getting the j quantum number for the large component.
        legend_labels_large.append(str(atom.orbital_info.n_qnums[i])+str(org_radwavefn.ltolsym(org_radwavefn.lktol(atom.orbital_info.k_qnums[i])))+"("+str(j_qnum(atom.orbital_info.k_qnums[i]))+")")
        legend_labels_small.append(str(atom.orbital_info.n_qnums[i])+str(org_radwavefn.ltolsym(org_radwavefn.lktol(atom.orbital_info.k_qnums[i])))+'_small')

        # Start writing out the raw data sets. They go at the beginning of the file.
        #   Write out the radial component name and then all of the radial values.
        g.write("ImportString(u'"+orb_text+"_radialcomp(numeric)','''\n")
        for rad_comp in atom.orbital_info.radial_comps[i]:
          g.write(str(rad_comp)+'\n')
        g.write("''')\n")

        #   Write out the large component name and then all of the large values.
        g.write("ImportString(u'"+orb_text+"_largecomp(numeric)','''\n")
        for l_r in atom.orbital_info.large_div_r[i]:
          g.write(str(l_r)+'\n')
        g.write("''')\n")

        #   Write out the small component name and then all of the small values.
        g.write("ImportString(u'"+orb_text+"_smallcomp(numeric)','''\n")
        for s_r in atom.orbital_info.small_div_r[i]:
          g.write(str(s_r)+'\n')
        g.write("''')\n")

      # Write out the style for the vsz file.
      g.write("Set('colorTheme', u'default-latest')\n")
      g.write("Set('StyleSheet/axis-function/autoRange', u'next-tick')\n")

      # Zip together all of the names associated with each orbital as defined above
      #   and loop over them.
      for page,large_legend,small_legend,data_rad,data_lar,data_sma in zip(page_names,legend_labels_large,legend_labels_small,data_series_rad,data_series_lar,data_series_sma):

        # Add a page.
        g.write("Add('page', name=u'"+page.replace('n','-')+"', autoadd=False)\n")
        g.write("To(u'"+page.replace('n','-')+"')\n")

        # Add in a graph.
        g.write("Add('graph', name=u'"+page.replace('n','-')+"', autoadd=False)\n")
        g.write("To(u'"+page.replace('n','-')+"')\n")

        # Add in an x axis.
        g.write("Add('axis', name=u'x', autoadd=False)\n")
        g.write("To(u'x')\n")
        g.write("Set('label', u'r (angstroms)')\n")
        g.write("Set('min', u'Auto')\n")
        g.write("Set('max', u'Auto')\n")
        g.write("Set('Label/font', u'LMRoman10')\n")
        g.write("Set('TickLabels/font', u'LMRoman10')\n")
        g.write("To('..')\n")

        # Add in a y axis.
        g.write("Add('axis', name=u'y', autoadd=False)\n")
        g.write("To(u'y')\n")
        g.write("Set('label', u'\\chi(r)/r')\n")
        g.write("Set('min', u'Auto')\n")
        g.write("Set('max', u'Auto')\n")
        g.write("Set('direction', u'vertical')\n")
        g.write("Set('Label/font', u'LMRoman10')\n")
        g.write("Set('TickLabels/font', u'LMRoman10')\n")
        g.write("To('..')\n")

        # Add in the large component data set as an xy type.
        g.write("Add('xy', name=u'"+data_lar+"', autoadd=False)\n")
        g.write("To(u'"+data_lar+"')\n")
        g.write("Set('markerSize', u'1.5pt')\n")
        g.write("Set('xData', u'"+data_rad+"')\n")
        g.write("Set('yData', u'"+data_lar+"')\n")
        g.write("Set('key', u'"+large_legend+"')\n")
        g.write("To('..')\n")

        # Add in the small component data set as an xy type.
        g.write("Add('xy', name=u'"+data_sma+"', autoadd=False)\n")
        g.write("To(u'"+data_sma+"')\n")
        g.write("Set('markerSize', u'1.5pt')\n")
        g.write("Set('xData', u'"+data_rad+"')\n")
        g.write("Set('yData', u'"+data_sma+"')\n")
        g.write("Set('key', u'"+small_legend+"')\n")
        g.write("Set('PlotLine/color', u'auto')\n")
        g.write("Set('MarkerLine/color', u'auto')\n")
        g.write("Set('MarkerFill/color', u'auto')\n")
        g.write("To('..')\n")

        # Add in the y=0 line for the graph.
        g.write("Add('function', name=u'yax', autoadd=False)\n")
        g.write("To(u'yax')\n")
        g.write("Set('function', u'0')\n")
        g.write("Set('variable', u'y')\n")
        g.write("Set('Line/color', u'black')\n")
        g.write("To('..')\n")

        # Add in the x=0 line for the graph.
        g.write("Add('function', name=u'xax', autoadd=False)\n")
        g.write("To(u'xax')\n")
        g.write("Set('function', u'0')\n")
        g.write("Set('Line/color', u'black')\n")
        g.write("To('..')\n")

        # Add in the legend.
        g.write("Add('key', name=u'key1', autoadd=False)\n")
        g.write("To(u'key1')\n")
        g.write("Set('title', u'')\n")
        g.write("Set('horzPosn', u'right')\n")
        g.write("Set('vertPosn', u'top')\n")
        g.write("Set('horzManual', 0.0)\n")
        g.write("Set('vertManual', 0.0)\n")

        # Go back to the root to add in the next page.
        g.write("To('..')\n")
        g.write("To('..')\n")
        g.write("To('..')\n")
 
  # If there's an equations.dat file which has the form:
  #   orb1 A*(x**l)*e**(-alpha*x**2)+...
  #
  #   orb2 A*(x**l)*e**(-alpha*x**2)+...
  elif mode=='plot_orbital_equations':

    # Create holding lists.
    orbital=[]
    equation=[]

    # Make a list of files in the directory. If there isn't
    #  an equations.dat file then don't run.
    files=os.listdir()

    # If there is an equations.dat file then create the .vsz
    if 'equations.dat' in files:

      # Read in the data for the equations.dat file.
      with open('equations.dat') as g:
        for line in g:
          if '\n' in line:
            t=line.strip()
          else:
            t=line
          l=t.split()
          
          if len(l)==0:
            continue
 
          orbital.append(l[0])
          equation.append(l[1])
 

      # Open the vsz file for the orbital equations graphs.
      with open(atom.element+'_orbital_equations_graph.vsz','w') as g:

        # Loop over the orbital, equation pairs.
        for orb,eqn in zip(orbital,equation):

          # Create the page associated with the orbital.
          g.write("Add('page', name=u'"+orb+"', autoadd=False)\n")
          g.write("To(u'"+orb+"')\n")

          # Create the graph.
          g.write("Add('graph', name=u'"+orb+"', autoadd=False)\n")
          g.write("To(u'"+orb+"')\n")

          # Add an x axis.
          g.write("Add('axis', name=u'x', autoadd=False)\n")
          g.write("To(u'x')\n")
          g.write("Set('label', u'r (angstroms)')\n")
          g.write("Set('min', u'Auto')\n")
          g.write("Set('max', u'Auto')\n")
          g.write("Set('Label/font', u'LMRoman10')\n")
          g.write("Set('TickLabels/font', u'LMRoman10')\n")
          g.write("To('..')\n")

          # Add a y axis.
          g.write("Add('axis', name=u'y', autoadd=False)\n")
          g.write("To(u'y')\n")
          g.write("Set('label', u'\\chi(r)/r')\n")
          g.write("Set('min', u'Auto')\n")
          g.write("Set('max', u'Auto')\n")
          g.write("Set('direction', u'vertical')\n")
          g.write("Set('Label/font', u'LMRoman10')\n")
          g.write("Set('TickLabels/font', u'LMRoman10')\n")
          g.write("To('..')\n")

          # Add in the orbital description as a function type.
          g.write("Add('function', name=u'"+orb+"', autoadd=False)\n")
          g.write("To(u'"+orb+"')\n")
          g.write("Set('function', u'"+eqn+"')\n")
          g.write("Set('key', u'"+orb+"')\n")
          g.write("Set('steps',10000)\n") 
          g.write("Set('Line/color', u'black')\n")
          g.write("To('..')\n")
 
          # Add in the y=0 line for the graph.
          g.write("Add('function', name=u'yax', autoadd=False)\n")
          g.write("To(u'yax')\n")
          g.write("Set('function', u'0')\n")
          g.write("Set('variable', u'y')\n")
          g.write("Set('Line/color', u'black')\n")
          g.write("To('..')\n")
 
          # Add in the y=0 line for the graph. 
          g.write("Add('function', name=u'xax', autoadd=False)\n")
          g.write("To(u'xax')\n")
          g.write("Set('function', u'0')\n")
          g.write("Set('Line/color', u'black')\n")
          g.write("To('..')\n")
          g.write("Add('key', name=u'key1', autoadd=False)\n")
          g.write("To(u'key1')\n")
          g.write("Set('title', u'')\n")
          g.write("Set('horzPosn', u'right')\n")
          g.write("Set('vertPosn', u'top')\n")
          g.write("Set('horzManual', 0.0)\n")
          g.write("Set('vertManual', 0.0)\n")
          g.write("To('..')\n") 
          g.write("To('..')\n")
          g.write("To('..')\n") 
    else:
      print('No equations file present.\n')

if __name__=='__main__':
  stime=time.time()

  print('Running Module File')

  # Read in the data from the isodata and rwfn files.
  #  This consitutes pretty much all of the needed
  #    information to fit and graph the orbitals.
  element=org_radwavefn.atomic_system('isodata','rwfn.out')

  create_veusz_file('plot_grasp_orbitals', element)
  create_veusz_file('plot_orbital_equations', element)

  print('Time to run: '+str(time.time()-stime)+'s.')
