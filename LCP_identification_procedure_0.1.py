# -*- coding: utf-8 -*-
"""
authors: Stefan Paporakis, Jack Binns, Andrew Martin
"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import seaborn as sns
import math
sns.set()

########################################################################################################################
########################################################################################################################


#datapath
datapath = "path to your files here"
outpath = datapath


# q limits
qmin = 1
qmax = 2

# Prominence
prominence =10

# Match tolerance
match_tolerance=0.005

# Plot parameters
# Plot title on/off?
plot_title = True
# File name as title of the plot?
File_name_as_title = True
# Custom plot title
custom_title = "insert custom plot title here"
# Title font size
title_font_size = 16
# Log axis?
Logx = False
Logy = False
# Axis titles
Xaxis = "insert x axis title here"
Yaxis = "insert y axis title here"
font_size = 12

########################################################################################################################
########################################################################################################################


#Below are the functions used in the code

##############################################################################
####### Find the peaks, trim the q range, list the peaks
############### OBJECTS: PEAKHUNT, QTRIM, LOAD_AND_PEAK_FIND, WRITE PEAK LIST TO OUTPATH
# Find peaks.
  

def load_and_peakfind(file_path,qmin,qmax,prominence): 
    profile = np.loadtxt(file_path,skiprows=2)
    trim_profile = QTrim(profile,qmin,qmax)
    print("Finding peaks...")
    peaks_in_q = PeakHunt(trim_profile, prominence)
    print("number of peaks found is:", len(peaks_in_q))
    print("Peaks found at :", peaks_in_q)
    return trim_profile, peaks_in_q

def PeakHunt(profile, prom):
    peak_find, _ = find_peaks(profile[...,1], prominence = prom)
    peaks_in_q = []
    for i, datap in enumerate(profile):
        if i in peak_find:
            peaks_in_q.append(profile[i][0])
    peaks_in_q = np.array(peaks_in_q)
    return peaks_in_q

# Trim scattering vector (q)
def QTrim(profile,qmin,qmax):
    #trims to qmin qmax
    qdata = []
    for point in profile:
        if qmin <= point[0] <= qmax:
            qdata.append(point)
    qdata = np.array(qdata)
    return qdata

# Pick up data and fit peak finding program    


 # Write peak list to outpath    
def write_peaks(peak_list, outpath):
    out_file = open(outpath, "w+") 
    for i in peak_list:
        out_string = ""
        out_string +=str(i)
        out_string += "\n"
        out_file.write(out_string)
    out_file.close()

##############################################################################
###### Generate LCPs, match the peaks, calculate the LP
############## OBJECTS: GENMODELPEAKS, PEAKMATCH, LATTICSPAR    
# Generate model peaks based off lattice ratios
def GenModelPeaks(mdl_flag, root,match_tolerence):
    mdl_peaks = []
    
    if mdl_flag == 'pn3m':
        a_pn3m = (pn3m[0] / root)
        for hkl in pn3m:
            model_peak = hkl / a_pn3m
            if qmin < model_peak < qmax:
                mdl_peaks.append(model_peak)
        return mdl_peaks

    if mdl_flag == 'hex':
        a_hex = ((2*hex[0])/(root*np.sqrt(3)))
        for hkl in hex:
            model_peak = ((2 / np.sqrt(3)) * hkl) / a_hex
            if qmin < model_peak < qmax:
                mdl_peaks.append(model_peak)
        return mdl_peaks
    
    
    if mdl_flag == 'lam':
        a_lam= (lam[0] / root)
        for hkl in lam:
            model_peak = hkl / a_lam
            if qmin < model_peak < qmax:
                mdl_peaks.append(model_peak)     
        return mdl_peaks
    
    
    if mdl_flag == 'ia3d':
        a_ia3d= (ia3d[0] / root)
        for hkl in ia3d:
            model_peak = hkl/ a_ia3d
            if qmin < model_peak < qmax:
                mdl_peaks.append(model_peak)
        return mdl_peaks
    
    if mdl_flag == 'im3m':
        a_im3m= (im3m[0] / root)
        for hkl in im3m:
            model_peak = hkl /a_im3m
            if qmin < model_peak < qmax:
                mdl_peaks.append(model_peak)
        return mdl_peaks
    
    if mdl_flag == 'ignore':
        a_ignore=(ignore[0]/root)
        for hkl in ignore:
            model_peak=hkl/a_ignore
            if qmin < model_peak < qmax:
                mdl_peaks.append(model_peak)
        return mdl_peaks
    
# Match found peaks to model peaks    
def PeakMatching(model_peaks, exp_peaks, mdl):
    tolerance = 0.002
    n_matches = 0
    for epeak in exp_peaks:
        for mpeak in model_peaks:
            diff = abs(epeak - mpeak)
            if diff <= tolerance:
                n_matches = n_matches + 1
    return n_matches 

# Calculate lattice parameters
def LatticePar(plot_mdl,  root):

    if plot_mdl == 'pn3m':
        a_pn3m = (2*math.pi*pn3m[0] / root)
        print ("pn3m lattice parameter is", a_pn3m)
        lp_out = datapath + file_tag + "_pn3m_lp.txt"
        y= np.array([a_pn3m])
        LatticeP.append(a_pn3m)
        
    if plot_mdl == 'hex':
        a_hex = ((4*math.pi*hex[0])/(root*np.sqrt(3)))
        print ("hex lattice parameter is", a_hex)
        lp_out = datapath + file_tag + "_hex_lp.txt"
        y= np.array([a_hex])
        LatticeP.append(a_hex)
        
    if plot_mdl == 'lam':
        a_lam= (2*math.pi*lam[0] / root)
        print("lam lattice parameter is", a_lam)
        lp_out = datapath + file_tag + "_lam_lp.txt"
        y= np.array([a_lam])
        LatticeP.append(a_lam)
        
    if plot_mdl == 'ia3d':
        a_ia3d= (2*math.pi*ia3d[0] / root)
        print("ia3d lattice parameter is", a_ia3d)
        lp_out = datapath + file_tag + "_ia3d_lp.txt"
        y=np.array([a_ia3d])
        LatticeP.append(a_ia3d)

    if plot_mdl == 'im3m':
        a_im3m= (2*math.pi*im3m[0] / root)
        print("im3m lattice parameter is", a_im3m)
        lp_out = datapath + file_tag + "_im3m_lp.txt"
        y= np.array([a_im3m])
        LatticeP.append(a_im3m)
        
    if plot_mdl == 'ignore':
        LatticeP.append("--")
        
        
##############################################################################        
######CLASS: Plot found to matched peaks, choose compare model
############# OBJECTS: INITIALPLOT, PLOTFITS, MATCHPEAKS, SINGLEFITINSPECTION


# Plot profile with all theoretical LCP models before analysis
def InitialPLot(labels, root, prf, piq):
        phases_to_test = [lam, hex, ia3d, pn3m, im3m, ignore]
        refpeaks = [0, 0, 0, 0, 0, 0]
        legend = ["Found peaks", "Lam", "hex",  "ia3d", "Pn3m", "im3m", "ignore"]
        colours = ["pale red", "purple", "green", "orange", "blue", "black", "brown"]
        fig, ax = plt.subplots(figsize=(7, 5))
        plt.plot(prf[..., 0], prf[..., 1])
        for i in np.arange(len(piq)):
            plt.axvline(x=piq[i], ymin=0, ymax=0.2, color=sns.xkcd_rgb[colours[0]])
            plt.text(0.05, 0.1, legend[0], transform=ax.transAxes, color=sns.xkcd_rgb[colours[0]], fontsize=12, verticalalignment='top')
        for j, p in enumerate(phases_to_test):
            for peak in p:
               plt.axvline(x=piq[refpeaks[j]]*peak/p[0], ymin=0.1+(j+1)*0.1, ymax=0.2+(j+1)*0.1, color=sns.xkcd_rgb[colours[j+1]])   
               plt.text(0.05, 0.15+ 0.1*(j+1), legend[j+1], transform=ax.transAxes, color=colours[j+1], fontsize=12, verticalalignment='top')
        plt.xlim(qmin,qmax)
        if Logx == True:
            plt.xscale("log")
        if Logy ==True:
            plt.yscale("log")
        plt.xlabel(Xaxis, fontsize = font_size)
        plt.ylabel(Yaxis, fontsize = font_size)
        plt.grid(False)
        if plot_title is True and File_name_as_title is True:
            plt.title(split_path[-1],fontsize = title_font_size)
        if plot_title is True and File_name_as_title is False:
            plt.title(custom_title,fontsize = title_font_size)
        plt.show()
        
 # Plot chosen model over SAXS profile   
def PlotFits(prf,piq,plot_mdl_peaks,plot_mdl):
    fig = plt.figure(figsize=(7, 5))
    plt.plot(prf[..., 0], prf[..., 1])
    for i in np.arange(len(piq)):
        plt.axvline(x=piq[i], ymin=0, ymax=0.2, color=sns.xkcd_rgb["pale red"])
    for j in np.arange(len(plot_mdl_peaks)):
            plt.axvline(x=plot_mdl_peaks[j], ymin=0.2, ymax=0.4, color=sns.xkcd_rgb["purple"])
    plt.xlabel(Xaxis, fontsize = font_size)
    plt.ylabel(Yaxis, fontsize = font_size)
    if Logx == True:
        plt.xscale("log")
    if Logy ==True:
        plt.yscale("log")
    plt.xlim(qmin,qmax)
    plt.grid(False)
    if plot_title is True and File_name_as_title is True:
        plt.title(split_path[-1],fontsize = title_font_size)
    if plot_title is True and File_name_as_title is False:
        plt.title(custom_title,fontsize = title_font_size)
    plt.savefig( outpath+file_tag+"_"+plot_mdl+'_plot.png')
    plt.show()
    
 # Match found peaks to model peaks within tolerence value   
def Matchpeaks(match_tolerance,piq,plot_mdl_peaks,labels,root,prf):
    unmatched_peaks = []
    for e_peak in np.array(piq):
        match_flag=False
        for m_peak in np.array(plot_mdl_peaks):
            diff = abs(e_peak - m_peak)
            if diff < match_tolerance:
                match_flag=True
                break 
            else:
                match_flag = False
        if match_flag is False:
            unmatched_peaks.append(e_peak)
        elif match_flag is True:
            continue
    return unmatched_peaks




# Choose model to compare with found peaks
def SingleFitInspection(labels,root,prf,piq):
    InitialPLot(labels, root, prf, piq)    
    flag = input(" inspect? [y/n]")
    if flag == "y":
        inspection_flag = True
        while inspection_flag is True:
            print(labels)
            lbl = int(input("Plot which model ? [0:5]?"))
            plot_mdl = labels[lbl]
            call_LP = LatticePar(plot_mdl,  root)  
            plot_mdl_peaks = GenModelPeaks(plot_mdl, root, match_tolerance)
            PlotFits(prf, piq, plot_mdl_peaks, plot_mdl)
            print (plot_mdl+" model peaks are",plot_mdl_peaks)
            Phases.append(plot_mdl)
            with open(outpath+'Data_summary.csv', 'a+') as summary:
                Current_file = File_title[-1]
                Current_phase  = Phases[-1]
                Current_latticeP = LatticeP[-1]
                string = (f'{Current_file}, {Current_phase}, {Current_latticeP}')
                summary.write(string+'\n')
               
            unmatch=Matchpeaks(match_tolerance, piq, plot_mdl_peaks,labels,root,prf)
            print("unmatched_peaks:")
            print(unmatch)
            return unmatch, inspection_flag
    elif flag == "n":
        inspection_flag = False
        unmatch= []
        return unmatch, inspection_flag

###############################################################################
######## MAIN CODE BELOW
###############################################################################
if __name__ == "__main__":
   
    File_title = ['file name']
    Phases = ['identified phase']
    LatticeP = ['lattice parameter']
     
    # Lattice ratios/d-spacings:
    
    ignore= np.array([np.sqrt(1),np.sqrt(1.000001)])
    pn3m = np.array([np.sqrt(2), np.sqrt(3), np.sqrt(4), np.sqrt(6),
                     np.sqrt(8),np.sqrt(9), np.sqrt(10),np.sqrt(11),np.sqrt(12)] )
    hex = np.array( [1, np.sqrt(3), np.sqrt(4)])
    im3m = np.array( [np.sqrt(2),np.sqrt(4),np.sqrt(6),np.sqrt(8),np.sqrt(10),np.sqrt(12),np.sqrt(14),np.sqrt(16),np.sqrt(18),np.sqrt(20),np.sqrt(22),
                      np.sqrt(22),np.sqrt(24),np.sqrt(26)] )
    ia3d = np.array([np.sqrt(6), np.sqrt(8), np.sqrt(14), np.sqrt(16), np.sqrt(20), np.sqrt(22),np.sqrt(24),np.sqrt(26)] )
    lam = np.array([1,2,3,4])
   
# Main script
    print("\n###############################################")
    print("multiphase analysis")
    print("###############################################")
    print("")
    print("")
    
    files = sorted(glob.glob(datapath+"*.dat"))
    print("Dataset : ",)
    print("Total number of files in data set : ",len(files))
    startup = input("start at which file? [0:"+str(len(files)-1)+"]")
    counter = int(startup)
    for data in files[int(startup):]:
        split_path = data.split("\\") 
        file_name = split_path[-1]
        File_title.append(file_name)
        file_tag = file_name[:-4] 
        counter = counter+1
        print("\n##################################################################################")
        print("looking at file: "+split_path[-1]+". File "+str(counter)+ " of "+str(len(files)))
        print("##################################################################################")
       
        #found peaks by algorithm
        profile, peaks_in_q = load_and_peakfind(data,qmin,qmax,prominence)
        if len(peaks_in_q)==0:
            print("no peaks found in dataset, here is the profile")
            fig, ax=plt.subplots()
            ax.set_facecolor("white")
            ax.plot(profile[..., 0], profile[..., 1])
            if Logx == True:
                plt.xscale("log")
            if Logy ==True:
                plt.yscale("log")
            plt.xlim(qmin,qmax)
            plt.xlabel(Xaxis, fontsize = font_size)
            plt.ylabel(Yaxis, fontsize = font_size)
            plt.grid(False)
            if plot_title is True and File_name_as_title is True:
                plt.title(split_path[-1],fontsize = title_font_size)
            if plot_title is True and File_name_as_title is False:
                plt.title(custom_title,fontsize = title_font_size)
            plt.show()
            
            have_a_look = input("continue? [y/n]")
            if have_a_look =="y":
                continue
            elif have_a_look == "n":
                break
        np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
        

         # Define first peak in data
        peak_zero = peaks_in_q[0]
        print("Model is generating peaks using peak: ", peak_zero)
        
        # Generate models and display how many model peaks match found peaks
        model_labels = ['pn3m', 'hex', 'lam', 'im3m', 'ia3d','ignore' ] 
        results_dict = {}
        for mdl in model_labels:
            model_peaks = GenModelPeaks(mdl, peak_zero,match_tolerance)
            matches = PeakMatching(model_peaks, peaks_in_q, mdl)
            fit = matches / len(model_peaks)
            results_dict[str(mdl)] = fit
           
            if fit >0.6:
                print( mdl, "found ", str(matches), " matching peaks out of ", str(len(model_peaks)),
                  " model peaks with a single phase fit of",100.*fit,"%")
            else: 
                print("no match for "+mdl+" single phase")
        
            
                
        # Apply new model on stripped profile if peaks are still found by algorithm,
        # repeats until all peaks are matched with LCP model.
        
       # open(datapath+file_tag+"_found_phases.txt", 'w').close()
        unmatch, insflag =SingleFitInspection(model_labels, peak_zero, profile, peaks_in_q)
        if insflag is True and len(unmatch)>0:
            newfit = input("fit new phase? [y/n]")
            if newfit =="y":
                newfitflag = True
            elif newfit == "n":
                newfitflag = False
            while newfitflag is True and len(unmatch)>0:
                question = input("Model new phase from which unmatched peak?[0:"+str(len(unmatch)-1)+"]")
                peak_zero=unmatch[int(question)]
                print('new peak zero is', peak_zero)
                unmatch, insflag= SingleFitInspection(model_labels, peak_zero, profile, unmatch)
            while newfitflag is True and len(unmatch)==0:
               # phases=np.loadtxt(datapath+file_tag+"_found_phases.txt",dtype=str,delimiter='\n')
               # phasearray=str(phases)
                cont=input("no peaks left, "
                          # ", matched phases for "+split_path[-1]+" were "+phasearray+".
                           "Move on to next dataset? [y,n]")
                if cont=="y":  
                    break
                elif cont=="n":
                    contflag= False
                while contflag is False:
                    print ("code has stopped running at "+split_path[-1]+", file "+str(counter) +" of "+str(len(files)))
                    sys.exit()
            while newfitflag is False and len(unmatch)>0:
                #phases=np.loadtxt(datapath+file_tag+"_found_phases.txt",dtype=str,delimiter='\n')
               # phasearray=str(phases)
                print( split_path[-1] + ' has '+ str(len(unmatch))+' unmatched peaks remaining!' )
              #  print ('matched phases for '+split_path[-1]+' are '+phasearray+'. with ' +str(len(unmatch))+' unmatched peaks remaining!')
                newcont=input("continue?[y,n]")
                if newcont=="y":
                    break
                elif newcont=="n":
                    print ("code has stopped running at file name: "+split_path[-1]+", file "+str(counter) +" of "+str(len(files)))
                    sys.exit()  
        elif insflag is True and len(unmatch)==0:
            print('no more peaks to match')
            continue
        elif insflag is False:
            
            continue
            
       