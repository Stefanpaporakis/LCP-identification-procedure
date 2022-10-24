# LCP-identification-procedure
Code for identifying LCPs from 1D SAXS patterns

A few tips to using this code

datapath:
This is where you have your data stored. Make sure that the file format is ".dat" or ".npy"
Make sure you have the "" and the // at the end

outpath:
This is where the code will save the outcome of your data to.
The codes output is in the Data_summary.csv file, it should include Dataset file title,
the lattice parameter and found phases.
typically you'll keep this the same as the datapath

q limits:
"q" is the scattering vector, and basically your x-axis value for the plots.
Essentially the qmin and qmax values determine where the code looks for peaks 
and plots them on the profile.

Prominence:
This is the peak finding algorithm sensitivity. If you want the code to find
more peaks (the red lines on the graphs), decrease the number, if you want less peaks, increase the number.
This value should help you filter out background noise from actual peaks.
Have a play around see what works for your data.

Match tolerance:
Peaks can sometime be close together or not perfect fits for the model LCPs.
The modelled LCPs (lam, hex, pn3m etc) are all based of a series of lattice 
ratios in the main script of the code. The match tolerance is how close a found peak
is to a model peak. Lets say the match tolerance is set to 0.005, if the code 
finds a peak at q=0.110, and there's a LCP model peak at q=0.108, then it's
within the match tolerance and the code will identify that peak as a match. 
If the values are outside the match_tolerance, then it won't say it's a match. 

Plot parameters:
These are a bunch of different quality of life parameters. You can turn log scales off and on for the axes by setting
Logx or Logy to False or True. Same as the title of the plot, you can remove it by setting plot_title = False. Code is 
set to call each plot by the name of the file, this can be changed to a custom file name by setting File_name_as_title = False

Realistically, lines 19-50 are the only things that should be changed upon using the identification procedure.
