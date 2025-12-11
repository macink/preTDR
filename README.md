To reproduce the plots in the preTDR:

1) All plot macros are located in the 'plot_macros' folder
2) Include header files in the 'header_files' folder
3) To generate a plot, you must give the necessary filename corresponging to the plot that is being produced. The files used are in the 'root_files' folder.
   
For example, to generate the |t| distribution including all methods (i.e. method E, method L, and projection method) run: root -b 'preTDR_tDistribution_all_methods.C("filename.root")'

Where "filename.root" is the corresponding root file that the plot is being generated from (i.e. for this analysis, the filename = Coherent_phi_sartre_25_10_2_eAu_10x100.root).
	
	
