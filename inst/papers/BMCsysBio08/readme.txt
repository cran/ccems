The purpose of this directory is to provide an actively maintained version of the R script  
Additional File in Radivoyevitch BMC Systems Biology 2008.

The origSupFiles directory holds the original files of this paper. It includes the data RNR.RData, 
the function definitions in fRt.r, and the script that runs them Rt.r. These codes
may continue to run on new versions of R by chance, but they will not be maintained
and are thus included only for reference purposes. 

The files RNR.RData and fRt.r were absorbed into ccems. The idea then is to 
maintain RtCcems.r and autoRt.r as scripts that functionally/conceptually, replaces Rt.r 
but only depend on ccems. This leaves out all E-shaped models as these are not supported. 
It also leaves out n-shaped models deemed too unlikely to be plausible 
and thus not automatically created by ems. 
 
