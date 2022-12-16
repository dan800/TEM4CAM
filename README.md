# TEM4CAM
Code to post-process CESM-CAM/WACCM transformed Eulerian mean output to useful circulation diagnostics (e.g., v*, w*).
Based on code from Isla Simpson.

The code follows the 'TEM recipe' from Appendix A of 
Gerber, E. P. and Manzini, E.: The Dynamics and Variability Model Intercomparison Project (DynVarMIP) for CMIP6: assessing the stratosphere–troposphere system, Geosci. Model Dev., 9, 3413–3425, https://doi.org/10.5194/gmd-9-3413-2016, 2016.

pdf available here: https://gmd.copernicus.org/articles/9/3413/2016/gmd-9-3413-2016.pdf

Output
```
Table A1. Momentum budget variable list (2-D monthly / daily zonal means, YZT).

Name      Long name [unit]
epfy      northward component of the Eliassen–Palm flux [m3 s−2]
epfz      upward component of the Eliassen–Palm flux [m3 s−2]
vtem      Transformed Eulerian mean northward wind [m s−1] wtem Transformed Eulerian mean upward wind [m s−1]
psitem    Transformed Eulerian mean mass stream function [kg s−1]
utendepfd tendency of eastward wind due to Eliassen–Palm flux divergence [m s−2]
utendvtem tendency of eastward wind due to TEM northward wind advection and the Coriolis term [m s−2] 
utendwtem tendency of eastward wind due to TEM upward wind advection [m s−2]
```
