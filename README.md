# Health-facility-catchment-area-generation
Generate health facility catchment areas from accessibility mapping in R and r.cost aglorithm in QGIS. Heath facility catchment area is the area from which a health facility attracts patients.
This script requires the two user supplied datasets:
 (a) The friction surface, which is available here:  http://www.map.ox.ac.uk/accessibility_to_cities/
 (b) A user-supplied .csv of points (i.e., known geographic coordinates) 

Notes:
 (a) All file paths and names should be changed as needed.
 (b) Important runtime details can be found in the comments.
 (c) This script is suitable only for analyses of moderately sized areas for most (e.g., up to 10 million km^2 in lower latitude settings - GLOBAL RUNS WILL NOT WORK).
     We recommend using Google Earth Engine for larger areas, with the exception of high-latitude areas where custom approaches are typically required.

 Citation: D.J. Weiss, A. Nelson, H.S. Gibson, W. Temperley, S. Peedell, A. Lieber, M. Hancher, E. Poyart, S. Belchior, N. Fullman, B. Mappin, U. Dalrymple, J. Rozier, 
 T.C.D. Lucas, R.E. Howes, L.S. Tusting, S.Y. Kang, E. Cameron, D. Bisanzio, K.E. Battle, S. Bhatt, and P.W. Gething. A global map of travel time to cities to assess 
 inequalities in accessibility in 2015. (2018). Nature. doi:10.1038/nature25181.
