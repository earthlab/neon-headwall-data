# Georegistering Headwall Radiance Cubes
to produce georegistered radiance data from the headwall system, need to load the flightlines into ENVI 

 

Example: 


Input files:  

1) raw_7296_rd (.hdr) 

2) raw_7296_rd_igm (.hdr) 

 

In ENVI... use the Georeference from IGM tool in the Geometric Correction group 


 

Select the input data file in the dialog. Here, raw_7296_rd. Leave all defaults. 


 

In the next dialog, select the X Geometry values by selecting band 1 of the IGM file. Similarly, select Band 2 for the Y geometry values. 


 

Select the projection for the output. The first example shows it as UTM (can change Datum to WGS–84 if desired) 


 

Whereas this example shows choosing Geographic Lat/Lon as the output (EPSG:4326). Either is fine, just be aware of what was chosen. 


 

Finally, BE SURE TO SET OUTPUT ROTATION TO ZERO!!!!! It is not 0 by default, but rather, some number to minimize storage. If the data is not north-up, then most software I have come across to open the file will do so incorrectly. 

 

Choose an output file name (no extension) 


 
