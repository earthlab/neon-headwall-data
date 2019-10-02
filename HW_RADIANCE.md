to produce georegistered radiance data from the headwall system, need to load the flightlines into ENVI 

 

Example: 


Input files:  

1) raw_7296_rd (.hdr) 

2) raw_7296_rd_igm (.hdr) 

 

In ENVI... use the Georeference from IGM tool in the Geometric Correction group 


 
![GetImageAttachment (5)](https://user-images.githubusercontent.com/4762214/64807141-d094f100-d551-11e9-96a9-01f8c58ab885.png)


Select the input data file in the dialog. Here, raw_7296_rd. Leave all defaults. 

![GetImageAttachment (4)](https://user-images.githubusercontent.com/4762214/64807123-c83cb600-d551-11e9-859b-2defb5a0b65c.png)
 

In the next dialog, select the X Geometry values by selecting band 1 of the IGM file. Similarly, select Band 2 for the Y geometry values. 



 
![GetImageAttachment (3)](https://user-images.githubusercontent.com/4762214/64807101-bb1fc700-d551-11e9-86c4-35e0f69614bc.png)

Select the projection for the output. The first example shows it as UTM (can change Datum to WGS –84 if desired) 



![GetImageAttachment (2)](https://user-images.githubusercontent.com/4762214/64807061-a7746080-d551-11e9-9b39-50d30b93dd27.png)

Whereas this example shows choosing Geographic Lat/Lon as the output (EPSG:4326). Either is fine, just be aware of what was chosen. 


 ![GetImageAttachment (1)](https://user-images.githubusercontent.com/4762214/64807046-a0e5e900-d551-11e9-986b-970d79bab6ef.png)

 

Finally, BE SURE TO SET OUTPUT ROTATION TO ZERO!!!!! It is not 0 by default, but rather, some number to minimize storage. If the data is not north-up, then most software I have come across to open the file will do so incorrectly. 

 
![GetImageAttachment](https://user-images.githubusercontent.com/4762214/64806999-87dd3800-d551-11e9-9195-efdaed77ee7b.png)


Choose an output file name (no extension) 


 
