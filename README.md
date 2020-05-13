# Point-Set-Registration-of-CIM-on-Fengyun-3-MWRI-Data


>`The code of proposed geolocation error estimation method is written in MATLAB Version R2018b (Math Works, Inc.) and the process of error correction is running on the CentOS Linux operating system on the VMware virtual machine. All codes are provided publicly` [online](https://github.com/Jiazheng-Liu/Point-Set-Registration-of-CIM)
>
>`Detailed information of the experimental data set can be downloaded from the` [Fengyun Satellite Data Center](http://data.nsmc.org.cn/portalsite/default.aspx) 



### bulePoints.m

> Extract the coordinate information of the GSHHS corresponding to the detected coastline from the  .hdf into .mat 
> 
> This is the second paragraph in the blockquote.
>


### red_points.m
> Extract the coordinate information of the detected coastline from the .hdf into .mat
>

### haianxian.m
>> #### LLT_init.m
>> #### MR.m
>> #### norm_ind.m

>Match the red point set (detected coastline) and blue point set (Coastline truth point) obtained previously
>
>save the matching data in the target domain
>

### paintRegion.m
>>#### Plot_Color.m
>
> The main function of drawing the local brightness temperature map
>


### myHDF.m
>>#### hdf_create.m
>
> The calculated new matching relationship is written to the HDF file. The results will be corrected by satellite in orbit simulation
>


### pixelError.m
> Extract the coordinate information of the detected coastline from the .hdf into .mat
>

### gshhs_land_f.mat
> Global Self-consistent, Hierarchical, High-resolution Shoreline (GSHHS) data version 2.3.7 (June 15, 2017).
>
