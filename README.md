# main function:
**start_icorre.m**


# parameters:

**Batch Directory (optional)** : the root directory for imaging data

**Batch Processing**: if False, only process one single session given by Batch Directory

**Save Settings As**: path for saving the user_settings.mat
	
**Scan Image Version**: should be 5
	
**Number of Frames to Average for Reverence Image**: self-expalnatory
	
**Seed Correction**: Initial rigid-correction for following rigid/non-rigid corrections
	
**Rigid Correction**: Max. Number of Repeats: if larger than 1, run motion correction repeatedly until criterion is satisfied
	
**Non-Rigid Correction: Max. Number of Repeats**: if larger than 1, run motion correction repeatedly until criterion is satisfied
	
**Rigid Correction: Max. Shift(pixels)**: maximum possible frame shifts for rigid correction
	
**Non-Rigid Correction Max. Deviation from Rigid Shift(pixels)**: maximum possibel shifts for non-rigid correction
	
**Non-Rigid Corection Patch Width**: patch width for non-rigid correction
	
**Error Tolerance(pixels)**: criterion for iterative motion corrections
	
**Reference Channel (0 for 1-color imaging)**: two-channel imaging (not tested)
	
**Follower Channel (0 for 1-color imaging)**: two-channel imaging (not tested)
	
**Save Concatenated Copy of Corrected Stacks? (T/F)**: set to T if a session is recorded in separate files

**Downsample Factor for Concatenated Stack**: downsample the session if it is too large for computer memory
	
**Delete MAT files at Completion? (T/F)**

# note
1. converting tiffs to mats helps align the imaging data to behavior data\
2. if running both rigid and non-rigid correction, the non-rigid correction will be run on the basis of rigid motion correction.

# File struture
## before
```
root_dir
└───session 1
│   │   file1.tif
│   │   file2.tif
│   │	...
└───session 2
    │   file1.tif
    │   file2.tif
    |   ...
└───session 3  
    │   file1.tif
    │   file2.tif
    |   ...
...
```
## after
```
root_dir
└───session 1
|   |   reg_info,mat
|   |   stack_info.mat
|   |   file_0DS4.tif  (downsampled)
|   └───raw
│       │   file1.tif
│       │   file2.tif
│       │   ...
|   └───mat
|       |    file_00001.mat
|       |    file_00002.mat
|       |    ...
|   └───registered (motion corrected)
|       |    _file_00001.tif
|       |    _file_00002.tif
|       |    ...
|   └───QualityControl (correlation figures and data)
|       |    metrics.png
|       |    QC_NRM.mat
└───session 2
    |   ...
└───session 3  
    |   ...
```
