# Selfcal
Written by Joseph Smidt (I think), second author on Cooray et al 2012 and Zemcov et al 2014, based on the least squares self-calibration method of Arendt et al 2000. Archiving for future use. The code was originally written for Spitzer IRAC, so I am not entirely sure it can be extended to other types of data if it is doing IRAC specific things, but may be worth a try. I don't know what most of the files in this code do and there is almost no documentation, but I have been able to get it to mosaic data together, so maybe it will work for others too.

To use:

1. Download the selfcal package to the directory of your choice.

2. The file that runs everything is called 'escargot_main.pro'. Before running, open this file up in a text editor. You will need to change three things in the file before you can run it. The first is on line 65, where you will have to edit the variable 'tag_save'. This is just the directory name and/or file prefix where you would like the output of the code (i.e., the mosaic) to be saved to. Next is on line 92, where you need to edit the variable 'files' to the name of the directory where your input data tiles are located. I'll discuss inputs in the next bullet point. Last, on lines 97 and 98, you need to state the dimensions of your input data tiles via 'XPIX' and 'YPIX'. The default is for IRAC images, which are 256 pixels by 256 pixels.

3. The code takes .fits files as input. These fits files must all be located in a directory specified by the 'files' variable in escargot_main.pro.

4. Fire up IDL in the directory containing the Selfcal package and run the command ".r escargot_main.pro" Hopefully it runs. Depending on how many tiles you are mosaicking together it can take up to 12 hours (which it does to mosaic the 11,000 Spitzer IRAC data tiles that I am working on) or 10-15 minutes if there aren't too many tiles to mosaic.

A few other things to note: At certain times while it is running, some plots and partial images may pop up. I don't know what they are for, but it is ok to close them without affecting the code.
Also, running it over ssh does not seem to work. I think it might be that the plots/images that pop up cause some sort of X11 forwarding error, but I haven't gone deep enough into the code to figure out how to turn these pop ups off.
