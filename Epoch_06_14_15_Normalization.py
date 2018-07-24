from pylab import * 
from numpy import *
from astropy.io import fits
from PIL import Image
from numpy import concatenate
from numpy import interp
from numpy import set_printoptions
from numpy import arange 
from numpy import nan
import matplotlib.pyplot as plt
import os
import matplotlib.ticker as mtick
from numpy import sqrt
from numpy import ones
from numpy import convolve
import matplotlib.patches as mpatches
from numpy import polyfit
from numpy import polyval
#------------------------------------------------------------------------------
#Functions
def badpix(list_name,xmin_indice,xmax_indice,replace_val):
    for i in range(xmin_indice,xmax_indice):
        list_name[i]=replace_val
def smooth(y, box_pts):
    box = ones(box_pts)/box_pts
    y_smooth = convolve(y, box, mode='same')
    return y_smooth
def find_index(your_list,your_value):#this finds the INDEX of the element of the list CLOSEST to what you put in.
    x_1=list(your_list)
    return(x_1.index(min(x_1, key=lambda x:abs(x-your_value))))
def closest_value(your_list,value):# this returns the VALUE of the element of the list CLOSEST to what you put in.
    return(your_list[find_index(your_list,value)])
def make_polyfit_point(xlist,ylist,start,end):#give this function some data and a range and it spits out the average x and y
    x_ave=closest_value(xlist,(start+end)/2)
    start_index=find_index(xlist,closest_value(xlist,start))
    end_index=find_index(xlist,closest_value(xlist,end))
    y_sum=0
    for i in range(start_index,end_index,1):
        y_sum+=ylist[i]
    y_ave=closest_value(ylist,y_sum/(end_index-start_index))       
    return [x_ave,y_ave]

# These two functions below are dangerous: they change values...!
def remove_zero_error(error):
    TME1=[]
    for i in range (0,len(error)):
        if error[i] == 0:
            TME1.append(1)
        else:
            TME1.append(error[i])  
    return TME1
def remove_small(mylist,threshold):
    for i in range(0,len(mylist)):
        if mylist[i]<=threshold:
            mylist[i]=0
    return mylist
def find_index_from_value(mylist, value):
    return find_index(mylist,closest_value(mylist,value))   
#------------------------------------------------------------------------------
#Fits files Opening
#File says that the date obs was 6/14/15 7:15:17
Data1 = fits.open("../Data/lcn701010_x1dsum.fits")
#File says that the date obs was 6/14/15 7:47:56
Data2 = fits.open("../Data/lcn701020_x1dsum.fits")
#File says that the date obs was 6/14/15 9:12:04
Data3 = fits.open("../Data/lcn701030_x1dsum.fits")
#------------------------------------------------------------------------------
#Data Extraction
#headNa  = DataN[0].header headNb  = DataN[1].header where N is an int; was removed to save space 
#Headers are used to give more information and details about what is in the Fits files
TbData1 = Data1[1].data
Data1.close()
TbData2 = Data2[1].data
Data2.close()
TbData3 = Data3[1].data
Data3.close()
#------------------------------------------------------------------------------
#Extracting data into useable lists
wavelength1=TbData1['wavelength']
flux1=TbData1['flux']
error1=TbData1['error']
#
wavelength2=TbData2['wavelength']
flux2=TbData2['flux']
error2=TbData2['error']
#
wavelength3=TbData3['wavelength']
flux3=TbData3['flux']
error3=TbData3['error']
#------------------------------------------------------------------------------
#Concatenating data sets
#1162.9A - 1479.5A
TWL1=concatenate((wavelength1[1],wavelength1[0]),axis=0)
TF1=concatenate((flux1[1],flux1[0]),axis=0)
ER1=concatenate((error1[1],error1[0]),axis=0)
#892.4A-1208A
TWL2=concatenate((wavelength2[1],wavelength2[0]),axis=0)
TF2=concatenate((flux2[1],flux2[0]),axis=0)
ER2=concatenate((error2[1],error2[0]),axis=0)
#1374.4A-1763.1A
TWL3=concatenate((wavelength3[1],wavelength3[0]),axis=0)
TF3=concatenate((flux3[1],flux3[0]),axis=0)
ER3=concatenate((error3[1],error3[0]),axis=0)
#------------------------------------------------------------------------------
#SWAPPING VARIBLE NAMES SO THEY ARE IN ORDER OF INCREASING WAVELENGTH
TWL1,TWL2=TWL2,TWL1
TF1,TF2=TF2,TF1
ER1,ER2=ER2,ER1
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Smooth the data for easy visualization
TFS1=smooth(TF1,7)   # Use an odd number between 3 and 11, never more than 15 
TFS2=smooth(TF2,7)
TFS3=smooth(TF3,7)


#------------------------------------------------------------------------------
# Begining of Normalization
#------------------------------------------------------------------------------
#Creating point selection for first order polynomial
#Search for flat regious with no absorption and emission points for normalization
Grating_1_polyfitpoints=[make_polyfit_point(TWL1,TFS1,1119,1121),
                         make_polyfit_point(TWL1,TFS1,1161,1162)]
Grating_2_polyfitpoints=[make_polyfit_point(TWL2,TFS2,1178,1182),
                         make_polyfit_point(TWL2,TFS2,1350,1360),
                         make_polyfit_point(TWL2,TFS2,1395,1402),
                         make_polyfit_point(TWL2,TFS2,1435,1445)]
Grating_3_polyfitpoints=[make_polyfit_point(TWL3,TFS3,1425,1460),
                         make_polyfit_point(TWL3,TFS3,1500,1525),
                         make_polyfit_point(TWL3,TFS3,1590,1605),
                         make_polyfit_point(TWL3,TFS3,1675,1720)]
# For different epochs, input different ranges that represent the continuum. 
#------------------------------------------------------------------------------
#Determines the y (flux) value of the x values chosen above
x_poly_1=[item[0]for item in Grating_1_polyfitpoints]
y_poly_1=[item[1]for item in Grating_1_polyfitpoints]
x_poly_2=[item[0]for item in Grating_2_polyfitpoints]
y_poly_2=[item[1]for item in Grating_2_polyfitpoints]
x_poly_3=[item[0]for item in Grating_3_polyfitpoints]
y_poly_3=[item[1]for item in Grating_3_polyfitpoints]

#Creates first order polynomial fitting those two points
best_fit_poly_1=(polyfit(x_poly_1,y_poly_1,1))
best_fit_poly_2=(polyfit(x_poly_2,y_poly_2,1))
best_fit_poly_3=(polyfit(x_poly_3,y_poly_3,1))
#------------------------------------------------------------------------------
#Graphical Tool
#Used to create a list of points that represent a normalization curves

# What is this and why is it commented? ???
"""
xping1=arange(TWL1[0],TWL1[-1],.05)
yping1=[]
for i in range(len(xping1)):
    yping1.append(best_fit_poly_1[0]*(xping1[i])+ best_fit_poly_1[1])

xping2=arange(TWL2[0],TWL2[-1],.05)
yping2=[]
for i in range(len(xping2)):
    yping2.append(best_fit_poly_2[0]*(xping2[i])+ best_fit_poly_2[1])

xping3=arange(TWL3[0],TWL3[-1],.05)
yping3=[]
for i in range(len(xping3)):
    yping3.append(best_fit_poly_3[0]*(xping3[i])+ best_fit_poly_3[1])
"""

#plt.plot(TWL1,TFS1)

#------------------------------------------------------------------------------
#Second Order poly fit point selection for the PV region between gratings
Grating_1_2ndOrderPolyfit_Points=[  make_polyfit_point(TWL1,TFS1,1174.4,1175.4),
                                    make_polyfit_point(TWL1,TFS1,1191.2,1192.4),
                                    make_polyfit_point(TWL1,TFS1,1194.3,1196.0),
                                    make_polyfit_point(TWL1,TFS1,1184.0,1186.0),
                                    make_polyfit_point(TWL1,TFS1,1181.0,1182.0)]
Grating_2_2ndOrderPolyfit_Points=[  make_polyfit_point(TWL2,TFS2,1174.4,1175.4),
                                    make_polyfit_point(TWL2,TFS2,1191.2,1192.4),
                                    make_polyfit_point(TWL2,TFS2,1194.3,1196.0),
                                    make_polyfit_point(TWL2,TFS2,1184.0,1186.0),
                                    make_polyfit_point(TWL2,TFS2,1181.0,1182.0),
                                    make_polyfit_point(TWL2,TFS2,1204.0,1205.0)]
#------------------------------------------------------------------------------
#Determines the y (flux) value of the x values chosen above
x_2ndpoly_1=[item[0]for item in Grating_1_2ndOrderPolyfit_Points]
y_2ndpoly_1=[item[1]for item in Grating_1_2ndOrderPolyfit_Points]
x_2ndpoly_2=[item[0]for item in Grating_2_2ndOrderPolyfit_Points]
y_2ndpoly_2=[item[1]for item in Grating_2_2ndOrderPolyfit_Points]

#Creates Second Order Polynomial from seleted points
best_fit_2ndpoly_1=(polyfit(x_2ndpoly_1,y_2ndpoly_1,2))
best_fit_2ndpoly_2=(polyfit(x_2ndpoly_2,y_2ndpoly_2,2))

#------------------------------------------------------------------------------
#Normalizes for each Grating, special care taken for Intersection of Grating 1 and 2

# The two points used are to switch from fitting a first-order polynomium to a second order polynomium

# PRH: Typically, we divide first by first order and later correct with a second order...

#From mathematica, the polynomial intersects around 1167.21 and 1199.84

#Remember to change this every time the point selection changes
#------------------------------------------------------------------------------

# Just dividing by first order:
#For Grating 1
Normal1_TF1=[]
Normal1_ER1=[]
for i in range(len(TWL1)):
    Normal1_TF1.append(TF1[i]/(best_fit_poly_1[0]*(TWL1[i])+ best_fit_poly_1[1]))
    Normal1_ER1.append(ER1[i]/(best_fit_poly_1[0]*(TWL1[i])+ best_fit_poly_1[1]))

#For Grating 2
Normal1_TF2=[]
Normal1_ER2=[]
for i in range(len(TWL2)):
    Normal1_TF2.append(TF2[i]/(best_fit_poly_2[0]*(TWL2[i])+ best_fit_poly_2[1]))
    Normal1_ER2.append(ER2[i]/(best_fit_poly_2[0]*(TWL2[i])+ best_fit_poly_2[1]))

#For Grating 3
Normal1_TF3=[]
Normal1_ER3=[]
for i in range(len(TWL3)):
    Normal1_TF3.append(TF3[i]/(best_fit_poly_3[0]*(TWL3[i])+ best_fit_poly_3[1]))
    Normal1_ER3.append(ER3[i]/(best_fit_poly_3[0]*(TWL3[i])+ best_fit_poly_3[1]))


#  Dividing by first and second order around overlap range of gratings 1 and 2

#For Grating 1:
Normal_TF1=[]
Normal_ER1=[]
for i in range(0,find_index(TWL1,closest_value(TWL1,1167.23))):                #First part, First order poly  
    Normal_TF1.append(TF1[i]/(best_fit_poly_1[0]*(TWL1[i])+ best_fit_poly_1[1]))
    Normal_ER1.append(ER1[i]/(best_fit_poly_1[0]*(TWL1[i])+ best_fit_poly_1[1]))

for i in range(find_index(TWL1,closest_value(TWL1,1167.23)),len(TWL1)):        #Second part, Second order poly
    Normal_TF1.append(TF1[i]/(best_fit_2ndpoly_1[0]*(TWL1[i])**2+ best_fit_2ndpoly_1[1]*(TWL1[i])+ best_fit_2ndpoly_1[2]))
    Normal_ER1.append(ER1[i]/(best_fit_2ndpoly_1[0]*(TWL1[i])**2+ best_fit_2ndpoly_1[1]*(TWL1[i])+ best_fit_2ndpoly_1[2]))

#For Grating 2
Normal_TF2=[]
Normal_ER2=[]
for i in range(0,find_index(TWL2,closest_value(TWL2,1199.84))):                #First part, Second Order Polynomial
    Normal_TF2.append(TF2[i]/(best_fit_2ndpoly_2[0]*(TWL2[i])**2+ best_fit_2ndpoly_2[1]*(TWL2[i])+ best_fit_2ndpoly_2[2]))
    Normal_ER2.append(ER2[i]/(best_fit_2ndpoly_2[0]*(TWL2[i])**2+ best_fit_2ndpoly_2[1]*(TWL2[i])+ best_fit_2ndpoly_2[2]))
for i in range(find_index(TWL2,closest_value(TWL2,1199.84)),len(TWL2)):        #Second part, First Order Polynomial
    Normal_TF2.append(TF2[i]/(best_fit_poly_2[0]*(TWL2[i])+ best_fit_poly_2[1]))
    Normal_ER2.append(ER2[i]/(best_fit_poly_2[0]*(TWL2[i])+ best_fit_poly_2[1]))

#For Grating 3
Normal_TF3 = Normal1_TF3
Normal_ER3 = Normal1_ER3

#-----------------------------------------------------------------------------
    #Begining of combining Gratings
#------------------------------------------------------------------------------
# all Flux below 2% of the normalization line is equated to 0 : Wrong
# PRH: This doesn't seem right: you cannot change actual values, except for plotting purposes.

#remove_small(Normal_TFS1,0.02)
#remove_small(Normal_TFS2,0.02)
#remove_small(Normal_TFS3,0.02)

#------------------------------------------------------------------------------
# Interpolation:

#Finding stepsizes for each grating; then sets the wavelength array to the largest stepsize of the two
TWL1_STEPSIZE=TWL1[1]-TWL1[0]
TWL2_STEPSIZE=TWL2[1]-TWL2[0]
TWL3_STEPSIZE=TWL3[1]-TWL3[0]

if(TWL1_STEPSIZE<=TWL2_STEPSIZE):
    TWL_1_2_STEPSIZE=TWL2_STEPSIZE
else:
    TWL_1_2_STEPSIZE=TWL1_STEPSIZE
    
if(TWL_1_2_STEPSIZE<=TWL3_STEPSIZE):
    TWL_1_2_3_STEPSIZE=TWL3_STEPSIZE
else:
    TWL_1_2_3_STEPSIZE=TWL_1_2_STEPSIZE


# Create arrays 12 and 123 with the initial and ending points as before, but the largest stepsize of the ones combined
TWL1_TWL2= arange(TWL1[0],TWL2[-1],TWL_1_2_STEPSIZE)
TWL1_TWL2_TWL3= arange(TWL1[0],TWL3[-1],TWL_1_2_3_STEPSIZE)

#Interpolation of flux and error

Iflux1=interp(TWL1_TWL2,TWL1,Normal_TF1,left=0,right=0)
Ierror1=interp(TWL1_TWL2,TWL1,Normal_ER1,left=0,right=0)
Iflux2=interp(TWL1_TWL2,TWL2,Normal_TF2,left=0,right=0)
Ierror2=interp(TWL1_TWL2,TWL2,Normal_ER2,left=0,right=0)


Iflux3=interp(TWL1_TWL2_TWL3,TWL3,Normal_TF3,left=0,right=0)
Ierror3=interp(TWL1_TWL2_TWL3,TWL3,Normal_ER3,left=0,right=0)


#------------------------------------------------------------------------------
# Edges of spectra go to zero and affect the combination later, so we cut the edges.

# Inspect visually the regions using plot(TWL1[31950:32000],TF1[31950:32000]), for example,
# to decide the wavelength regions that need to be cut. 

CUT_GRATINGS=True
if( CUT_GRATINGS):# Enter the wavelength of the areas you would like to cut off, or "None"(no quotes)
    TWL1_Start_Cut=None  # Will cut off all values before this 
    TWL1_End_Cut=1200.0  # Will cut off all values after this 
    TWL2_Start_Cut=1170.8
    TWL2_End_Cut=1467.3 
    TWL3_Start_Cut=1390.0
    TWL3_End_Cut=None
   

#Selecting Index for usage when combining Gratings. <-- remove, left for teaching purposes    
# Assigns indexes to the start and end wavelength cuts, if set. If not, None. 
    
    if not TWL1_Start_Cut==None:TWL1_Start_Cut_index=find_index_from_value(TWL1,TWL1_Start_Cut)
    else:
        TWL1_Start_Cut_index=None
    if not TWL1_End_Cut==None:  TWL1_End_Cut_index=find_index_from_value(TWL1,TWL1_End_Cut)
    else:
        TWL1_End_Cut_index=None
    if not TWL2_Start_Cut==None:TWL2_Start_Cut_index=find_index_from_value(TWL2,TWL2_Start_Cut)
    else:
        TWL2_Start_Cut_index=None
    if not TWL2_End_Cut==None:  TWL2_End_Cut_index=find_index_from_value(TWL2,TWL2_End_Cut)
    else:
        TWL2_End_Cut_index=None
    if not TWL3_Start_Cut==None:TWL3_Start_Cut_index=find_index_from_value(TWL3,TWL3_Start_Cut)
    else:
        TWL3_Start_Cut_index=None
    if not TWL3_End_Cut==None:  TWL3_End_Cut_index=find_index_from_value(TWL3,TWL3_End_Cut)
    else:
        TWL3_End_Cut_index=None
#------------------------------------------------------------------------------
# Use the indexes above to the normalized flux and the error spectrum 
#In this Region Michael is using If and else statements to determing where to cut, from given selection above <-- remove, left for teaching purposes
# The sentence above is very wrong. He is not doing that, and we can tell already he uses if/then/else statements by looking. Do not comment that.     
    
# Never change values of data...! Create a new array where you do that, but not the main arrays. 

    Normal_TF1_plot = Normal_TF1
    Normal_TF2_plot = Normal_TF2
    Normal_TF3_plot = Normal_TF3

# The ER array was not normalized...! That is huge mistake. Corrected above.
# Assigning random numbers to the original error value is not recommended
# --> I removed all the ER=1 assignations.     
    
    if not TWL1_Start_Cut==None and not TWL1_End_Cut==None :
        for i in range(len(TWL1)):  # This Loop sets the cut portion of the flux values to zero, not used in the code but it makes plotting the gratings not show the cut portion <-- PRH: What do you mean "not used in the code"? It is changing the Normal_TFS1 and the ERR values!!!
# <-- remove comment above, left for teaching purposes
            if  i<=TWL1_Start_Cut_index:
                Normal_TF1_plot[i]=0
            if i>TWL1_End_Cut_index:
                Normal_TF1_plot[i]=0
    elif not TWL1_Start_Cut==None and TWL1_End_Cut==None:
        for i in range(TWL1_Start_Cut_index):
                Normal_TF1_plot[i]=0
    elif TWL1_Start_Cut==None and not TWL1_End_Cut==None:
        for i in range(TWL1_End_Cut_index,len(TWL1)):
                Normal_TF1_plot[i]=0
                
    if not TWL2_Start_Cut==None and not TWL2_End_Cut==None :
        for i in range(len(TWL2)): # This Loop sets the cut portion of the flux values to zero for plotting purposes
            if  i<=TWL2_Start_Cut_index:
                Normal_TF2_plot[i]=0
            if i>TWL2_End_Cut_index:
                Normal_TF2_plot[i]=0
    elif not TWL2_Start_Cut==None and TWL2_End_Cut==None:
        for i in range(TWL2_Start_Cut_index):
                Normal_TF2_plot[i]=0
    elif TWL2_Start_Cut==None and not TWL2_End_Cut==None:
        for i in range(TWL2_End_Cut_index,len(TWL2)):
                Normal_TF2_plot[i]=0
    if not TWL3_Start_Cut==None and not TWL3_End_Cut==None :
        for i in range(len(TWL3)):#This Loop removes the cut portion of the flux values you have set, not used in the code but it makes plotting the gratings not show the cut portion
            if  i<=TWL3_Start_Cut_index:
                Normal_TF3_plot[i]=0
            if i>TWL3_End_Cut_index:
                Normal_TF3_plot[i]=0
    elif not TWL3_Start_Cut==None and TWL3_End_Cut==None:
        for i in range(TWL3_Start_Cut_index):
                Normal_TF3_plot[i]=0
    elif TWL3_Start_Cut==None and not TWL3_End_Cut==None:
        for i in range(TWL3_End_Cut_index,len(TWL3)):
                Normal_TF3_plot[i]=0

# Correct badpixels for plotting purposes
badpix(Normal_TF2_plot,5180,5405,0)


#------------------------------------------------------------------------------
#Remove Zero Error
# I think this is dangerous. Removed. Check how much it is needed and whether we can do something different. 
#TME1=remove_zero_error(Ierror1)
#TME2=remove_zero_error(Ierror2)
#TME3=remove_zero_error(Ierror3)
#------------------------------------------------------------------------------
#Reconstructing Continuum using data and error  
       
Combined_TF1_TF2=[]
for i in range(0,len(TWL1_TWL2)):
    if (Iflux1[i]==0 and Iflux2[i]==0).all():
            Combined_TF1_TF2.append(0.0)
    elif(Iflux1[i]==0 and Iflux2[i]!=0):
        Combined_TF1_TF2.append(Iflux2[i])
    elif(Iflux1[i]!=0 and Iflux2[i]==0):
        Combined_TF1_TF2.append(Iflux1[i])
        
    elif(Ierror1[i]==0 and Ierror2[i]!=0):
        Combined_TF1_TF2.append(Iflux2[i])
    elif(Ierror1[i]!=1 and Ierror2[i]==1):
        Combined_TF1_TF2.append(Iflux1[i])
    elif(Ierror1[i]==1 and Ierror2[i]==1):
        Combined_TF1_TF2.append(0.0)
    else:  
        
        if Ierror1[i]!=0:
            weight1=(1./Ierror1[i])
        else:
            weight1=0
        if Ierror2[i]!=0:
            weight2=(1./Ierror2[i])
        else:
            weight2=0   
        
        spec1=Iflux1[i]
        spec2=Iflux2[i]
        Combined_TF1_TF2.append((weight1*spec1+weight2*spec2)/(weight1+weight2))  

        print(spec1,spec2)
        print((weight1*spec1+weight2*spec2)/(weight1+weight2))
        

# Check this below. Bad pixels should have removed from the original ones and only for plotting. 
badpix(Combined_TF1_TF2,32323,32550,0)



# It first combines and then interpolates to the TWL1_TWL2_TWL3 wavelength.        
Iflux_1_2=interp(TWL1_TWL2_TWL3,TWL1_TWL2,Combined_TF1_TF2,left=0,right=0)

Error_TF1_TF2=[]
for i in range(len(TWL1_TWL2)):
    if (Ierror1[i]!=0 and Ierror2[i]==0).all():
        Error_TF1_TF2.append(Ierror1[i])
    elif (Ierror1[i]==0 and Ierror2[i]!=0).all():
        Error_TF1_TF2.append(Ierror2[i])
    elif (Ierror1[i]==0 and Ierror2[i]==0).all():
        Error_TF1_TF2.append(1)            
    else:#if neither of the error values are one then take the weighted average of the error values.
        Error_TF1_TF2.append(sqrt((Ierror1[i])**2+(Ierror2[i])**2))   

Ierror_1_2=interp(TWL1_TWL2_TWL3,TWL1_TWL2,Error_TF1_TF2,left=0,right=0)

#TME_1_2=remove_zero_error(Ierror_1_2)

Combined_TF1_TF2_TF3=[]
for i in range(0,len(TWL1_TWL2_TWL3)):  
    if(Iflux_1_2[i]==0 and Iflux3[i]==0):
        Combined_TF1_TF2_TF3.append(0.0)
    elif(Iflux_1_2[i]==0 and Iflux3[i]!=0):
        Combined_TF1_TF2_TF3.append(Iflux3[i])
    elif(Iflux_1_2[i]!=0 and Iflux3[i]==0):
        Combined_TF1_TF2_TF3.append(Iflux_1_2[i])
        
    elif(Ierror_1_2[i]==0 and Ierror3[i]!=0):
        Combined_TF1_TF2_TF3.append(Iflux3[i])
    elif(Ierror_1_2[i]!=0 and Ierror3[i]==0):
        Combined_TF1_TF2_TF3.append(Iflux_1_2[i])
    elif(Ierror_1_2[i]==0 and Ierror3[i]==0):
        Combined_TF1_TF2_TF3.append(0.0)
    
    else:
        
        if Ierror_1_2[i]!=0:
            weight1=(1./Ierror_1_2[i])
        else:
            weight1=0
        if Ierror2[i]!=0:
            weight2=(1./Ierror_1_2[i])
        else:
            weight2=0  
        
        spec1=Iflux_1_2[i]
        spec2=Iflux3[i]
        Combined_TF1_TF2_TF3.append((weight1*spec1+weight2*spec2)/(weight1+weight2))
        
        print(spec1,spec2)
        print((weight1*spec1+weight2*spec2)/(weight1+weight2))
        
        
Error_TF1_TF2_TF3=[]
for i in range(len(TWL1_TWL2_TWL3)):
    if (Ierror_1_2[i]!=0 and Ierror3[i]==0).all():
        Error_TF1_TF2_TF3.append(Ierror_1_2[i])
    elif (Ierror_1_2[i]==0 and Ierror3[i]!=0).all():
        Error_TF1_TF2_TF3.append(Ierror3[i])
    elif (Ierror_1_2[i]==0 and Ierror3[i]==0).all():
        Error_TF1_TF2_TF3.append(1)            
    else:#if neither of the error values are one then take the weighted average of the error values.
        Error_TF1_TF2_TF3.append(sqrt((Ierror_1_2[i])**2+(Ierror3[i])**2))  
#------------------------------------------------------------------------------
    #Name change for final normalized and combined data set
#------------------------------------------------------------------------------

Final_e_spectrum=Error_TF1_TF2_TF3
Final_x_spectrum=TWL1_TWL2_TWL3
Final_y_spectrum=Combined_TF1_TF2_TF3

SFinal_e_spectrum=smooth(Final_e_spectrum,7)
SFinal_y_spectrum=smooth(Final_y_spectrum,7)


a=20000
b=71163

plt.plot(Final_x_spectrum[a:b],SFinal_y_spectrum[a:b])
plt.plot(Final_x_spectrum[a:b],SFinal_e_spectrum[a:b])


