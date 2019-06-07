#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import datetime
import fnmatch
import shutil



####################################### User Input #####################################################
image_datapath = '/Users/aita2652/Cell Profiler Material/15MayLiveCell/2019-05-15/93'
print('Image Datapath: ',image_datapath)
output_datapath = '/Users/aita2652/Cell Profiler Material/15MayLiveCell/Python'
print('Image Destination path: ',output_datapath)
os.chdir(image_datapath)

 #example file name: "MCF7-40X-LIVE-pilot_C03_s26B47385C-DFA7-4039-ABBB-BCF4F2EC1260". Here 'C03' is the well number
    #  and 's2' is the site. Data after 's2' is not useful for us. 
    #U2OS-20X-BpAHDDB-60uM30uM15uM-LIVE_C03_s1_w101F8A2CB-878F-4158-AB72-2E26EC96FB92
    
unique_separator = '_' ##The unique separator present in the image names

#####################################---------XXX----------###########################################



################################# Metadata Text File Generator ########################################

text_filehandle = os.path.join(output_datapath,'Metadata file.txt') ##Generating the text file
f = open(text_filehandle, 'w')
f.write('------------------------The Metadata File------------------------\n\n\n')
f.write('Generated on: ')
f.write(str(datetime.datetime.now()))
f.write('\n\n\n\n')

#####################################---------XXX----------###########################################



################################ Main Program Body ###################################################

for filename in os.listdir(os.getcwd()): #looping over the files present in the provided dir
    
    imp = 0 #counter reset
    if(filename == '.DS_Store'): ##removes .DS_Store
            continue
            
    else:
        for imagefile in os.listdir(filename): #looping over the directory of the filenames!
            
            if(imagefile == '.DS_Store'): ##removes .DS_Store
                continue
                
            elif(fnmatch.fnmatch(imagefile, '*thumb*')): ##removes image thumbnails
                continue
                
            elif(fnmatch.fnmatch(imagefile, '*w1*')): ##ignores nucleus channel images
                continue
                
            else:
                print("File name:",imagefile)
                f.write('1. File name: ')
                f.write(imagefile)
                f.write('\n')
                
            
                metadata, well, site, channel = imagefile.split(unique_separator) #Spliting using the separator provided
                site = site[:2] #to slice out the unwanted text after the site
                channel = channel[:2] #to slice out the unwanted text after the site
                
                print("Metadata:",metadata)
                f.write('2. Metadata: ')
                f.write(metadata)
                f.write('\n')
                
                print("Well number:",well)
                f.write('3. Well number: ')
                f.write(well)
                f.write('\n')
                
                print("Site:",site)
                f.write('4. Site: ')
                f.write(site)
                f.write('\n')
                
                print("Channel:",channel)
                f.write('5. Channel: ')
                f.write(channel)
                f.write('\n')
                
                bb1, file_ext = os.path.splitext(imagefile) ##to extract the extension
                bb2, timepoint = filename.split('_') ##bb are temporary storage 
                                                    ##this extracts the time-point
                
                imp+=1 # counter to check if all the files are processed within the Time_Point folder! 
                print(timepoint)
                f.write('6. Time Point: ')
                f.write(timepoint)
                f.write('\n')
                    
                timepoint = timepoint.zfill(3) ##to make 3-digit timepoint~(001, 055, 100)
                new_image_name = '{}{}'.format(timepoint, file_ext) #appending file extension
            
                img_path_orig = os.path.join(os.getcwd(), filename, imagefile) ##original image path
                
                print("Original Image Path:", img_path_orig)
                f.write('7. Original Image Path: ')
                f.write(img_path_orig)
                f.write('\n')
                

            
            
           
                path_new = os.path.join(output_datapath, metadata, well, site) #new image path
                
                if(os.path.exists(path_new)):
                    img_path_new = os.path.join(path_new, new_image_name)
                    shutil.copy2(img_path_orig, img_path_new)
                    
                    print("New Image Path:", img_path_new)
                    f.write('8. New Image Path: ')
                    f.write(img_path_new)
                    f.write('\n')
                    
                    ####To avoid system errors!
                
                else:
                    os.makedirs(os.path.join(path_new)) #makes the new required directory
                    img_path_new = os.path.join(path_new, new_image_name)
                    shutil.copy2(img_path_orig, img_path_new)
                    
                    print("New Image Path:", img_path_new)
                    f.write('8. New Image Path: ')
                    f.write(img_path_new)
                    f.write('\n')
                    
                print("Iteration number:",imp)
                f.write('9. Iteration number: ')
                f.write(str(imp))
                f.write('\n')
                print("  ----   ")
                f.write('----------------\n')
                    
                    
                
                
f.close()


# In[ ]:




