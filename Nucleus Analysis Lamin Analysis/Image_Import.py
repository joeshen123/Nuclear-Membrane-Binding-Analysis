from nd2reader.reader import ND2Reader
import numpy as np
import time
from tkinter import ttk
from tkinter import simpledialog
from tkinter import filedialog
from tkinter import messagebox
import tkinter as tk 
import matplotlib.pyplot as plt
from skimage import io
import warnings
import h5py
from colorama import Fore
from tqdm import tqdm
import os 
import glob

#Ignore warnings issued by skimage through conversion to uint8
warnings.simplefilter("ignore",UserWarning)
warnings.simplefilter("ignore",RuntimeWarning)

# Use tkinter to interactively select files to import
root = tk.Tk()
root.withdraw()

root.directory = filedialog.askdirectory()
Directory_name = root.directory

answer = messagebox.askyesno("Question","Does the image contain multiple z stacks?")

def find_first_value_by_key(d, key_name):
    def recurse(subdict):
        if isinstance(subdict, dict):
            if key_name in subdict:
                return subdict[key_name]
            for value in subdict.values():
                result = recurse(value)
                if result is not None:
                    return result
        elif isinstance(subdict, list):
            for item in subdict:
                result = recurse(item)
                if result is not None:
                    return result
        return None

    return recurse(d)



def Z_Stack_Images_Extractor(Image_Sequence, fields_of_view, z_answer):
    Channel_list = Image_Sequence.metadata['channels']
    print(Channel_list)
    # Default assignments
    DsRed_Channel = 0
    GFP_Channel = 0
    BFP_Channel = None
    if len(Channel_list) == 1:
        pass  # Defaults already set

    elif len(Channel_list) == 2:
        if 'DSRED' in Channel_list[0] or '561' in Channel_list[0]:
            DsRed_Channel = 0
            GFP_Channel = 1
        else:
            DsRed_Channel = 1
            GFP_Channel = 0

    elif len(Channel_list) == 3:
        dsred_idx = next((i for i, ch in enumerate(Channel_list) if 'DSRED' in ch or '561' in ch), None)
        gfp_idx   = next((i for i, ch in enumerate(Channel_list) if 'GFP' in ch or '488' in ch), None)

        if dsred_idx is not None and gfp_idx is not None:
            DsRed_Channel = dsred_idx
            GFP_Channel = gfp_idx
        elif dsred_idx is not None:
            DsRed_Channel = dsred_idx
            GFP_Channel = (set([0, 1, 2]) - set([dsred_idx])).pop()
        else:
            DsRed_Channel = (set([0, 1, 2]) - set([gfp_idx])).pop()
            GFP_Channel = gfp_idx

    elif len(Channel_list) == 4:
        bfp_idx   = next((i for i, ch in enumerate(Channel_list) if 'DAPI' in ch or '405' in ch), None)
        dsred_idx = next((i for i, ch in enumerate(Channel_list) if 'DSRED' in ch or '561' in ch), None)
        gfp_idx   = next((i for i, ch in enumerate(Channel_list) if 'GFP' in ch or '488' in ch), None)

        DsRed_Channel = dsred_idx if dsred_idx is not None else 0
        GFP_Channel = gfp_idx if gfp_idx is not None else 1
        BFP_Channel = bfp_idx if bfp_idx is not None else 2

        print(DsRed_Channel)
        print(GFP_Channel)
        print(BFP_Channel)

    
    time_series = Image_Sequence.sizes['t']
   
    if z_answer:
      z_stack = Image_Sequence.sizes['z']
   
    GFP_Stack = []
    DsRed_Stack = []
    BFP_Stack = []

    n = 0

    for time in range(time_series):

     if z_answer:
       z_stack_dsRed = [] 
       z_stack_GFP = []
       z_stack_BFP = []

       for z_slice in range(z_stack):
          dsRed_slice = Image_Sequence.get_frame_2D(c=DsRed_Channel, t=time, z=z_slice, v=fields_of_view)
          GFP_slice = Image_Sequence.get_frame_2D(c=GFP_Channel, t=time, z=z_slice, v=fields_of_view)
          z_stack_dsRed.append(dsRed_slice)
          z_stack_GFP.append(GFP_slice)
          
          if BFP_Channel != None:
             BFP_slice = Image_Sequence.get_frame_2D(c=BFP_Channel, t=time, z=z_slice, v=fields_of_view)
             z_stack_BFP.append(BFP_slice)
          
          
     else:
        z_stack_dsRed = Image_Sequence.get_frame_2D(c=DsRed_Channel, t=time, v=fields_of_view)
        z_stack_GFP = Image_Sequence.get_frame_2D(c=GFP_Channel, t=time, v=fields_of_view)
 
        if BFP_Channel != None:
           z_stack_BFP= Image_Sequence.get_frame_2D(c=BFP_Channel, t=time, z=z_slice, v=fields_of_view)

     DsRed_Stack.append(z_stack_dsRed)
     GFP_Stack.append(z_stack_GFP)

     if BFP_Channel != None:
       BFP_Stack.append(z_stack_BFP)

     
    DsRed_Stack = np.array(DsRed_Stack)
    GFP_Stack = np.array(GFP_Stack)

    if BFP_Channel != None:
       BFP_Stack = np.array(BFP_Stack)

       return (DsRed_Stack, GFP_Stack, BFP_Stack)
    
    else:
       
       return (DsRed_Stack, GFP_Stack)




os.chdir(Directory_name)

df_filenames = glob.glob('*.nd2' )

# create progress bar
pb = tqdm(range(len(df_filenames)), bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.RED, Fore.RESET))

for img_num in pb:
   pb.set_description("Converting nd2 to hdf5 files")
   
   img = df_filenames[img_num]
   
   Image_Sequence = ND2Reader(img)
   #print(Image_Sequence.parser._raw_metadata.image_metadata)
   #Get z_step and x and y pixels
   z_step = find_first_value_by_key(Image_Sequence.parser._raw_metadata.image_metadata,b'dZStep')
   x_pixel= Image_Sequence.metadata['pixel_microns']
   y_pixel= Image_Sequence.metadata['pixel_microns']

   voxel_info = (x_pixel, y_pixel, z_step)

   print(voxel_info)
   FOV_list = Image_Sequence.metadata['fields_of_view']



   File_save_names = '.'.join(img.split(".")[:-1])

   for fov in range(len(FOV_list)):
      
      try:
        Channel_561, Channel_488 = Z_Stack_Images_Extractor(Image_Sequence,fields_of_view=fov,z_answer=answer)

        Image_Name='{File_Name}_{num}.hdf5'.format(File_Name = File_save_names, num = fov + 1)

        with h5py.File(Image_Name, "w") as f:
           f.create_dataset('488 Channel', data = Channel_488, compression = 'gzip')
           f.create_dataset('561 Channel', data = Channel_561, compression = 'gzip')
           f.create_dataset('voxel_info', data=voxel_info, compression='gzip')
      
      except:
         Channel_561, Channel_488, Channel_405 = Z_Stack_Images_Extractor(Image_Sequence,fields_of_view=fov,z_answer=answer)

         Image_Name='{File_Name}_{num}.hdf5'.format(File_Name = File_save_names, num = fov + 1)

         with h5py.File(Image_Name, "w") as f:
           f.create_dataset('488 Channel', data = Channel_488, compression = 'gzip')
           f.create_dataset('561 Channel', data = Channel_561, compression = 'gzip')
           f.create_dataset('405 Channel', data = Channel_405, compression = 'gzip')
           f.create_dataset('voxel_info', data=voxel_info, compression='gzip')
      

'''my_filetypes = [('all files', '.*'),('Movie files', '.nd2')]

Image_Stack_Path = filedialog.askopenfilename(title='Please Select a Movie', filetypes = my_filetypes)
'''


# Define a function to convert time series of ND2 images to a numpy list of 
# images (t,z,y,x).





'''

Image_Sequence = ND2Reader(Image_Stack_Path)
FOV_list = Image_Sequence.metadata['fields_of_view']

GUV_Image_list = []
Intensity_list = []

for fov in range(len(FOV_list)):
   GUV_Images, Image_Intensity = Z_Stack_Images_Extractor(Image_Stack_Path,fields_of_view=fov,z_answer=answer)
   GUV_Image_list.append(GUV_Images)
   Intensity_list.append(Image_Intensity)


File_save_names = '.'.join(Image_Stack_Path.split(".")[:-1])

for n in range(len(FOV_list)):
   GUV_Image_Name='{File_Name}_{num}.hdf5'.format(File_Name = File_save_names, num = n + 1)
   
   GUV_Images = GUV_Image_list[n]
   Image_Intensity = Intensity_list[n]

   with h5py.File(GUV_Image_Name, "w") as f:
      f.create_dataset('488 Channel', data = Image_Intensity, compression = 'gzip')
      f.create_dataset('561 Channel', data = GUV_Images, compression = 'gzip')
'''