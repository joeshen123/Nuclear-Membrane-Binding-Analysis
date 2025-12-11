import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure
from aicssegmentation.core.pre_processing_utils import intensity_normalization, image_smoothing_gaussian_3d, \
    edge_preserving_smoothing_3d, image_smoothing_gaussian_slice_by_slice
from aicssegmentation.core.MO_threshold import MO
from scipy import ndimage as ndi
from skimage.morphology import remove_small_objects, binary_erosion, binary_dilation, ball, disk, erosion, dilation
from aicssegmentation.core.utils import get_middle_frame
from skimage import transform, measure
import h5py
from tqdm import tqdm
import pandas as pd
from pandas import DataFrame
import warnings
import h5py
from skimage.filters import sobel, scharr, gaussian, median
from skimage.segmentation import watershed
from skimage.morphology import binary_closing, ball
from skimage.measure import label, regionprops, regionprops_table
from aicssegmentation.core.utils import hole_filling
import trackpy as tp
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import napari
from matplotlib.pyplot import cm
from tkinter import simpledialog
from skimage.restoration import rolling_ball
from skimage.restoration import ellipsoid_kernel
from colorama import Fore
from skimage.segmentation import clear_border
import os
from skimage.exposure import rescale_intensity
from skimage.filters import threshold_otsu
from skimage.morphology import binary_opening
from aicssegmentation.core.vessel import filament_2d_wrapper

os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'


# Make a function to get middle frame based on segmented area
def get_middle_frame_area(labelled_image_stack):
    max_area = 0
    max_n = 0
    for z in range(labelled_image_stack.shape[0]):
        img_slice = labelled_image_stack[z, :, :]
        area = np.count_nonzero(img_slice)

        if area >= max_area:
            max_area = area
            max_n = z

    return max_n


def get_3dseed_from_mid_frame(bw, stack_shape, mid_frame, hole_min, bg_seed=True):
    from skimage.morphology import remove_small_objects
    out = remove_small_objects(bw > 0, hole_min)

    out1 = label(out)
    stat = regionprops(out1)

    # build the seed for watershed
    seed = np.zeros(stack_shape)
    seed_count = 0
    if bg_seed:
        seed[0, :, :] = 1
        seed_count += 1

    for idx in range(len(stat)):
        py, px = np.round(stat[idx].centroid)
        seed_count += 1
        seed[mid_frame, int(py), int(px)] = seed_count

    return seed


# Clean out non-track objects and color by tracks of segmented image stack
def clean_segmentation(table, seg_image):
    original_label = table.label.to_list()
    original_label = [x * 1000 for x in original_label]
    new_label = table.particle.to_list()
    # new_label = [x+1 for x in new_label]

    label_convert = list(zip(original_label, new_label))

    labelled_image = seg_image.copy().astype(np.uint16) * 1000

    for num in np.unique(labelled_image):
        if num != 0 and num not in original_label:
            labelled_image[labelled_image == num] = 0

    for element in label_convert:
        original = element[0]
        new = element[1]
        labelled_image[labelled_image == original] = new

    return labelled_image


# Plot mid section area, volume and Protein binding
# Define a function to plot both section area and protein intensity from panda df
def plotting(df, marker=None,saving=False):
    static_canvas = FigureCanvas(Figure(figsize=(4, 1)))

    axes = static_canvas.figure.subplots(4, sharex=True)

    # list all remaining tracks
    all_labels = df.particle.unique()

    color = cm.tab20b(np.linspace(0, 1, len(all_labels)))

    axes[0].set_ylabel('Mid Section', fontsize=16, fontweight='bold')
    axes[1].set_ylabel('Volume', fontsize=16, fontweight='bold')
    axes[2].set_ylabel('Protein Binding', fontsize=16, fontweight='bold')
    axes[3].set_ylabel('NM Folds', fontsize=16, fontweight='bold')
    axes[3].set_xlabel('Frame', fontsize=16, fontweight='bold')

    num = 0

    for n in all_labels:
        df_subset = df[df.particle == n]
        frame_list = df_subset['frame'].tolist()
        area_list = df_subset['mid section area'].tolist()
        volume_list = df_subset['normalized volume'].tolist()
        # binding_list = df_subset['intensity ratio'].tolist()
        binding_list = df_subset['normalized intensity ratio'].tolist()
        Folds = df_subset['normalized folds'].tolist()
        axes[0].plot(frame_list, area_list, color=color[num], label=str(n), marker=marker)
        axes[1].plot(frame_list, volume_list, color=color[num], label=str(n), marker=marker)
        axes[2].plot(frame_list, binding_list, color=color[num], label=str(n), marker=marker)
        axes[3].plot(frame_list, Folds, color=color[num], label=str(n), marker=marker)
        num += 1

    axes[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    if saving == True:

       fig, axes = plt.subplots(2, figsize=(2, 1), sharex=True)
       fig.set_size_inches(14, 18)
 
       # list all remaining tracks
       all_labels = df.particle.unique()

       color = cm.tab20b(np.linspace(0, 1, len(all_labels)))

       
       axes[0].set_ylabel('normalized volume', fontsize=16, fontweight='bold')
       axes[1].set_ylabel('Protein Binding', fontsize=16, fontweight='bold')
       axes[1].set_xlabel('Frame', fontsize=16, fontweight='bold')

       num = 0

       for n in all_labels:
        df_subset = df[df.particle == n]
        frame_list = df_subset['frame'].tolist()
        #area_list = df_subset['mid section area'].tolist()
        volume_list = df_subset['normalized volume'].tolist()
        # binding_list = df_subset['intensity ratio'].tolist()
        binding_list = df_subset['normalized intensity ratio'].tolist()
        #Folds = df_subset['normalized folds'].tolist()
        axes[0].plot(frame_list, area_list, color=color[num], label=str(n), marker=marker)
        axes[0].plot(frame_list, volume_list, color=color[num], label=str(n), marker=marker)
        axes[1].plot(frame_list, binding_list, color=color[num], label=str(n), marker=marker)
        axes[3].plot(frame_list, Folds, color=color[num], label=str(n), marker=marker)
        num += 1

        #axes[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

        #plt.show()


    return static_canvas


# Function to view the segmentation intermediates and final results in Napari Viewer
def Napari_Viewer(raw_image,structure_img_smooth,intensity_image_stack_raw,intensity_image_stack,segmented_object_image,bw_img,filament_image,skeleton_image, seed_map_list,branch_point_list,contour_image_crop,segmented_image_crop,track_points,track_table,x,y,z):
    # Visualize in Napari
    viewer = napari.Viewer()
    viewer.add_image(raw_image, name='Raw LMNB1 Image', colormap='green', scale=[1, z, y, x],
                     blending='additive')
    viewer.add_image(structure_img_smooth, name='Smooth Image', colormap='green', scale=[1, z, y, x],
                     blending='additive')
    viewer.add_image(intensity_image_stack_raw, name='Intensity Image', colormap='red',
                     scale=[1, z, y, x], blending='additive',visible=False)
    viewer.add_image(intensity_image_stack, name='Intensity Image smooth', colormap='red',
                     scale=[1, z, y, x], blending='additive', visible=False)
    viewer.add_labels(segmented_object_image, name='Segmented Object', scale=[1, z, y, x], blending='additive',
                      visible=False)
    viewer.add_image(bw_img, name='Thresholded', scale=[1, z, y, x], blending='additive', visible=False)
    viewer.add_image(filament_image, name='Filament of NM wrinkles', colormap='inferno',
                     scale=[1, z, y, x], blending='additive', visible=False)
    viewer.add_image(skeleton_image, name='Skeleton of NM wrinkles', colormap='inferno',
                     scale=[1, z, y, x], blending='additive')
    peaks = np.nonzero(seed_map_list)
    viewer.add_points(np.array(peaks).T, name='peaks', size=5, face_color='red', scale=[1, z, y, x],
                      blending='additive', visible=False)
    folds = np.nonzero(branch_point_list)
    viewer.add_points(np.array(folds).T, name='NM Folds/Branch Points', size=20000, face_color= 'red', scale=[1, z, y, x],
                      blending='additive')
    viewer.add_image(contour_image_crop, name='contours crop', scale=[1, z, y, x], blending='additive',visible=False)
    viewer.add_labels(segmented_image_crop, name='Segmented Object Crop', scale=[1, z, y, x],
                      blending='additive', visible=False)
    viewer.add_tracks(track_points, name='tracks', scale=[1, z, y, x], blending='additive',
                      visible=False)

    # Plotting
    canvasobj = plotting(track_table, marker='o')
    viewer.window.add_dock_widget(canvasobj, area='right', name='Analysis Plot')
    
    napari.run()

## Make the 3D Cell Segmentation Pipeline

# Initiate the cell segmentation class
from skimage import filters


class cell_segment:
    def __init__(self, Time_lapse_image,sigma=1,mask=None, scale=1,scale_z=1):
        self.image = Time_lapse_image.copy()
        self.Time_pts = Time_lapse_image.shape[0]

        self.smooth_param = sigma

        self.bw_img = np.zeros(self.image.shape)
        self.structure_img_smooth = np.zeros(self.image.shape)
        self.segmented_object_image = np.zeros(self.image.shape, dtype=np.uint8)
        self.seed_map_list = np.zeros(self.image.shape, dtype=np.uint8)
        self.mask = mask
        self.scale = scale
        self.scale_z = scale_z
    
       # define function to crop the images
    def cropping_image (self):

        viewer = napari.Viewer()
        viewer.add_image(self.image, blending='additive', name="Raw Images", colormap='green', visible=True,
                         scale=(1, self.scale_z, self.scale, self.scale))
        rec_layer = viewer.add_shapes(ndim=4, shape_type='rectangle', edge_width=1.5,
                                       edge_color='b')

        napari.run()

        y1, _, y2, _ = rec_layer.data[0][:, 2]
        x1, x2, _, _ = rec_layer.data[0][:, 3]
        y1 = int(y1/self.scale)
        y2 = int(y2/self.scale)

        x1 = int(x1/self.scale)
        x2 = int(x2/self.scale)

        if y1 <=0:
            y1 = 0

        if y1 >= self.image.shape[2]:
            y1 = self.image.shape[2]

        if y2 <= 0:
            y2 = 0

        if y2 >= self.image.shape[2]:
            y2 = self.image.shape[2]

        if x1 <= 0:
            x1 = 0

        if x1 >= self.image.shape[3]:
            x1 = self.image.shape[3]

        if x2 <= 0:
            x2 = 0

        if x2 >= self.image.shape[3]:
            x2 = self.image.shape[3]

        temp_image = np.zeros(self.image.shape)
        temp_image[:, :, y1:y2, x1:x2] = self.image[:, :, y1:y2, x1:x2]

        print(self.image[:, :, y1:y2, x1:x2].shape)
        self.image = temp_image
        
    # define a function to apply normalization and smooth on Time lapse images
    def img_norm_smooth(self):
        pb = tqdm(range(self.Time_pts), bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.RED, Fore.RESET))
        for t in pb:
            pb.set_description("Smooth and background substraction")
            img = self.image[t].copy()
            struct_img = img

            self.structure_img_smooth[t] = image_smoothing_gaussian_3d(img, sigma=self.smooth_param)
            


    def threshold_Time(self):
        pb = tqdm(range(self.Time_pts), bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET))
        for t in pb:
            pb.set_description("Thresholding and Watershed Segmentation")
            # MO threhodling
            # Check if there is NAN data
            img_sum = np.sum(self.structure_img_smooth[t])
            Nan_check = np.isnan(img_sum)

            # If there is NAN in the data, skip analysis for this image
            if Nan_check == True:
                print("Nan image found")
                continue

            bw, object_for_debug = MO(self.structure_img_smooth[t], global_thresh_method='ave', extra_criteria=True,
                                      object_minArea=4000, return_object=True)
            
            #threshold = threshold_otsu(self.structure_img_smooth[t])
            #bw = self.structure_img_smooth[t] > threshold

            # Morphological operations to fill holes and remove small/touching objects
            bw = binary_closing(bw, footprint =np.ones((4, 4, 4)))
            bw = hole_filling(bw, 1, 40000, fill_2d=True)
            bw = remove_small_objects(bw > 0, min_size=8000, connectivity=1)
            bw = clear_border(bw, mask=self.mask)

            self.bw_img[t] = bw

            # Get middle frame
            mid_z = get_middle_frame_area(bw)

            if t > 0:
                mid_z_temp = get_middle_frame_area(bw)

                if abs(mid_z_temp - mid_z) < 4:
                    mid_z = mid_z_temp

            bw_mid_z = bw[mid_z, :, :]

            # Get seed map
            seed = get_3dseed_from_mid_frame(bw_mid_z, bw.shape, mid_z, 2000, bg_seed=False)
            
            edge = scharr(self.image[t])

            seg = watershed(edge, markers=label(seed), mask=bw, watershed_line=True)
            seg = clear_border(seg, mask=self.mask)

            seg = remove_small_objects(seg > 0, min_size=8000, connectivity=1)
            seg = hole_filling(seg, 1, 40000, fill_2d=True)
            final_seg = label(seg)

            self.segmented_object_image[t] = final_seg
            self.seed_map_list[t] = seed


## Track and quantify nucleus geometry (e.g. volume) and protein binding over time
class cell_tracking:
    def __init__(self, segmented_image_seq,smooth_image_stack, intensity_image_stack, smooth_sigma, x_vol, y_vol, z_vol):
        self.labelled_stack = segmented_image_seq
        self.smooth_image_stack = smooth_image_stack
        self.segmented_image_crop = self.labelled_stack.copy()
        self.contour_image_crop = np.zeros(self.labelled_stack.shape)
        self.filament_image_crop = np.zeros(self.labelled_stack.shape)
        self.skeleton_image_crop = np.zeros(self.labelled_stack.shape)
        self.branch_point_crop = np.zeros(self.labelled_stack.shape)
        self.t = segmented_image_seq.shape[0]
        self.positions_table = None
        self.intensity_image_stack_raw = intensity_image_stack
        self.intensity_image_stack = np.zeros(intensity_image_stack.shape)
        self.smooth_sigma = smooth_sigma
        self.x_vol = x_vol
        self.y_vol = y_vol
        self.z_vol = z_vol

    # Function to  smooth intesnity image stack
    def intensity_stack_smooth(self):
        # pb = tqdm(range(self.intensity_image_stack.shape[0]), bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.YELLOW, Fore.RESET))
        pb = tqdm(range(self.intensity_image_stack.shape[0]),
                  bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.YELLOW, Fore.RESET))
        for t in pb:
            pb.set_description("Intensity stack smooth")
            img = self.intensity_image_stack_raw[t].copy()
        
            self.intensity_image_stack[t] = image_smoothing_gaussian_3d(img, sigma=self.smooth_sigma)
            
    # function to create pandas table of cell attributes without tracking info
    def create_table_regions(self):
        from skimage.morphology import skeletonize
        positions = []
        pb = tqdm(range(self.t), bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.CYAN, Fore.RESET))
        for n in pb:
            pb.set_description("Create table")
            labelled_slice = self.labelled_stack[n]
            smooth_Image = self.smooth_image_stack[n]
            for region in measure.regionprops(labelled_slice):
                position = []

                z_pos = region.centroid[0]
                y_row = region.centroid[1]
                x_col = region.centroid[2]

                volume = region.area * (self.x_vol * self.y_vol * self.z_vol)

                nucleus_image = labelled_slice == region.label
                
                mid_z = get_middle_frame_area(nucleus_image)
               
                mid_nucleus_image = nucleus_image[mid_z, :, :]
                filament_image = np.zeros(mid_nucleus_image.shape)
                skeleton_image = np.zeros(mid_nucleus_image.shape)
                branch = np.zeros(mid_nucleus_image.shape)
                #plt.imshow(label_image)
                #plt.show()

                min_z, min_y, min_x, max_z, max_y, max_x = region.bbox

                cropped_image = smooth_Image[:, min_y:max_y, min_x:max_x]
                
                cropped_image_middle = cropped_image[mid_z]

                f2_param = [[1.6,0.2]]

                cropped_image_filament= filament_2d_wrapper(cropped_image_middle, f2_param)
                
                cropped_image_filament = hole_filling(cropped_image_filament, 1, 300, fill_2d=True)
                cropped_image_filament = remove_small_objects(cropped_image_filament > 0, min_size=50, connectivity=1)
                filament_image[min_y:max_y, min_x:max_x] = cropped_image_filament

                self.filament_image_crop[n,mid_z,:,:] += filament_image
                #plt.imshow(cropped_image_filament)
                #plt.show()
                
                from skimage.morphology import skeletonize
                from scipy.ndimage import convolve
        

                skeleton = skeletonize(cropped_image_filament)
                skeleton_image[min_y:max_y, min_x:max_x] = skeleton
                self.skeleton_image_crop[n,mid_z,:,:] += skeleton_image

                label_image = label(skeleton)
                props = regionprops(label_image)
                if not props:
                   return None

                kernel = np.ones((3, 3), dtype=int)
                neighbor_count = convolve(skeleton.astype(int), kernel, mode='constant', cval=0)
                neighbor_count = neighbor_count * skeleton  # restrict to skeleton pixels
                neighbor_count-= skeleton  # remove self count

                branch_points = (neighbor_count >= 3)

                labels = label(branch_points, connectivity=2)
                
                # Initialize empty array for condensed branch points
                condensed = np.zeros_like(branch_points, dtype=bool)
                
                for region1 in regionprops(labels):
                   coords = region1.coords              # array of [[row, col], ...]
                   centroid = region1.centroid          # (row_float, col_float)
                   # compute squared distances to float-centroid
                   d2 = np.sum((coords - centroid)**2, axis=1)
                   # pick the coordinate with minimal distance
                   r, c = coords[np.argmin(d2)]
                   condensed[r, c] = True
                
                 # 4) overwrite your branch_points mask
                branch_points = condensed
                #print(branch_points)
                #print(branch_points.shape)
                num_branch_points = np.count_nonzero(branch_points)
                branch[min_y:max_y, min_x:max_x] = branch_points
                self.branch_point_crop[n, mid_z,:,:] += branch
               
                mid_section_area = max([region.area * (self.x_vol * self.y_vol) for region in
                                        measure.regionprops(label(mid_nucleus_image))])
               

                # Draw contour
                segmented_image_shell = np.logical_xor(erosion(mid_nucleus_image, footprint=disk(2)),
                                                       erosion(mid_nucleus_image, footprint=disk(5)))
                

                
                                                       
                bg_segmented_image_shell = mid_nucleus_image
                bg_segmented_image_shell = np.logical_xor(erosion(mid_nucleus_image, footprint=disk(12)),
                                                          erosion(mid_nucleus_image, footprint=disk(15)))

                labels = region.label
                intensity_single = self.intensity_image_stack[n]
                intensity_image = intensity_single[mid_z]
                intensity_median = np.median(intensity_image[segmented_image_shell == True])

                intensity_background = np.median(intensity_image[bg_segmented_image_shell == True])
                 

                intensity_median_ratio = intensity_median / intensity_background
                
                
        
                position.append(x_col)
                position.append(y_row)
                position.append(z_pos)
                position.append(int(n))
                position.append(labels)

                position.append(volume)
                position.append(mid_section_area)
                position.append(intensity_median)
                position.append(intensity_median_ratio)
                position.append(num_branch_points)
                positions.append(position)

        self.positions_table = DataFrame(positions,
                                         columns=['x', 'y', 'z', "frame", 'label', 'volume', 'mid section area',
                                                  'intensity', 'intensity ratio','Number of Folds'])

    # function to track subsequent frame
    def tracking(self, s_range=130, stop=2, step=0.95, gap=1, pos=['x', 'y', 'z']):
        self.track_table = tp.link_df(self.positions_table, s_range, adaptive_stop=stop, adaptive_step=step, memory=gap,
                                      pos_columns=pos)
        print(self.track_table)
        # Filter out tracks with low number of frames. Add 1 to particle number to avoid 0.
        self.track_table = tp.filter_stubs(self.track_table, self.t - 5)
        self.track_table['particle'] = self.track_table['particle'] + 1

        # Extract track points for visualization
        track_df = self.track_table[['particle', 'frame', 'z', 'y', 'x']]
        track_df.index.names = ['Data']
        track_df.sort_values(by=['particle', 'frame'], inplace=True)
        self.track_points = track_df.values

    # function to crop out broken frame based on tracks (only on segmentation image)
    def crop_segmentation(self):
        tracks = np.unique(self.track_points[:, 0])

        remain_track_list = []
        for tk in tracks:
            remain_tracks = self.track_points[self.track_points[:, 0] == tk][:, 1]
            remain_track_list.append(remain_tracks)

        print(remain_track_list)
        # Find all frames (union) among different tracks
        if len(remain_track_list) > 1:
            union_tracks = set.union(*map(set, remain_track_list))
            union_tracks = list(map(int, union_tracks))
        elif len(remain_track_list) == 1:
            union_tracks = remain_track_list[0]
            union_tracks = list(map(int, union_tracks))
        else:
            union_tracks = list(np.arange(self.t))

        print(union_tracks)
        # loop through all time, only keep objects that are in tracks
        pb = tqdm(range(self.t), bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.GREEN, Fore.RESET))
        for t in pb:
            pb.set_description("Remove broken objects from segmentation")
            if t in union_tracks:
                table_t = self.track_table.loc[self.track_table.frame == t]
                seg_t = self.labelled_stack[t]
                self.segmented_image_crop[t] = clean_segmentation(table_t, seg_t)

            else:
                if t == 0:
                    table_t = self.track_table.loc[self.track_table.frame == min(union_tracks)]
                    seg_t = self.labelled_stack[min(union_tracks)]
                    self.segmented_image_crop[t] = clean_segmentation(table_t, seg_t)

                else:
                    self.segmented_image_crop[t] = self.segmented_image_crop[t - 1]

    # function to crop out broken frame according to tracks (on mid plane contours)
    def crop_contour(self):
        # loop through all time, only keep objects that are in tracks
        pb = tqdm(range(self.t), bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.MAGENTA, Fore.RESET))
        for t in pb:
            pb.set_description("Remove broken object from mid contour")
            slice_t = self.segmented_image_crop[t]

            for region in regionprops(slice_t):
                y_row = region.centroid[1]
                x_col = region.centroid[2]

                nucleus_image = slice_t == region.label
                mid_z = get_middle_frame_area(nucleus_image)

                mid_nucleus_image = nucleus_image[mid_z, :, :]

                # Draw contours

                segmented_image_shell = np.logical_xor(erosion(mid_nucleus_image, footprint=disk(2)),
                                                       erosion(mid_nucleus_image, footprint=disk(5)))
                bg_segmented_image_shell = np.logical_xor(erosion(mid_nucleus_image, footprint=disk(12)),
                                                          erosion(mid_nucleus_image, footprint=disk(15)))

                self.contour_image_crop[t, mid_z, :, :] += segmented_image_shell
                self.contour_image_crop[t, mid_z, :, :] += bg_segmented_image_shell

    # function to normalize binding intensity to 1st frame
    def binding_normalize(self):
        self.track_table['normalized intensity ratio'] = self.track_table['intensity ratio']
        self.track_table['normalized folds'] = self.track_table['Number of Folds']
        self.track_table['normalized volume'] = self.track_table['volume']

        all_labels = self.track_table.particle.unique()

        for label in all_labels:
            norm_value = self.track_table.loc[self.track_table.particle == label, 'intensity ratio'].iloc[0]
            self.track_table.loc[self.track_table.particle == label, 'normalized intensity ratio'] = \
            self.track_table.loc[self.track_table.particle == label, 'normalized intensity ratio'] / norm_value


            norm_value_folds = self.track_table.loc[self.track_table.particle == label, 'Number of Folds'].iloc[0]
            self.track_table.loc[self.track_table.particle == label, 'normalized folds'] = \
            self.track_table.loc[self.track_table.particle == label, 'normalized folds'] / norm_value_folds

            norm_value_volume = self.track_table.loc[self.track_table.particle == label, 'volume'].iloc[0]
            self.track_table.loc[self.track_table.particle == label, 'normalized volume'] = \
            self.track_table.loc[self.track_table.particle == label, 'normalized volume'] / norm_value_volume

