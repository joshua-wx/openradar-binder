#########################################################################

import os
from glob import glob
from datetime import datetime, timedelta
import tempfile #used to create temporary folders to store data
import zipfile #used to extract tar files

import urllib
import pandas
from matplotlib import pyplot as plt
import cartopy
import cartopy.crs as ccrs # A toolkit for map projections
import numpy as np
from tqdm import tqdm
from adjustText import adjust_text
from IPython.display import display, FileLink

import pyart

import warnings
warnings.filterwarnings('ignore')  
warnings.simplefilter('ignore')

#########################################################################

def get_date_from_filename(filename, delimiter='_', date_fmt='%Y%m%d_%H%M%S'):
    """
    INPUT:
    filename (str):
        can contain path, has no impact on this function.
        filename must use the convention IDDDD_DATE. delimiter + everything else 
        DATE must have same format as date
    OUTPUT:
    radar_id (int)
    """
    if not isinstance(filename, str):
        raise ValueError(f"get_id_from_filename: filename is not a string: {filename}")
        return None
    if delimiter not in filename:
        raise ValueError(f"get_id_from_filename: Delimiter not found in filename: {filename}")
        return None
    fn = os.path.basename(filename)
    fn_parts = fn.split(delimiter)
    try:
        dtstr = fn_parts[1] + '_' + fn_parts[2].split('.')[0]
        dt = datetime.strptime(dtstr, date_fmt)
        return dt
    except:
        raise ValueError(f"get_id_from_filename: Failed to extract radar if from: {filename}")
        return None
        
def get_id_from_filename(filename, delimiter='_'):
    """
    INPUT:
    filename (str):
        can contain path, has no impact on this function.
        filename must use the convention IDDDD + delimiter + everything else 
    OUTPUT:
    radar_id (int)
    """
    if not isinstance(filename, str):
        raise ValueError(f"get_id_from_filename: filename is not a string: {filename}")
        return None
    if delimiter not in filename:
        raise ValueError(f"get_id_from_filename: Delimiter not found in filename: {filename}")
        return None
    fn = os.path.basename(filename)
    fn_parts = fn.split(delimiter)
    try:
        radar_id = int(fn_parts[0])
        return radar_id
    except:
        raise ValueError(f"get_id_from_filename: Failed to extract radar if from: {filename}")
        return None

def make_gif(files, output, delay=100, repeat=True,**kwargs):
    """
    Uses imageMagick to produce an animated .gif from a list of
    picture files.
    INPUT:
        files: list of full file paths
        output: full filename for output gif
        delay: delay in ms between animation frames
        repeat: Set to infinite loop
    OUTPUT:
        None
    """
    loop = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s'
              %(delay,loop," ".join(files),output))  

def _read_csv(csv_ffn, header_line):
    """
    CSV reader used for the radar locations file (comma delimited)
    """
    df = pandas.read_csv(csv_ffn, header=header_line)
    as_dict = df.to_dict(orient='list')
    return as_dict

# Function which generates a plot for each sweep 
def _fast_plot(odim_ffn, cdict, img_path): 

    #open figure 
    fig = plt.figure(figsize=(10, 8), facecolor='w') 

    #load radar object 
    my_radar = pyart.aux_io.read_odim_h5(odim_ffn)

    #find limits 
    radar_lat = my_radar.latitude['data'][0]
    radar_lon = my_radar.longitude['data'][0]
    
    plot_range = cdict['prng'] * 1000 #convert to m
    min_lon, min_lat = pyart.core.cartesian_to_geographic_aeqd(-plot_range, -plot_range, radar_lon, radar_lat)
    max_lon, max_lat = pyart.core.cartesian_to_geographic_aeqd(plot_range, plot_range, radar_lon, radar_lat)

    # Set up the GIS projection 
    projection = ccrs.Mercator( 
                    central_longitude=radar_lon, 
                    min_latitude=min_lat, max_latitude=max_lat) 

    #load diplay class for radar using cartopy 
    display = pyart.graph.RadarMapDisplay(my_radar) 

    #set plotting options
    if cdict['field']=='reflectivity':
        vmin = cdict['vmin_ref']
        vmax = cdict['vmax_ref']
        cmap = cdict['cmap_ref']
    elif cdict['field']=='velocity':
        vmin = cdict['vmin_vel']
        vmax = cdict['vmax_vel']
        cmap = cdict['cmap_vel']
    else:
        vmin = cdict['vmin_pol']
        vmax = cdict['vmax_pol']
        cmap = cdict['cmap_pol']
    
    #create plot 
    display.plot_ppi_map(cdict['field'], cdict['tilt'], 
                           projection=projection, colorbar_flag=True,
                           min_lon=min_lon, max_lon=max_lon, min_lat=min_lat, max_lat=max_lat, 
                           vmin=vmin, vmax=vmax, cmap=cmap, 
                           resolution='10m', mask_outside=True)
    #add city markers
    ax = plt.gca()
    fname = cartopy.io.shapereader.natural_earth(resolution='10m', category='cultural', name='populated_places')
    reader = cartopy.io.shapereader.Reader(fname)
    city_list = list(reader.records())
    texts_list = []
    for city in city_list:
        if (((city.attributes['LATITUDE'] >= min_lat) and (city.attributes['LATITUDE'] <= max_lat))
            and ((city.attributes['LONGITUDE'] >= min_lon) and (city.attributes['LONGITUDE'] <= max_lon))
            and (city.attributes['NAME'] not in cdict['hide_place_names'])):
            ax.scatter(city.attributes['LONGITUDE'], city.attributes['LATITUDE'], s=6, color='black',
                       transform=ccrs.PlateCarree(), zorder=5)
            texts_list.append(ax.text(city.attributes['LONGITUDE']+0.01, city.attributes['LATITUDE']+0.01, 
                    city.attributes['NAME'], fontsize=10, transform=ccrs.PlateCarree()))
    #optimise the location of text
    adjust_text(texts_list)

    #add range rings
    display.plot_range_rings(np.array(cdict['range_ring_list'],dtype=np.float),
                             ax=ax, col='k', ls='-', lw=0.5)
            
    #Now we add lat lon lines 
    gl = display.ax.gridlines(draw_labels=True, 
                             linewidth=1, color='gray', alpha=0.5, 
                             linestyle='--') 
    gl.xlabel_style = {'size': 12} 
    gl.ylabel_style = {'size': 12} 
    gl.xlabels_top = False 
    gl.ylabels_right = False 

    #fix title
    current_title = ax.get_title()
    new_title = cdict['nice_radar_name'] + ' ' + current_title
    ax.set_title(new_title)
   
    out_ffn = f'{img_path}/{os.path.basename(odim_ffn)[:-3]}_{cdict["field"]}.png' 

    plt.savefig(out_ffn, dpi=100)  # Saving figure. 
    plt.close()  # Release memory 
    fig.clf()  # Clear figure 

    #todo
    #-fix title

def build_animation(cdict): 

    #config
    base_url     = 'http://dapds00.nci.org.au/thredds/fileServer/rq0' #base url for NCI dataset  
    
    #parse inputs 
    radar_id     = cdict['rid']
    start_dt     = datetime.strptime(cdict['rdate'] + ' ' + cdict['rtime'], '%Y/%m/%d %H:%M') 
    end_dt       = start_dt + timedelta(hours = cdict['rdur']) 

    #build request filename url 
    zip_fn       = str(radar_id) + '_' + start_dt.strftime('%Y%m%d') + '.pvol.zip'
    zip_ffn      = '/tmp/' + zip_fn
    request_url  = '/'.join([base_url, str(radar_id), start_dt.strftime('%Y'), 'vol', zip_fn]) 

    #download the zip file 
    if not os.path.isfile(zip_ffn): 
        print('Fetching:', request_url) 
        urllib.request.urlretrieve(request_url, zip_ffn) 
    else:
        print('File already downloaded:', request_url) 
    
    #extract the zip file to a temporary directory 
    temp_dir = tempfile.mkdtemp() 
    zip_fh = zipfile.ZipFile(zip_ffn) 
    zip_fh.extractall(path = temp_dir) 
    zip_fh.close() 

    #list all the volumes extracted from the zip file 
    file_list = sorted(glob(temp_dir + '/*')) 

    #now let's read the datetime numbers of all the volumes for comparision 
    file_dt_list = [] 
    for i, fname in enumerate(file_list): 
        file_dt_list.append(get_date_from_filename(fname))

    #find the index of volumes within our start and end times 
    file_dt_array    = np.array(file_dt_list) 
    filter_file_list = [] 
    for i, file_dt in enumerate(file_dt_array): 
        if file_dt >= start_dt and file_dt <= end_dt: 
            filter_file_list.append(file_list[i]) 

    #generate images
    print('Generating Images')
    img_path = tempfile.mkdtemp()
    n_img    = len(filter_file_list)
    for odim_ffn in tqdm(filter_file_list, total=n_img): 
        _fast_plot(odim_ffn, cdict, img_path) 

    #create image folder if requires
    if not os.path.exists('images'):
        os.mkdir('images')
    
    #create animation
    now = datetime.now()
    gif_ffn = './images/' + f'{radar_id}_{start_dt.strftime("%Y-%m-%d")}_{cdict["field"]}.animation.{now.strftime("%d-%m-%Y_%H:%M:%S")}.gif'
    files_to_animate = sorted(glob(img_path + '/*'))
    make_gif(files_to_animate, gif_ffn, delay=100, repeat=True)
    
    #zip files and move to local directory
    now = datetime.now()
    img_zip_fn = f'{radar_id}_{start_dt.strftime("%Y-%m-%d")}_{cdict["field"]}.image_request.{now.strftime("%d-%m-%Y_%H:%M:%S")}.zip'
    img_zip_ffn = './images/' + img_zip_fn
    #zip up
    zipf = zipfile.ZipFile(img_zip_ffn, 'w', zipfile.ZIP_DEFLATED)
    files_to_zip = glob(img_path + '/*')
    for item in files_to_zip:
        zipf.write(item)
    zipf.close()
    
    #remove temp images
    os.system('rm -rf ' + img_path)
    
    local_file = FileLink(img_zip_ffn, result_html_prefix="Click here to download zip of images: ")
    display(local_file)
    
    local_file = FileLink(gif_ffn, result_html_prefix="Click here to download animation: ")
    display(local_file)