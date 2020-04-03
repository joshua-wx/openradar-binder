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
    fig = plt.figure(figsize=(10, 8)) 

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
   
    out_ffn = img_path + '/' + os.path.basename(odim_ffn)[:-3] + '.png' 

    plt.savefig(out_ffn, dpi=100)  # Saving figure. 
    plt.close()  # Release memory 
    fig.clf()  # Clear figure 

    #todo
    #-fix title

def build_animation(cdict): 

    #config
    base_url     = 'http://dapds00.nci.org.au/thredds/fileServer/rq0' #base url for NCI dataset  
    
    #parse inputs 
    radar_id_str = str(cdict['rid']).zfill(2) #convert radar id to a string and fill with a leading 0 if only one digit 
    start_dt     = datetime.strptime(cdict['rdate'] + ' ' + cdict['rtime'], '%Y/%m/%d %H:%M') 
    end_dt       = start_dt + timedelta(hours = cdict['rdur']) 

    #build request filename url 
    zip_fn       = radar_id_str + '_' + start_dt.strftime('%Y%m%d') + '.pvol.zip'
    zip_ffn      = '/tmp/' + zip_fn
    request_url  = '/'.join([base_url, 'odim_pvol', radar_id_str, start_dt.strftime('%Y'), 'vol', zip_fn]) 

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
        file_dt_list.append(datetime.strptime(os.path.basename(fname)[3:18],'%Y%m%d_%H%M%S')) 

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
        
    #zip files and move to local directory
    now = datetime.now()
    img_zip_fn = radar_id_str + '_' + start_dt.strftime('%Y-%m-%d') + '.image_request.' + now.strftime("%d-%m-%Y_%H:%M:%S") + '.zip'
    img_zip_ffn = './images/' + img_zip_fn
    #zip up
    zipf = zipfile.ZipFile(img_zip_ffn, 'w', zipfile.ZIP_DEFLATED)
    files_to_zip = glob(img_path + '/*')
    for item in files_to_zip:
        zipf.write(item)
    zipf.close()
    
    #remove temp images
    os.system('rm -rf ' + img_path)
    
    local_file = FileLink(img_zip_ffn, result_html_prefix="Click here to download: ")
    display(local_file)