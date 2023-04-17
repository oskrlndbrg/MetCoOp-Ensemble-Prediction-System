'''
This script uses the Hyposometric equation in order to calculate the wind speed at hub height.
The background is that estimates of the hub-height wind speed, e.g., the log-law (https://en.wikipedia.org/wiki/Log_wind_profile)
do not take into account the different wind speed variability profiles at the different height.
'''

import pandas as pd
import pyproj
import numpy as np
import netCDF4 # Requires v1.3.1 (pip3 install netCDF4==1.3.1). Newer throws error 68
import time
import datetime
import os

PATH_SAVE = r'C:\Users\oskli230\Box\Repositories\Other\MEPS\Vasakronan\WS_HH'

start, end = 0, 24 # Select forecast horizons
lat, lon = 59.854510, 17.703780
hub_height = 80 # Hub height to be interpolated to (in meters)

strtdt = datetime.datetime.strptime('2022-01-16', "%Y-%m-%d") 
enddt = datetime.datetime.strptime('2023-03-09', "%Y-%m-%d") 
df_dates = pd.date_range(start=strtdt, end=enddt,freq='D')

for date in df_dates:
    temp_time = time.time()
    filename = "https://thredds.met.no/thredds/dodsC/meps25epsarchive/" + str(date.year) + "/" + '{:02d}'.format(date.month) + "/" + '{:02d}'.format(date.day) + "/meps_det_2_5km_" + str(date.year) + '{:02d}'.format(date.month) + '{:02d}'.format(date.day) +"T00Z.nc"
    
    try:
        strt_time = time.time()

        file = netCDF4.Dataset(filename,"r")
        
        print('\n #############################################')
        print('Initialising downloading of: ' + str(date.date()) + '... \n')
                
        # See information of the projection by: print(ncfile.variables["projection_lambert"])
        crs = pyproj.CRS.from_cf(
            {
                "grid_mapping_name": "lambert_conformal_conic",
                "standard_parallel": [63.3, 63.3],
                "longitude_of_central_meridian": 15.0,
                "latitude_of_projection_origin": 63.3,
                "earth_radius": 6371000.0,
            }
            )
        
        # Transformer to project from ESPG:4368 (WGS:84) to their lambert_conformal_conic
        proj = pyproj.Proj.from_crs(4326, crs, always_xy=True)
        
        # Compute projected coordinates of lat/lon points
        X,Y = proj.transform(lon,lat)
        
        # Find nearest neighbour
        x = file.variables['x'][:]
        y = file.variables['y'][:]
        
        Ix = np.argmin(np.abs(x - X))
        Iy = np.argmin(np.abs(y - Y))       
        
        # Hybrid variables
        ap = file.variables['ap'][:]
        b = file.variables['b'][:]
        
        # Air temperature at different model levels
        ta = file.variables['air_temperature_ml'][start:end,-10:,Iy,Ix]
        
        # Wind speeds at different model levels
        u = file.variables['x_wind_ml'][start:end,-10:,Iy,Ix]
        v = file.variables['y_wind_ml'][start:end,-10:,Iy,Ix]
        ws = (u**2 + v**2)**(1/2) # Note: opposite order (high->low)
        ws = np.flip(ws, axis=1) # Therefore, flip it
        
        # Surface pressure
        ps = file.variables['surface_air_pressure'][start:end,0,Iy,Ix]
        
        # Calculate pressure at each model level and time step
        p = np.repeat(ap,len(list(range(start,end)))).reshape(65,len(list(range(start,end)))) + np.outer(b,ps)
        dp = np.flip(np.diff(p,axis=0),axis=0)
        dp = dp[0:10,:]
        
        R = 287
        rho = np.flip(p[-10:,:].transpose()/R/ta,axis=1)
        
        dz = dp.transpose()/rho/9.81
        z = np.cumsum(dz,axis=1)
        
        # Wind speed at hub height
        z_diff = abs(z-hub_height)
        wshh = [] # Wind speed at hub height
        for i in range(z_diff.shape[0]):
            idx = sorted(range(len(z_diff[i,:])), key = lambda sub: z_diff[i,:][sub])[:2] # Find index of the two closest heights
            wshh = np.append(wshh, np.interp(hub_height, z[i,idx], ws[i,idx])) # Interpolate to the correct hub height wind speed
        
        df_wshh = pd.DataFrame(data=wshh,columns=['ws_' + str(hub_height) + 'm']).round(2)
        
        # WS at 10 meters
        x_10m_wind = file.variables["x_wind_10m"][start:end,0,Iy,Ix]
        y_10m_wind = file.variables["y_wind_10m"][start:end,0,Iy,Ix]
        wind_speed_10m = (x_10m_wind[:]**2 + y_10m_wind[:]**2)**(1/2)
        df_ws = pd.DataFrame(data=wind_speed_10m, columns=['ws_10m']).round(2)
        
        # Get reference times of forecasts
        ref_time = file.variables["forecast_reference_time"][:]
        ref_time = pd.to_datetime(datetime.datetime.utcfromtimestamp(ref_time+0), utc=True) # .tz_convert('Europe/Stockholm')
        ref_times = [ref_time]*(end-start)
        df_ref_times = pd.DataFrame(data=ref_times,columns=['ref_time'])
    
        # Forecast valid times
        valid_time = file.variables["time"][start:end]
        valid_time = np.vstack([pd.to_datetime(datetime.datetime.utcfromtimestamp(np.asarray(i)+0), utc=True) for i in valid_time]) # .tz_convert('Europe/Stockholm')
        df_valid_time = pd.DataFrame(data=valid_time,columns=['valid_time'])
    
        file.close() # Close server, sometimes throws error otherwise
        
        df_complete = pd.concat([df_ref_times, df_valid_time, df_wshh, df_ws], axis = 1)
    
        fname = os.path.join(PATH_SAVE, ref_time.strftime("%Y%m%dT%H%M")+".txt")
        
        df_complete.to_csv(fname, sep="\t")
        
        print('The NWP data were saved to: ' + fname) 
              
        print('\n It took:' + str(time.time() - strt_time), ' seconds to download the data. \n')
        
    except (KeyError, OSError, RuntimeError) as err:
        print(err)
        fname = os.path.join(PATH_SAVE, 'Failed', date.strftime("%Y%m%dT%H%M")+".txt")
        # Store the error
        df_complete = pd.DataFrame()
        df_complete.to_csv(fname, sep="\t")

#%% 
dir_list = os.listdir(os.path.expanduser(PATH_SAVE))

df_nwp = pd.DataFrame()

for filename in dir_list:
    print(filename)
    nwp = pd.read_csv(os.path.join(PATH_SAVE, filename), sep='\t', parse_dates = True).iloc[24:48, 2:]
    df_nwp = df_nwp.append(nwp, ignore_index=True)

df_nwp["valid_time"] = pd.to_datetime(df_nwp["valid_time"])
df_nwp = df_nwp.set_index("valid_time")
df_nwp.index = (df_nwp.index).tz_localize(None)

df_nwp.to_csv(os.path.join(PATH_SAVE,"{}.{}".format("compiled_nwp","txt")),sep="\t")
