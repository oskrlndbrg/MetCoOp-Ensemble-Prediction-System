import pandas as pd
import pyproj
import numpy as np
import netCDF4 # Requires v1.3.1 (pip3 install netCDF4==1.3.1). Newer throws error 68
import time
import datetime
import os

PATH_SAVE = r'C:\Users\oskli230\Box\Repositories\Other\MEPS\Vasakronan'

df_dates = pd.date_range(start='11/09/2021', end='11/10/2022',freq='D')

# Remove certain dates (if neccecary)
#date = df_dates[0] #12/05/2016, start: 07/26/2017
#remove_list = pd.to_datetime(pd.Series(['12/05/2016',
#                                        '12/16/2016']))
#df_dates = df_dates.drop(remove_list)

start, end = 0, 67 # Select forecast horizons
lat, lon = 59.854510, 17.703780

strtdt = datetime.datetime.strptime('2021-03-10', "%Y-%m-%d") 
enddt = datetime.datetime.strptime('2023-03-09', "%Y-%m-%d") 
df_dates = pd.date_range(start=strtdt, end=enddt,freq='D')

'''
# After download, check missing dates (in download folder) and download these
# Otherwise, comment out this part of the code
# From here...
dir_list = os.listdir(os.path.expanduser(PATH_SAVE))
strtdt = datetime.datetime.strptime('20190101T0000.txt', "%Y%m%dT0000.txt") # datetime.datetime.strptime(dir_list[0], "%Y%m%dT0000.txt")
enddt = datetime.datetime.strptime('20201101T0000.txt', "%Y%m%dT0000.txt") # datetime.datetime.strptime(dir_list[-1], "%Y%m%dT0000.txt")

df_dates = pd.date_range(start=strtdt, end=enddt,freq='D')

for filename in dir_list:
    dt_temp = datetime.datetime.strptime(filename, "%Y%m%dT0000.txt")
    if dt_temp in df_dates:
        df_dates = df_dates[df_dates != dt_temp]
# ...to here.
'''

for date in df_dates:
    filename = "https://thredds.met.no/thredds/dodsC/meps25epsarchive/" + str(date.year) + "/" + '{:02d}'.format(date.month) + "/" + '{:02d}'.format(date.day) + "/meps_det_2_5km_" + str(date.year) + '{:02d}'.format(date.month) + '{:02d}'.format(date.day) +"T00Z.nc"

    try:
        strt_time = time.time()

        file = netCDF4.Dataset(filename,"r")
        
        print('\n #############################################')
        print('Initialising downloading of: ' + str(date.date()) + '... \n')
        
        proj = pyproj.Proj(file.variables["projection_lambert"].proj4)
        
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
        
        # Get reference times of forecasts
        ref_time = file.variables["forecast_reference_time"][:]
        ref_time = pd.to_datetime(datetime.datetime.utcfromtimestamp(ref_time+0), utc=True) # .tz_convert('Europe/Stockholm')
        ref_times = [ref_time]*(end-start)
        df_ref_times = pd.DataFrame(data=ref_times,columns=['ref_time'])
    
        # Forecast valid times
        valid_time = file.variables["time"][start:end]
        valid_time = np.vstack([pd.to_datetime(datetime.datetime.utcfromtimestamp(np.asarray(i)+0), utc=True) for i in valid_time]) # .tz_convert('Europe/Stockholm')
        df_valid_time = pd.DataFrame(data=valid_time,columns=['valid_time'])
        
        # Surface pressure
        print('Downloading \Surface pressure/...')
        s_time = time.time()
        ps = file.variables['surface_air_pressure'][start:end,0,Iy,Ix] # in Pa
        df_ps = pd.DataFrame(data=ps, columns=['surface_pressure']).round(2)
        print('...\Surface pressure/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # Temperature at 2 meters
        print('Downloading \Surface Temperature/...')
        s_time = time.time()
        temperatures = file.variables["air_temperature_2m"][start:end,0,Iy,Ix] # in degrees K
        df_temp = pd.DataFrame(data=temperatures, columns=['air_temp_2m']).round(2)
        print('...\Surface Temperature/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # Relative humidity at 2 meters
        print('Downloading \Relative humidity/...')
        s_time = time.time()
        rel = file.variables["relative_humidity_2m"][start:end,0,Iy,Ix] # in %
        df_rel = pd.DataFrame(data=rel, columns=['rel_hum_2m']).round(2)
        print('...\Relative humidity/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # Total cloud cover
        print('Downloading \Total cloud cover/...')
        s_time = time.time()
        tcc = file.variables["cloud_area_fraction"][start:end,0,Iy,Ix] # in %
        df_tcc = pd.DataFrame(data=tcc, columns=['tcc']).round(2)
        print('...\Total cloud cover/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # Total precipitation 
        print('Downloading \Total precipitation/...')
        s_time = time.time()
        prec = file.variables["precipitation_amount_acc"][start:end,0,Iy,Ix] # in kg/m^2
        df_prec = pd.DataFrame(data=prec, columns=['precipitation']).round(2)
        print('...\Total precipitation/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # WS U-component at 10 meters
        print('Downloading \ws u at 10 meters/...')
        s_time = time.time()
        x_10m_wind = file.variables["x_wind_10m"][start:end,0,Iy,Ix] # in m/s
        df_wsx = pd.DataFrame(data=x_10m_wind,columns=['wind_x_10m']).round(2)
        print('...\ws u at 10 meters/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # WS V-component at 10 meters
        print('Downloading \ws v at 10 meters/...')
        s_time = time.time()
        y_10m_wind = file.variables["y_wind_10m"][start:end,0,Iy,Ix] # in m/s
        df_wsy = pd.DataFrame(data=y_10m_wind,columns=['wind_y_10m']).round(2)
        print('...\ws v at 10 meters/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # GHI at surface
        print('Downloading \GHI/...')
        s_time = time.time()
        ghi_accumulated = file.variables['integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time'][:,0,Iy,Ix] / 3600 # units: Ws/m2
        ghi = np.diff(ghi_accumulated,axis=0,prepend=0) # Prepend with 0 because it is the integral of shortwave flux.
        ghi = ghi[start:end] # Select the forecast period.
        df_ghi = pd.DataFrame(data=ghi,columns=['ghi']).round(2)
        print('...\GHI/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        # GHI clear sky at surface
        print('Downloading \TOA GHI/...')
        s_time = time.time()
        toa_accumulated = file.variables['integral_of_toa_downwelling_shortwave_flux_wrt_time'][start:end,0,Iy,Ix] / 3600 # units: Ws/m2
        toa = np.diff(toa_accumulated,axis=0,prepend=0) # Prepend with 0 because it is the integral of shortwave flux.
        df_toa = pd.DataFrame(data=toa,columns=['toa_ghi']).round(2)
        print('...\TOA GHI/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')
                
        # DNI (direct normal irradiance) at surface
        print('Downloading \DNI/...')
        s_time = time.time()
        dni_accumulated = file.variables['integral_of_surface_direct_normal_irradiance_wrt_time'][start:end,0,Iy,Ix] / 3600 # units: Ws/m2
        dni = np.diff(dni_accumulated,axis=0,prepend=0) # Prepend with 0 because it is the integral of shortwave flux.
        df_dni = pd.DataFrame(data=dni,columns=['dni']).round(2)
        print('...\DNI/ was downloaded for: \n ' + str(round(time.time() - s_time, 2)) + ' seconds. \n')

        file.close() # Close server, sometimes throws error otherwise
             
        df_complete = pd.concat([df_ref_times, df_valid_time, df_ps, df_temp, df_rel, df_tcc, df_prec, df_wsx, df_wsy, df_ghi, df_toa, df_dni], axis = 1)
    
        fname = os.path.join(PATH_SAVE, ref_time.strftime("%Y%m%dT%H%M")+".txt")
        
        df_complete.to_csv(fname, sep="\t")
        
        print('The NWP data were saved to: ' + fname) 
              
        print('\n It took:' + str(time.time() - strt_time), ' seconds to download the data. \n')
        
    except (KeyError, OSError, RuntimeError) as err:
        print(err)
        fname = os.path.join(PATH_SAVE, 'Failed', ref_time.strftime("%Y%m%dT%H%M")+".txt")
        # Store the error
        df_complete = pd.DataFrame(err)
        df_complete.to_csv(fname, sep="\t")
        
#%% Compile the downloaded data
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
