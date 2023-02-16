#!/glade/u/home/lien/miniconda3/envs/py39/bin/python
import numpy as np
from scipy.interpolate import griddata
from scipy.io import savemat
import xarray as xr
import os
import matplotlib.pyplot as plt
import pickle
import sys

#%%
def get_sc_lonlat(path_dat):
    # Open the file with read only permit
    f = open(path_dat, 'r')
    lon,lat = [],[]
    while True:
        # read line
        line = f.readline()
        # check if line is not empty
        if not line:
            break
        else:
            #print(line)
            lon.append(float(line[0:13]))
            lat.append(float(line[13:]))

    f.close()

    lon = np.array(lon)
    lat = np.array(lat)

    return lon,lat

def get_sc_t(path_dat):
    # Open the file with read only permit
    f = open(path_dat, 'r')
    t_sc1,t_sc2 = [],[]
    while True:
        # read line
        line = f.readline()
        # check if line is not empty
        if not line:
            break
        else:
            #print(line)
            t_sc1.append(float(line[0:13]))
            t_sc2.append(float(line[13:]))

    f.close()

    t_sc1 = np.array(t_sc1)
    t_sc2 = np.array(t_sc2)

    return t_sc1, t_sc2

def get_xyz_intp(glons, glats, alts):
    r0 = 6378.
    x2_obs, y2_obs, z2_obs = sph2cart(r0 + alts, glats, glons)
    x2_intp = x2_obs.flatten()
    y2_intp = y2_obs.flatten()
    z2_intp = z2_obs.flatten()

    return x2_intp, y2_intp, z2_intp

def sph2cart(r, lat, lon):
    phi = (-lat + 90.) * np.pi/180.
    thta = lon * np.pi/180.

    x = r * np.sin(phi) * np.cos(thta)
    y = r * np.sin(phi) * np.sin(thta)
    z = r * np.cos(phi)
    return x, y, z # = sph2cart(r, lat, lon)


def plt1(glonc_sc, glatc_sc,
    LON_model_p1, LAT_model_p1,
    LON_model_p2, LAT_model_p2,
    LON_model_p3, LAT_model_p3,
    LON_model_p4, LAT_model_p4):
    fig,ax = plt.subplots(figsize=(10.,5.))

    l0, = ax.plot( glonc_sc, glatc_sc, label='s/c', marker='x' )
    l1, = ax.plot( LON_model_p1, LAT_model_p1, label='p1', marker='$o$' )
    l2, = ax.plot( LON_model_p2, LAT_model_p2, label='p2', marker='$o$'  )
    l3, = ax.plot( LON_model_p3, LAT_model_p3, label='p3', marker='$o$'  )
    l4, = ax.plot( LON_model_p4, LAT_model_p4, label='p4', marker='$o$'  )

    ax.legend(handles=[l0, l1, l2, l3, l4])

    ax.set_xlabel('Geo. Longitude', fontsize=15)
    # ax.set_xlim( -180, 180 )
    # ax.set_xticks( np.arange(-135,135+45,45) )

    ax.set_ylabel('Geo. Latitude', fontsize=15)
    # ax.set_ylim( -90, 90 )
    # ax.set_yticks( np.arange(-60,60+3,30) )

    ax.tick_params(axis='both', which='major', labelsize=15)

    ax.grid(which='major', linestyle='-')

    str_title = 'test'
    ax.set_title( str_title, fontsize=15)

    fig.tight_layout()
    fn = 'test.png'
    fig.savefig( fn, dpi=300)
    print('figure save in:')
    print(' ', fn )

    plt.close()
    return

def plt2(c_model_p1, alts_model_p1,
    c_model_p2, alts_model_p2,
    c_model_p3, alts_model_p3,
    c_model_p4, alts_model_p4,
    c_intp, alts_prf):

    fig,ax = plt.subplots(figsize=(10.,5.))
    l1, = ax.plot( c_model_p1, alts_model_p1, 'r', label='p1' )
    l2, = ax.plot( c_model_p2, alts_model_p2, 'g', label='p2' )
    l3, = ax.plot( c_model_p3, alts_model_p3, 'b', label='p3' )
    l4, = ax.plot( c_model_p4, alts_model_p4, 'm', label='p4' )
    l0, = ax.plot( c_intp,          alts_prf, 'k', label='s/c', marker='o', fillstyle='none' )

    ax.legend(handles=[l1, l2, l3, l4, l0])

    ax.set_xlabel(var_name_intp + ' [' + unitvar_name_intp + ']', fontsize=15)
    # ax.set_xlim( -180, 180 )
    # ax.set_xticks( np.arange(-135,135+45,45) )

    ax.set_ylabel('Altitude [km], Zgm = Z3/(1-Z3/R_earth)', fontsize=15)
    # ax.set_ylim( -90, 90 )
    # ax.set_yticks( np.arange(-60,60+3,30) )

    ax.tick_params(axis='both', which='major', labelsize=15)

    ax.grid(which='major', linestyle='-')

    str_title = 'test'
    ax.set_title( str_title, fontsize=15)

    fig.tight_layout()
    fn = 'test2.png'
    fig.savefig( fn, dpi=300)
    print('figure save in:')
    print(' ', fn )

    plt.close()
    return

# def get_wx_lonlat(path_wx):

#     return lon,lat
#%% arginp
# inp_varname = 'U' # sys.argv[1] # O
# inp_unitvarname = 'm/s' # sys.argv[2] # O
# inp_alt = 150. #float( sys.argv[3]) # 150.
# inp_sc = 'sc1' #sys.argv[4] # sc1
# idx_t_sc_start = 0 #int( sys.argv[5] )
# idx_t_sc_end = 10 #int( sys.argv[6] )

inp_varname = sys.argv[1] # O
inp_unitvarname = sys.argv[2] # O
inp_alt = float( sys.argv[3] ) # 150.
inp_sc = sys.argv[4] # sc1
idx_t_sc_start = int( sys.argv[5] )
idx_t_sc_end = int( sys.argv[6] )
ifplt = int( sys.argv[7] )

idc_t_sc = range(idx_t_sc_start,idx_t_sc_end)
nt = len( idc_t_sc )

# 
# example:
# /glade/u/home/lien/miniconda3/envs/py39/bin/python sample_waccmx_jobs.py U m/s 150. sc1 0 10 0

#%%
#dir_wx_c = '/glade/campaign/univ/ucub0111/liuh/fx2000_ne120pg3L273.001br/atm/hist'
dir_wx_c = '/glade/campaign/univ/ucub0111/liuh/fx2000_ne120pg3L273.001_1min'

fn_wx_p1 = 'fx2000_ne120pg3L273.001_1min.cam.h1.0001-01-30-'
fn_wx_p3 = '.nc'

#tp_height_all = np.arange(100., 150.+10., 10.) # I assume this is altitude in km shown in plot_orbit_sim.pro
tp_height_all = np.arange(inp_alt, inp_alt+10., 10.) # I assume this is altitude in km shown in plot_orbit_sim.pro
ntp_height_all = tp_height_all.shape[0]

var_name_intp_all = [inp_varname]
unitvar_name_intp_all = [inp_unitvarname]
nvar_name_intp_all = len(var_name_intp_all)

#var_unit_intp_all = ['mol/mol']
#nvar_unit_intp_all = len(var_unit_intp_all)
sc_name = inp_sc

#waccmx_ts = np.arange(300., 86100. + 300., 300.)
waccmx_ts = np.arange(0., 21600. + 60., 60.)

alts_prf = np.arange(100.,250.,2.)
nalts_prf = alts_prf.shape[0]
#for idx_tp_height in range(0, ntp_height_all):
for idx_tp_height in range(0, 1):
    tp_height = tp_height_all[idx_tp_height]
    ## load txt for sc locations
    path_dat_sc_loc = '{1:s}_alt_{0:5.1f}km.dat'.format( tp_height, sc_name )
    print('loading '+path_dat_sc_loc)
    glons_sc,glats_sc = get_sc_lonlat( path_dat_sc_loc )
    #alts_sc = np.ones(glons_sc.shape) * tp_height

    # ## load txt for sc time in second
    # path_dat_sc12_t = 'tp12_{1:s}_t_{0:5.1f}km.dat'.format( tp_height, sc_name[-4:] )
    # print('loading '+path_dat_sc12_t)
    # ts_sc1, ts_sc2 = get_sc_t( path_dat_sc12_t )

    # if sc_name[0:3]=='tp1':
    #     ts_sc = ts_sc1
    # elif sc_name[0:3]=='tp2':
    #     ts_sc = ts_sc2

    ts_sc = np.arange(0., 86400., 60.) # all the same
    nts_sc = ts_sc.shape[0]


    var_name_intp = inp_varname
    unitvar_name_intp = inp_unitvarname
    #var_unit_intp = var_unit_intp_all[idx_var_name_intp]
    #print(var_name_intp)


    if 'vars_sampled' in locals():
        del vars_sampled
    vars_sampled = {}
    var_sampled = np.full( (nt,nalts_prf), np.nan, dtype=np.double )
    # ## interpolating for each satellite location & saving
    for count, idx_sc_t in enumerate(idc_t_sc):
    #for idx_sc_t in range( 100, 101 ):
        tc_sc = ts_sc[ idx_sc_t ]

        #finding nearest time
        idx_nearest_wacx_t = np.argmin( np.abs( tc_sc - waccmx_ts) )
        waccmx_tc = waccmx_ts[idx_nearest_wacx_t]
        #print( 'S/C time: {0:e} sec of a day, WACCMX nearest time: {1:e} sec of a day'.format(tc_sc, waccmx_tc ) )
        
        #determining time format in file names
        fn_wx_p2 = '{0:05d}'.format( int(waccmx_tc) )
        #print( fn_wx_p2 )

        fn_wx_c = fn_wx_p1 + fn_wx_p2 + fn_wx_p3
        path_wx_c = os.path.join( dir_wx_c, fn_wx_c)
        print('S/C time {1:5d}{0:5d}, loading '.format( int(idx_sc_t),int(tc_sc)) + path_wx_c)
        ds = xr.open_dataset(path_wx_c)

        glon_wx = ds['lon'].values
        # glon_wx[np.where(glon_wx>180.)] = glon_wx[np.where(glon_wx>180.)] - 360.
        # glon_wx = np.concatenate( (glon_wx[721:],glon_wx[0:721]), axis=0 )

        glat_wx = ds['lat'].values
        lev_wx = ds['lev'].values

        LAT_model, LON_model = np.meshgrid( glat_wx, glon_wx, indexing='ij')

        glonc_sc = glons_sc[idx_sc_t]
        #glonc_sc = 0
        if glonc_sc<0:
            glonc_sc = glonc_sc + 360.
        glatc_sc = glats_sc[idx_sc_t]
        #glatc_sc = 89
        
        idx_glon_r_raw = np.searchsorted(glon_wx, glonc_sc, side='right')
        idx_glon_l_raw = np.searchsorted(glon_wx, glonc_sc, side='right')-1

        idx_glon_r = np.remainder( idx_glon_r_raw, glon_wx.shape[0] )
        idx_glon_l = np.remainder( idx_glon_l_raw, glon_wx.shape[0] )

        idx_glat_r = np.searchsorted(glat_wx, glatc_sc, side='right')
        idx_glat_l = np.searchsorted(glat_wx, glatc_sc, side='right')-1

        LON_model_p1 = LON_model[idx_glat_r, idx_glon_r]
        LAT_model_p1 = LAT_model[idx_glat_r, idx_glon_r]

        LON_model_p2 = LON_model[idx_glat_r, idx_glon_l]
        LAT_model_p2 = LAT_model[idx_glat_r, idx_glon_l]

        LON_model_p3 = LON_model[idx_glat_l, idx_glon_r]
        LAT_model_p3 = LAT_model[idx_glat_l, idx_glon_r]

        LON_model_p4 = LON_model[idx_glat_l, idx_glon_l]
        LAT_model_p4 = LAT_model[idx_glat_l, idx_glon_l]

        # print( 'sc: glon = {0:12.8f}, glat = {1:12.8f}'.format( glonc_sc, glatc_sc ) )
        # print( 'p1: glon = {0:12.8f}, glat = {1:12.8f}'.format( LON_model_p1, LAT_model_p1 ) )
        # print( 'p2: glon = {0:12.8f}, glat = {1:12.8f}'.format( LON_model_p2, LAT_model_p2 ) )
        # print( 'p3: glon = {0:12.8f}, glat = {1:12.8f}'.format( LON_model_p3, LAT_model_p3 ) )
        # print( 'p4: glon = {0:12.8f}, glat = {1:12.8f}'.format( LON_model_p4, LAT_model_p4 ) )

        if ifplt==1:
            plt1(glonc_sc, glatc_sc,
                LON_model_p1, LAT_model_p1,
                LON_model_p2, LAT_model_p2,
                LON_model_p3, LAT_model_p3,
                LON_model_p4, LAT_model_p4)


        Z3 = ds['Z3'].values/1.e3
    #     #print(Z3.shape)
        R_earth = 6378.
        Zgm = Z3/(1-Z3/R_earth)

        c_model = ds[var_name_intp].values

        alts_model_p1 = np.squeeze( Zgm[0,:,idx_glat_r, idx_glon_r] )
        alts_model_p2 = np.squeeze( Zgm[0,:,idx_glat_r, idx_glon_l] )
        alts_model_p3 = np.squeeze( Zgm[0,:,idx_glat_l, idx_glon_r] )
        alts_model_p4 = np.squeeze( Zgm[0,:,idx_glat_l, idx_glon_l] )
        alts_model = np.concatenate( (alts_model_p1, 
                                      alts_model_p2,
                                      alts_model_p3,
                                      alts_model_p4), axis=0)


        c_model_p1 = np.squeeze( c_model[0,:,idx_glat_r, idx_glon_r] )
        c_model_p2 = np.squeeze( c_model[0,:,idx_glat_r, idx_glon_l] )
        c_model_p3 = np.squeeze( c_model[0,:,idx_glat_l, idx_glon_r] )
        c_model_p4 = np.squeeze( c_model[0,:,idx_glat_l, idx_glon_l] )
        c_model = np.concatenate( (c_model_p1, 
                                   c_model_p2,
                                   c_model_p3,
                                   c_model_p4), axis=0)


        glons_model_p1 = np.ones( c_model_p1.shape ) * LON_model_p1
        glons_model_p2 = np.ones( c_model_p2.shape ) * LON_model_p2
        glons_model_p3 = np.ones( c_model_p3.shape ) * LON_model_p3
        glons_model_p4 = np.ones( c_model_p4.shape ) * LON_model_p4
        # glons_model = np.concatenate( (glons_model_p1,
        #                                glons_model_p2,
        #                                glons_model_p3,
        #                                glons_model_p4), axis=0)        
        glons_model = np.array( [LON_model_p1,
                                 LON_model_p2,
                                 LON_model_p3,
                                 LON_model_p4] )

        glats_model_p1 = np.ones( c_model_p1.shape ) * LAT_model_p1
        glats_model_p2 = np.ones( c_model_p2.shape ) * LAT_model_p2
        glats_model_p3 = np.ones( c_model_p3.shape ) * LAT_model_p3
        glats_model_p4 = np.ones( c_model_p4.shape ) * LAT_model_p4
        # glats_model = np.concatenate( (glats_model_p1,
        #                                glats_model_p2,
        #                                glats_model_p3,
        #                                glats_model_p4), axis=0)
        glats_model = np.array( [LAT_model_p1,
                                 LAT_model_p2,
                                 LAT_model_p3,
                                 LAT_model_p4] )
        
        # c_intp = griddata( (glons_model, glats_model, alts_model), c_model,
        #                    (glonc_sc * np.ones( alts_prf.shape ), 
        #                     glatc_sc * np.ones( alts_prf.shape ), 
        #                     alts_prf), 
        #                    method='nearest')
        # step 1 interpolations @ the same altitude
        c_intp_p1 = griddata( (alts_model_p1), c_model_p1,
                              (alts_prf), 
                              method='linear')
        c_intp_p2 = griddata( (alts_model_p2), c_model_p2,
                              (alts_prf), 
                              method='linear')
        c_intp_p3 = griddata( (alts_model_p3), c_model_p3,
                              (alts_prf), 
                              method='linear')
        c_intp_p4 = griddata( (alts_model_p4), c_model_p4,
                              (alts_prf), 
                              method='linear')

        c_intp_step2 = np.full( (nalts_prf), np.nan, dtype=np.double )
        # step 2 interpolations @ the lon, lat plane
        for i_alt in range( 0, len(alts_prf) ):
            c_model = np.array( [c_intp_p1[i_alt],
                                 c_intp_p2[i_alt],
                                 c_intp_p3[i_alt],
                                 c_intp_p4[i_alt]] )
            c_intp_step2[i_alt] = griddata( (glons_model, glats_model), c_model,
                                            (glonc_sc, glatc_sc ), 
                                             method='linear')

        # if ifplt==1:
        #     plt2(c_model_p1, alts_model_p1,
        #         c_model_p2, alts_model_p2,
        #         c_model_p3, alts_model_p3,
        #         c_model_p4, alts_model_p4,
        #         c_intp, alts_prf)

        var_sampled[count, :] = c_intp_step2


    # # open a file, where you ant to store the data
    fnpkl = '{0:s}_{1:5.1f}km_{2:s}_tidx{3:05d}_{4:05d}.pkl'.format( var_name_intp,  inp_alt, inp_sc, idx_t_sc_start, idx_t_sc_end )
    filepkl = open( fnpkl, 'wb')
    vars_sampled['glons_sc'] = glons_sc[idx_t_sc_start:idx_t_sc_end]
    vars_sampled['glats_sc'] = glats_sc[idx_t_sc_start:idx_t_sc_end]
    vars_sampled['alts_sc']  = tp_height
    vars_sampled['alts_prf'] = alts_prf
    vars_sampled['ts_sc']    = ts_sc[idx_t_sc_start:idx_t_sc_end]
    vars_sampled[var_name_intp] = var_sampled

    # #saved
    # # dump information to that file
    pickle.dump(vars_sampled, filepkl)

    #close the file
    filepkl.close()


    fnpkl = '{0:s}_{1:5.1f}km_{2:s}_tidx{3:05d}_{4:05d}.pkl'.format( var_name_intp,  inp_alt, inp_sc, idx_t_sc_start, idx_t_sc_end )
    filepkl = open( fnpkl, 'rb')
    b = pickle.load(filepkl)

    print('pwd = ' + os.getcwd())
    print(fnpkl + ' saved!')
