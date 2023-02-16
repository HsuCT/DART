PRO gdc_ephem_write_ncdf, data, ncdf_fname, scid, gzip=gzip, output_path=output_path

 ;; keyword processing
 IF  ~KEYWORD_SET(output_path) THEN output_path = './'

 ;; error handling.
  catch, error
  IF (error NE 0) THEN BEGIN
      catch, /cancel
;;      void = error_message()
      obj_destroy, fileobj
      RETURN
  ENDIF

;; check keywords
  IF ~KEYWORD_SET(gzip) THEN gzip = 0
  IF (gzip GT 9) THEN gzip = 9
    
  tdim = N_ELEMENTS(data.time)
  vdim = 3
   
;; initial definition of objects.
  fileobj = obj_new()

;; open a new file for writing.
  print,'Writing file: ', output_path + PATH_SEP() + ncdf_fname

  fileobj = OBJ_NEW("ncdf_file", /create, /clobber, /netcdf4_format, output_path + PATH_SEP() + ncdf_fname)
  IF OBJ_VALID(fileobj) EQ 0 THEN $
      message, 'invalid file object returned from ncdf_file init.'
  
;; add global attributes to the file.
  fileobj -> writeglobalattr, datatype='char', 'title', 'GDC Design Reference Mission (DRM) data.'
  fileobj -> writeglobalattr, datatype='char', 'description', $
      'Ephemeris data for satellite ' + scid + ', covering all science phases from mission day 91 through day 1095.'
  fileobj -> writeglobalattr, 'Created', systime(), datatype='char'
  
;; add dimensions to the file.
  fileobj -> writedim, 'n_timestamps', tdim, object=tdimobj
  fileobj -> writedim, 'vec_comp', vdim, object=vdimobj
  
;; get the dimension names.
  dimnames_t = [tdimobj->getname()]
  dimnames_v = [vdimobj->getname(),tdimobj->getname()]

;; define variable fors the file.
  fileobj -> writevardef, 'time', dimnames_t, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'UTC Time in sec since 1970-01-01/00:00:00. ' $
                           + 'First element is set at a reference time, t0 = ' + TIME_STRING(data.time[0])
  fileobj -> writevardata, 'time', data.time

  fileobj -> writevardef, 'elapsed_time', dimnames_t, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Time in days since mission launch.'
  fileobj -> writevardata, 'elapsed_time', data.elapsed_time
  
  fileobj -> writevardef, 'scpos_j2k', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite position in J200 coordinates, in km.'
  fileobj -> writevardata, 'scpos_j2k', data.scpos_j2k

  fileobj -> writevardef, 'scvel_j2k', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite velocity in J200 coordinates, in km/s.'
  fileobj -> writevardata, 'scvel_j2k', data.scvel_j2k

  fileobj -> writevardef, 'scpos_gei', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite position in GEI coordinates, in km.'
  fileobj -> writevardata, 'scpos_gei', data.scpos_gei

  fileobj -> writevardef, 'scvel_gei', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite velocity in GEI coordinates, in km/s.'
  fileobj -> writevardata, 'scvel_gei', data.scvel_gei
  
  fileobj -> writevardef, 'scpos_geo', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite position in GEO coordinates, in km.'
  fileobj -> writevardata, 'scpos_geo', data.scpos_geo

  fileobj -> writevardef, 'scvel_geo', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite velocity in GEI coordinates, in km/s.'
  fileobj -> writevardata, 'scvel_geo', data.scvel_geo
  
  fileobj -> writevardef, 'scpos_gmag', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite position in GMag coordinates, in km.'
  fileobj -> writevardata, 'scpos_gmag', data.scpos_gmag

  fileobj -> writevardef, 'scvel_gmag', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Satellite velocity in GMag coordinates, in km/s.'
  fileobj -> writevardata, 'scvel_gmag', data.scvel_gmag

  fileobj -> writevardef, 'geod_lat', dimnames_t, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Geodetic latitude, in deg, from -90 --> 90.'
  fileobj -> writevardata, 'geod_lat', data.geo_lat
  
  fileobj -> writevardef, 'geod_lon', dimnames_t, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Geodetic longitude, in deg from 0 --> 360.'
  fileobj -> writevardata, 'geod_lon', data.geo_lon

  fileobj -> writevardef, 'geod_alt', dimnames_t, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Geodetic altitude, in km.'
  fileobj -> writevardata, 'geod_alt', data.geo_alt

  fileobj -> writevardef, 'gmag_lat', dimnames_t, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Geomagnetic latitude, in deg, from -90 --> 90.'
  fileobj -> writevardata, 'gmag_lat', data.gmag_lat
  
  fileobj -> writevardef, 'gmag_lon', dimnames_t, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Geomagnetic longitude, in deg from 0 --> 360.'
  fileobj -> writevardata, 'gmag_lon', data.gmag_lon
  
  fileobj -> writevardef, 'sunpos', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Solar position unit vector in J200 coordinates (dimensionless units).'
  fileobj -> writevardata, 'sunpos', data.sunpos_j2k

  fileobj -> writevardef, 'eqpm_pos', dimnames_v, datatype='float', object=dataobj, gzip=gzip
  fileobj -> writevarattr, dataobj, 'comment', 'Prime meridian at the equator unit vector in J200 coordinates (dimensionless units).'
  fileobj -> writevardata, 'eqpm_pos', data.pmpos_j2k
  
;; ; sync the file by writing memory to disk.
;;   fileobj -> sync
  
; browse the file.
;  fileobj -> browse, xoffset=250, yoffset=250, title='AMPERE'
  
; destroy this file object.
  obj_destroy, fileobj

return
end
