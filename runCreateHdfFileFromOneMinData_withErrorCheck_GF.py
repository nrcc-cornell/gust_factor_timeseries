'''Run automated error checks on one-minute data.

	Run automated QC on one-minute data. Primarily checking for physically impossible values.
	Additional extreme values that make it through this check are validated manually in subsequest script.

	- introduce automated QC checks on wind data in one-minute HDF5 files
	- write out new version of one-minute HDF5 files containing data with automated checks
	- write out minutes with extreme data that needs manual checking

	bnb2
'''
import os,sys
import h5py
import datetime
import requests
import numpy as np

try :
	import json
except ImportError :
	import simplejson as json

##############################
# list of states

stateList = ['PA']
# Project data directory
projectDataDir = '../data/'
# data resolution directory name
dataResDir = 'asos-onemin'
miss = '-99'
completedStns = []

def acis_ws(method,params) :
	# function to make ACIS web services call
	if method=='StnData':
		acis_url = 'https://data.nrcc.rcc-acis.org/'
	elif method=='StnMeta':
		acis_url = 'https://data.rcc-acis.org/'
	elif method=='GridData':
		acis_url = 'https://grid2.rcc-acis.org/'
	else:
		acis=url = ''
	r = requests.post(acis_url+method,json=params,headers={'Accept':'application/json'})
	try:
		return r.json()
	except:
		return None

# ACIS requests will occur through today
obsEndDate = datetime.datetime.now().strftime('%Y-%m-%d')

# File containing data to check manually
ofile_manual_dates = open('manual_data_dates_test.txt','a')
ofile_manual_check = open('manual_data_check_test.txt','a')

for state in stateList:
	# data source directory
	inDirBase = projectDataDir+dataResDir+'-hdf-files/'+state+'/'
	# HDF output directory
	outDirBase = projectDataDir+dataResDir+'-hdf-files-with-error-check/'+state+'/'
	# create output directory if it doesn't exist
	if not os.path.exists(outDirBase):
		res = os.system('mkdir -p '+outDirBase)
	# get station list for state
	fileList = os.listdir(inDirBase)
	stnList = [ f.split('-')[0] for f in fileList ]
	if len(stnList)==0:
		print ('No stations available for ',state)
		print ('END: ',state)
		continue
	for stn in stnList:
		if stn in completedStns: continue
		print ('START: ',stn)
		# untar src file
		cmd = 'tar xzf '+inDirBase + stn + '-'+dataResDir+'.h5.tgz --directory '+inDirBase
		res = os.system(cmd)
		# read src file
		h5in = h5py.File(inDirBase + stn + '-'+dataResDir+'.h5','r')
		dates = np.array(h5in['dates'])

		winddir_np = np.array(h5in['peakdir'])
		wind_np = np.array(h5in['peakspd'])
		#h5in.close()
		# delete untarred src file (tgz file remains)
		os.remove(inDirBase + stn + '-'+dataResDir+'.h5')

		################################
		### WIND DIRECTION
		################################
		print ('Peak wind direction missing minutes before QC : ', np.count_nonzero(winddir_np == float(miss)))
		winddir_mask = np.ma.mask_or(winddir_np<0, winddir_np>360)
		winddir_ma = np.ma.array(winddir_np, mask=winddir_mask)
		winddir_np = winddir_ma.filled(fill_value=float(miss))
		print ('Peak wind direction missing minutes after QC : ', np.count_nonzero(winddir_np == float(miss)))

		################################
		### WIND SPEED
		################################

		print ('Peak wind speed missing minutes before QC : ', np.count_nonzero(wind_np == float(miss)))
		# Get peak wind speed observations for QC of one-minute data
		# 1. Observed hourly peak wind speed data from ACIS.
		# 2. Calculate max peak wind speed each day.
		# 3. Assign daily max peak wind to each minute of that day. Use this array to compare to minute data.
		# NOTE: ACIS docs describe units as Knots, but they appear to actually be MPH.

		pData = {"sid":stn,"sdate":"2000-01-01","edate":obsEndDate,"elems":[{"vX":42,"vN":1}]}
		rData = acis_ws('StnData',pData)
		if rData is None:
			print ('ERROR problem getting data from ACIS for ',stn)
			continue
		if 'data' not in rData:
			print ('ERROR problem getting data from ACIS for ',stn)
			continue
		windspd_max24hr = []
		print ('ACIS wind dates : ', rData['data'][0][0], rData['data'][-1][0])
		for item in rData['data']:
			obsDataHourly = [v for v in item[1] if v!='M']
			obsData = int(max(obsDataHourly)) if len(obsDataHourly)!=0 else int(miss)
			windspd_max24hr += [obsData]*24*60
		if len(windspd_max24hr)==0: print ('WARNING: no ACIS wind observations for ',stn)
		windspd_max24hr_np = np.array(windspd_max24hr[:wind_np.shape[0]])

		### compare peakspd to ACIS daily wind max gust at this station.
		# check if one-minute peak speeds exceed reported peak speed for the day (+5).
		print (wind_np.shape, windspd_max24hr_np.shape)
		wind_mask_1 = np.logical_and(wind_np > windspd_max24hr_np+5, windspd_max24hr_np!=int(miss))

		# peak winds only reported if speeds > 30 knots. Therefore, test for one-minute
		# speeds greater than 35, when obs values do not include a peak speed reported.
		# In these instances, set minutes to missing.
		wind_mask_2 = np.logical_and(wind_np > 35, windspd_max24hr_np==int(miss))

		# wind speeds must be non-negative
		wind_mask_3 = wind_np < 0

		# wind speeds cannot exceed Mt Washington Record
		wind_mask_4 = wind_np >250

		## wind speeds greater than 100 cannot occur on a day when more than 20 hours observations are missing.

		wind_mask_5 = np.logical_and(wind_np > 100,len(obsDataHourly)<4)


		# combine wind masks to get full mask
		wind_mask = wind_mask_1 | wind_mask_2 | wind_mask_3 | wind_mask_4 | wind_mask_5


		# The above tests should catch erroneous extreme values. Just in case, manually check any remaining values above 80.
		# The following finds minutes that are not masked and have peak wind speed > 80 kts.
		manual_check_minutes_mask = np.logical_and(wind_np > 80, np.logical_not(wind_mask))
		datesToCheck_ma = np.ma.array(dates, mask=np.logical_not(manual_check_minutes_mask))
		datesToCheck = datesToCheck_ma.compressed()
		if len(datesToCheck)!=0:
			#print '*** MANUAL CHECK NEEDED, WIND : '+stn

			for d in datesToCheck:
				idx = list(dates).index(d)
				# file containing dates to check
				stringToOutput = ','.join([state,stn,'wind',str(dates[idx]),str(winddir_np[idx]),str(wind_np[idx]),'obs',str(windspd_max24hr_np[idx])])
				ofile_manual_dates.write(stringToOutput+os.linesep)
				# file containing dates to check, with surrounding days for verification
				for idx_rel in range(-6,7):
					stringToOutput = ','.join([state,stn,'wind',str(dates[idx+idx_rel]),str(winddir_np[idx+idx_rel]),str(wind_np[idx+idx_rel]),'obs',str(windspd_max24hr_np[idx+idx_rel]),str(idx_rel)])
					ofile_manual_check.write(stringToOutput+os.linesep)
			print ('***************************')

		wind_ma = np.ma.array(wind_np, mask=wind_mask)
		wind_np = wind_ma.filled(fill_value=float(miss))
		print ('Peak wind speed missing minutes after QC : ', np.count_nonzero(wind_np == float(miss)))

		pData = {"sids":stn}
		rData = acis_ws('StnMeta',pData)
		if rData is None:
			print ('ERROR problem getting data from ACIS for ',stn)
			continue
		if 'meta' not in rData:
			print ('ERROR problem getting data from ACIS for ',stn)
			continue
		loc = ",".join([str(item) for item in rData['meta'][0]['ll']])
		stn_meta = rData['meta'][0]
		
		##########################
		### create HDF5 file
		##########################
		h5file_out = h5py.File(outDirBase + stn + '-'+dataResDir+'.h5','w')

		dgrp = h5file_out.create_group('asos_one_min')

		dgrp.attrs['ID'] = stn
		dgrp.attrs['name'] = str(stn_meta['name'])
		dgrp.attrs['state'] = state
		dgrp.attrs['longitude'] = str(stn_meta['ll'][0])
		dgrp.attrs['latitude'] = str(stn_meta['ll'][1])
		dgrp.attrs['elevation'] = str(stn_meta['elev'])
		dgrp.attrs['src_wind'] = 'DSI-6405'
		dgrp.attrs['src_url'] = 'Raw text datasets from NCEP, https://www.ncei.noaa.gov/pub/data/asos-onemin/'
		dgrp.attrs['desc'] = 'ASOS one-minute peak wind speed/direction and precipitation. QC completed.'
		dgrp.attrs['created'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

		### dates
		data = list(dates)
		shapeOfHdf5Dataset = (len(data),)
		dset_dates = dgrp.create_dataset('dates', shapeOfHdf5Dataset, dtype='S12')
		dset_dates.attrs['desc'] = 'Date'
		dset_dates.attrs['units'] = 'minutes'
		dset_dates.attrs['format'] = 'YYYYMMDDHHMM'
		dset_dates.attrs['time_zone'] = 'LST'
		dset_dates[:] = data
		del data

		### peak wind direction for max 5-sec average over 1 minute
		data = list(winddir_np)
		shapeOfHdf5Dataset = (len(data),)
		dset_pdir = dgrp.create_dataset('peakdir', shapeOfHdf5Dataset, dtype='<f8')
		dset_pdir.attrs['desc'] = 'Peak Wind Direction'
		dset_pdir.attrs['units'] = 'degrees from North'
		dset_pdir.attrs['missing'] = '-99'
		dset_pdir.attrs['quality'] = 'post_QC'
		dset_pdir[:] = data
		del data

		### peak wind speed for max 5-sec average over 1 minute
		data = list(wind_np)
		shapeOfHdf5Dataset = (len(data),)
		dset_pspd = dgrp.create_dataset('peakspd', shapeOfHdf5Dataset, dtype='<f8')
		dset_pspd.attrs['desc'] = 'Peak Wind Speed'
		dset_pspd.attrs['units'] = 'knots'
		dset_pspd.attrs['missing'] = '-99'
		dset_pspd.attrs['quality'] = 'post_QC'
		dset_pspd[:] = data
		del data

		h5file_out.close()

		# tar and zip output file
		srcPath = os.getcwd()
		os.chdir(outDirBase)
		cmd = 'tar czf '+stn+'-'+dataResDir+'.h5.tgz '+stn+'-'+dataResDir+'.h5'
		res = os.system(cmd)
		os.remove(stn+'-'+dataResDir+'.h5')
		os.chdir(srcPath)

		print ('COMPETED: ',stn)

ofile_manual_dates.close()
ofile_manual_check.close()

