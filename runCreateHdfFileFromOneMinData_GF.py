'''Convert one-minute text file to HDF5 format.

	- primarily converts one-minute text files to HDF5 format
	- sets minutes with text file formatting errors to missing

	bnb2
'''
import os,sys
import h5py
import fileinput
import datetime
import numpy as np

##############################

stateList = ['PA']
# Project data directory
projectDataDir = '../data/'
outDataDir = '../data/'
# data resolution directory name
dataResDir = 'asos-onemin'
miss = '-99'
completedStations = []

def createHdf5(inDirBase,outDirBase,stn):
	################################
	### WIND DIRECTION AND SPEED
	################################

	### get file list for this station, and create fileinput object for looping lines of all files
	inDir = inDirBase + stn + '/'
	fileList = os.listdir(inDir)

	if len(fileList)==0:
		print ('No files available for ',state,':',stn)
		print ('END: ',state,':',stn)
		return 1
	fileList.sort()
	fileList = [inDir+f for f in fileList if 'pg1-'+stn in f]
	allfiles = fileinput.input(fileList)


	### loop all lines
	data_wind = []

	lastValidTime =  '199912312359'
	for line in allfiles:
		
		splitLine = line.strip().split()
		
		### four wind elements are in columns 68-90 (or, column indices 67-89)
		windElements = line[67:90].strip().split()
		### skip this time if elements are not as we expect
		# 1. line too short
		if len(splitLine)<8:
		#	print('too short ',splitLine[1][3:15])
			continue
		# 2. line too long
		if len(splitLine)>13:   # Changed from 12 to 13 2/10/24 atd
		#	print('too long ',splitLine[1][3:15],len(splitLine))
			continue
		# 3. unexpected number of wind elements
		if len(windElements)!=4:
		#	print('unexpected number ',splitLine[1][3:15])
			continue
		# 4. unexpected format
		if not (
			stn in splitLine[0] and
			stn[1:] in splitLine[1] and
			'.' in splitLine[2] and
			(splitLine[3]=='N' or splitLine[3]=='D') and
			windElements[-4].isdigit() and
			windElements[-3].isdigit() and
			windElements[-2].isdigit() and
			windElements[-1].isdigit()
		): 
		#	print('bad format ',splitLine[1][3:15])
			continue

		### extract variables of interest
		time = splitLine[1][3:15]


		winddir = windElements[-4]
		windspd = windElements[-3]
		peakdir = windElements[-2]
		peakspd = windElements[-1]



		if lastValidTime is not None:
			### determine how many minutes since last valid report
			dtTime = datetime.datetime.strptime(time,"%Y%m%d%H%M")
			dtLastValidTime = datetime.datetime.strptime(lastValidTime,"%Y%m%d%H%M")
			minutes = divmod((dtTime - dtLastValidTime).total_seconds(), 60)
			minutesSinceLastValidTime = int(minutes[0])

			### skip this line if it contains a repeated time
			if minutesSinceLastValidTime<=0: continue

			### fill missing minutes, if they exist since last valid time
			if minutesSinceLastValidTime!=1:
				for m in range(1,minutesSinceLastValidTime):
					dtThisTime = dtLastValidTime + datetime.timedelta(minutes=m)
					thisTime = dtThisTime.strftime("%Y%m%d%H%M")
					data_wind.append([thisTime]+[int(miss)]*len(windElements))
			if time[:4]!=lastValidTime[:4]: print (time[:4])
			sys.stdout.flush()
		else:
			print (time[:4])
			sys.stdout.flush()

		### fill data for current minute
		data_wind.append([time,int(winddir),int(windspd),int(peakdir),int(peakspd)])
		lastValidTime = time

	allfiles.close()



	dates_wind = [ item[0] for item in data_wind ]

	idxWindStart=0
	idxWindEnd=len(dates_wind)-1

	##########################
	### create HDF5 file
	##########################
	h5file_out = h5py.File(outDirBase+stn+'-asos-onemin.h5','w')

	### dates   chaged from strings to ints atd 2/10/24
	data = [ int(item[0]) for item in data_wind[idxWindStart:idxWindEnd+1] ]
	shapeOfHdf5Dataset = (len(data),)
	dset_dates = h5file_out.create_dataset('dates', shapeOfHdf5Dataset, dtype='S12')
#	dset_dates = h5file_out.create_dataset('dates', shapeOfHdf5Dataset,dtype='s12')
	
	dset_dates[:] = data
	del data

	### wind direction for 2-minute average
	data = [ item[1] for item in data_wind[idxWindStart:idxWindEnd+1] ]
	shapeOfHdf5Dataset = (len(data),)
	dset_wdir = h5file_out.create_dataset('winddir', shapeOfHdf5Dataset, dtype='<f8')
	dset_wdir[:] = data
	del data

	### wind speed for 2-minute average
	data = [ item[2] for item in data_wind[idxWindStart:idxWindEnd+1] ]
	shapeOfHdf5Dataset = (len(data),)
	dset_wspd = h5file_out.create_dataset('windspd', shapeOfHdf5Dataset, dtype='<f8')
	dset_wspd[:] = data
	del data

	### peak wind direction for max 5-sec average over 1 minute
	data = [ item[3] for item in data_wind[idxWindStart:idxWindEnd+1] ]
	shapeOfHdf5Dataset = (len(data),)
	dset_pdir = h5file_out.create_dataset('peakdir', shapeOfHdf5Dataset, dtype='<f8')
	dset_pdir[:] = data
	del data

	### peak wind speed for max 5-sec average over 1 minute
	data = [ item[4] for item in data_wind[idxWindStart:idxWindEnd+1] ]
	shapeOfHdf5Dataset = (len(data),)
	dset_pspd = h5file_out.create_dataset('peakspd', shapeOfHdf5Dataset, dtype='<f8')
	dset_pspd[:] = data
	del data

	h5file_out.close()

	del data_wind

	allfiles.close()

	# tar and zip
	srcPath = os.getcwd()
	os.chdir(outDirBase)
	cmd = 'tar czf '+stn+'-asos-onemin.h5.tgz '+stn+'-asos-onemin.h5'
	res = os.system(cmd)
	os.remove(stn+'-asos-onemin.h5')
	os.chdir(srcPath)

	#ofile.close()
	return 0

for state in stateList:
	# data source directory
	inDirBase = projectDataDir+dataResDir+'/'+state+'/'
	# HDF output directory
	hdfBase = outDataDir+dataResDir+'-hdf-files/'+state+'/'
	# create output directory if it doesn't exist
	if not os.path.exists(hdfBase):
		res = os.system('mkdir -p '+hdfBase)
	# get station list for state
	stnList = os.listdir(inDirBase)


	if len(stnList)==0:
		print ('No stations available for ',state)
		print ('END: ',state)
		continue

	for stn in stnList:
		# bnb2: these stns are already completed, skip
		if stn in completedStations: continue
		print ('START ',state,':',stn)
		sys.stdout.flush()

		res = createHdf5(inDirBase,hdfBase,stn)

		print ('COMPETED: ',stn)
