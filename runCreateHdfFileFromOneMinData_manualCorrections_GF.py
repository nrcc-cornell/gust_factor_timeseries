'''Apply manual corrections to one-minute HDF5 files.

	- read text file containing minutes that require correction after manual inspection

	bnb2
'''
import os,sys
import h5py
import numpy as np

try :
	import json
except ImportError :
	import simplejson as json

##############################
# Project data directory
projectDataDir = '../data/'
# data resolution directory name
dataResDir = 'asos-onemin'
miss = '-99'

# Get all elements to be corrected (set to missing).
# These are the items that did not pass manual testing.
# sample line of file
# state,ICAO,varname,minute,wdir,wspd,prcp,'obs',obsvalue
# FL,KPGD,prcp,200608161148,118.0,5.0,0.36,obs,0.8598633
# 'obsvalue' is either daily radar-guided precip total or max daily wind gust from hourly reports.
itemsToEditFile = open('manual_data_dates_long_POR_invalidSaved.txt','r')
linesOfFile = itemsToEditFile.readlines()
itemsToEditFile.close()

dataToEdit = {}
for l in linesOfFile:
	listOfElements = l.strip().split(',')
	state = listOfElements[0]
	if state not in dataToEdit: dataToEdit[state] = {}
	stn = listOfElements[1]
	if stn not in dataToEdit[state]: dataToEdit[state][stn] = []
	dataToEdit[state][stn].append(listOfElements[2:7])

# loop all states that have at least one data element to be corrected

for state in dataToEdit:

	# data source directory
	inDirBase = projectDataDir+dataResDir+'-hdf-files-with-error-check/'+state+'/'
	# HDF output directory
	outDirBase = projectDataDir+dataResDir+'-hdf-files-with-error-check-and-manual-corrections/'+state+'/'
	# create output directory if it doesn't exist
	if not os.path.exists(outDirBase):
		res = os.system('mkdir -p '+outDirBase)

	# loop all stations in this state that have at least one data element to be corrected
	for stn in dataToEdit[state]:

		# untar src file to output directory. We will edit the copy in the output directory.
		cmd = 'tar xzf '+inDirBase + stn + '-'+dataResDir+'.h5.tgz --directory '+outDirBase
		res = os.system(cmd)

		# open output file for reading and writing. We edit data in this section
		h5out = h5py.File(outDirBase + stn + '-'+dataResDir+'.h5','r+')
		dates = list(h5out['asos_one_min/dates'][:])

		# loop through each data element to be edited for this stn
		for varname,timeToEdit,wdirOld,wspdOld,prcpOld in dataToEdit[state][stn]:
			print (state,stn,varname,timeToEdit,wdirOld,wspdOld,prcpOld)

			idxToEdit = None
			idxToEdit = dates.index(timeToEdit.encode('utf-8'))
			if varname=='wind':
				if np.isclose( np.array([float(wspdOld)]), np.array([float(h5out['asos_one_min/peakspd'][idxToEdit])]) ):
					h5out['asos_one_min/peakspd'][idxToEdit] = float(miss)
					print ('    SUCCESS: peakspd set to missing at index '+str(idxToEdit))
				else:
					print ('    ERROR : expected peakspd values do not match, not removing data')

		h5out.close()

		# tar and zip output file
		srcPath = os.getcwd()
		os.chdir(outDirBase)
		cmd = 'tar czf '+stn+'-'+dataResDir+'.h5.tgz '+stn+'-'+dataResDir+'.h5'
		res = os.system(cmd)
		os.remove(stn+'-'+dataResDir+'.h5')
		os.chdir(srcPath)

