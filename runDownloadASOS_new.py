'''Download ASOS one-minute data for WDR project.

	Download ASOS one-minute wind/precip data from FTP source.
	All stations from the SE WDR project can be re-downloaded by using the
	stations listed in the GitHub station inventory.

	Execution time: about 6-7 minutes per station

	bnb2
'''
##
## install Hombrew https://docs.brew.sh
## then run brew install wget
## pip install does not work bc this NOT the python module wget
##
#program to get new ASOS data: 2022-present

import os
import requests
import json

stn_in = open('../all_stns.txt','r')

stnList = stn_in.readlines()

#url = "https://www.ncei.noaa.gov/data/automated-surface-observing-system-one-minute-pg1/access/2024/08/asos-1min-pg1-KBGM-202408.dat"
#for wind data only
url = "https://www.ncei.noaa.gov/data/automated-surface-observing-system-one-minute-pg1/access/"
#for precip data
#url = "https://www.ncei.noaa.gov/data/automated-surface-observing-system-one-minute-pg2/access/"

req = requests.get(url)
data = req.text
data = data.split('\n')
#print(data)

for stn_id in range(1,len(stnList)):
	stn = stnList[stn_id].strip()
	stn = stn.split('\t')
	print (stn[1])
	state = stn[0]
	#outDir = '/Volumes/exthd/data_new/asos-onemin/'+state+'/wind/'+stn[1]
	outDir = '../data/test/asos-onemin/'+state+'/'+stn[1]
	if not os.path.exists(outDir):
		res = os.system('mkdir -p '+outDir)
	#filePattern = '*'+stn[1]+'*.dat'
	
	years = list(range(1999,2025))
	
	for year in years:
		print(year)
		months = list(range(1,13))
		#months = ['01','02','03','04','05','06','07','08','09','10','11','12']
		for month in months:
			print(month)
			nurl = url+str(year)+'/'+str(month)+'/asos-1min-pg1-'+stn[1]+'-'+str(year)+str(month)+'.dat'
		#	nurl = url+year+'/'+month+'/asos-1min-pg2-'+stn[1]+'-'+year+month+'.dat'
			#wget downloads the file
			#-P followed by folder tells where to put the file
			cmd = 'wget ' +nurl+' -P '+outDir
			#res = response
			res = os.system(cmd)
	
	
	
	
	#mainOptions = '-r -q -nH --cut-dirs=4 -np'
	#cmd = 'wget '+mainOptions+' '+ftploc+' -P '+outDir+' -A '+filePattern
	#res = os.system(cmd)
	
	
#	ftploc = 'ftp://ftp.ncdc.noaa.gov/pub/data/asos-onemin'
#mainOptions = '-r -q -nH --cut-dirs=4 -np'