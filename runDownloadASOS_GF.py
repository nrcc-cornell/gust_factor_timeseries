'''Download ASOS one-minute data for WDR project.

	Download ASOS one-minute wind/precip data from two sources at NCEI (FTP and HTTP) check for duplicate files to create the most complete set of data files as possible .

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

url = "https://www.ncei.noaa.gov/data/automated-surface-observing-system-one-minute-pg1/access/"


ftploc = 'ftp://ftp.ncdc.noaa.gov/pub/data/asos-onemin'
mainOptions = '-r -q -nH --cut-dirs=4 -np'

req = requests.get(url)
data = req.text
data = data.split('\n')


for stn_id in range(1,len(stnList)):
	stn = stnList[stn_id].strip()
	stn = stn.split('\t')
	print (stn[1])
	state = stn[0]
	outDir = '../data/test/asos-onemin/'+state+'/'+stn[1]
	if not os.path.exists(outDir):
		res = os.system('mkdir -p '+outDir)

	years = list(range(1999,2025))
	
	for year in years:
		print(year)
		months = list(range(1,13))
		for month in months:
			print(month)
			nurl = url+str(year)+'/'+str(month).zfill(2)+'/asos-1min-pg1-'+stn[1]+'-'+str(year)+str(month).zfill(2)+'.dat'
			cmd = 'wget ' +nurl+' -P '+outDir
			res = os.system(cmd)
	
	outDir = '../data/test/asos-onemin/'+state+'/'+stn[1]
	if not os.path.exists(outDir):
		res = os.system('mkdir -p '+outDir)
	filePattern = '64050'+stn[1]+'*.dat'
	cmd = 'wget '+mainOptions+' '+ftploc+' -P '+outDir+' -A '+filePattern
	print (cmd)
	res = os.system(cmd)
	
	
	fileList = os.listdir(outDir)
	
	if len(fileList)!=0:
		for  filename in fileList:
			if filename[0:4]==stn[1]:
				continue
			elif filename[0:5]=='64050':
				if 'asos-1min-pg1-'+ filename[5:9]+'-'+filename[9:] in fileList:
					res = os.system('rm '+ outDir+'/'+filename)	
					continue
				else:
					cmd = 'mv '+outDir+'/'+filename+' '+outDir+'/asos-1min-pg1-'+ filename[5:9]+'-'+filename[9:]
					res = os.system(cmd)
