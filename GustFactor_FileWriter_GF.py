"""" 
program to compute gust factors based on methodology outlined in 
Masters, F.J., Vickery, P.J., Bacon, P. and Rappaport, E.N., 2010. 
Toward objective, standardized intensity estimates from surface wind speed observations. 
Bulletin of the American Meteorological Society, 91(12), pp.1665-1682.

""""

import numpy as np
import math
from scipy.optimize import least_squares

import datetime

import os,sys
import h5py

def vector_average_spd(spd2_list,dir2_list):
	u_list = []
	v_list = []
	for spd in range(len(spd2_list)):
		u_list.append(spd2_list[spd]*math.cos(math.radians(dir2_list[spd])))
		v_list.append(spd2_list[spd]*math.sin(math.radians(dir2_list[spd])))

	vector_mean = (np.asarray(u_list).mean()**2 + np.asarray(v_list).mean()**2)**0.5
	scalar_mean = np.asarray(spd2_list).mean()
 
	vector_mean_dir = math.degrees(np.arctan2(np.asarray(v_list).mean(),np.asarray(u_list).mean()))
	vector_mean_dir = vector_mean_dir%360
	
	return vector_mean,scalar_mean,vector_mean_dir

def yamartino(thetalist):
	s=0
	c=0
	n=0.0

	for theta in thetalist:
		s=s+math.sin(math.radians(theta))
		c=c+math.cos(math.radians(theta))
		n+=1
	s=s/n
	c=c/n
	if (1-(s**2+c**2))>0:
		eps=(1-(s**2+c**2))**0.5
	else:
		eps = 0.0
	sigma=(math.asin(eps)*(1+(2.0/3.0**0.5-1)*eps**3))

	return math.degrees(sigma)


def RA_test (x):

	N = len(x)
	h = np.zeros((N,N))

	count = 0
	for ii in range(N):
		for jj in range (N):
			while ii <jj:
				if x[ii]>x[jj]:
					h[ii,jj]=1
				else:
					h[ii,jj]=0
				count = count + 1
				if count > N:
					break

	A = np.zeros(N)
	
	for ii in range(N):
		A[ii]=sum(h[ii])
	
	A = sum(A)
	meanA = N*(N-1)/4.
	stdA = math.sqrt(N*(2*N+5)*(N-1)/72.)
	z = 1.96

	Asup = math.ceil(z*stdA+meanA-0.5)
	Ainf = math.floor(meanA-z*stdA+0.5)

	if A<Asup and A>Ainf:
		Result = 1
	else:
		Result= 0

	return Result



projectDataDir = '../data/'
dataResDir = 'asos-onemin'
CONUSDataDirRaw = 'asos-onemin-hdf-files'

#state_list = ['GA','NE','NH','NJ','NM','NY','NV','OK','OH','OR','PA','RI','SC','SD','TX','UT','VT','VA','WA','WV','WI','WY']
state_list = ['PA']

dir_list = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW']


Sonic_start = {}
aheight_dict = {}
meta_file = open('ASOS_anem_metadata.txt','r')

while 1:
	line = meta_file.readline()
	if len(line)==0:break
	line = line.strip()
	line = line.split('\t')
	Sonic_start[line[0]]=line[4]
	aheight_dict[line[0]]=int(line[6])

for state in state_list:
	CONUSDataDir = '-hdf-files-with-error-check-and-manual-corrections'
	manual_states = os.listdir(projectDataDir+dataResDir+CONUSDataDir)
	stn_list_all = os.listdir(projectDataDir+'asos-onemin-hdf-files-with-error-check/'+state)

	if state not in manual_states:
		CONUSDataDir = 'asos-onemin-hdf-files-with-error-check'
		stn_list_manual =[]
	else:
		stn_list_manual = os.listdir(projectDataDir+dataResDir+CONUSDataDir+'/'+state)


	for infile in stn_list_all:

		if infile not in stn_list_manual: CONUSDataDir = 'asos-onemin-hdf-files-with-error-check'

		gust_factor_dict = {}

	
		stn = infile[0:4]

		if stn in Sonic_start.keys():
			sonic_date = datetime.datetime(int(Sonic_start[stn][0:4]),int(Sonic_start[stn][5:7]),int(Sonic_start[stn][8:10]))
			aheight = aheight_dict[stn]
		else:
			print ('NO METADATA FOR ', stn)
			continue
				

		ofile_GF= open('../data/CONUS_gust_factors/'+stn+'_GF.txt','w')
		print (projectDataDir+dataResDir+'/'+ state + '/' + infile)

		cmd = 'tar xzf ' + projectDataDir+CONUSDataDirRaw+'/'+ state + '/' + infile
	
		res = os.system(cmd)
		h5out = h5py.File(stn + '-'+dataResDir+'.h5','r+')

		try:
			dates = list(h5out['dates'][:])
	
			windspd = list(h5out['windspd'][:])
			winddir = list(h5out['winddir'][:])
		except:
			try:
				dates = list(h5out['asos_one_min']['dates'][:])
				windspd = list(h5out['asos_one_min']['windspd'][:])
				winddir = list(h5out['asos_one_min']['winddir'][:])
			except:
				print ('PROBLEM WITH RAW HD5 FILE FROM JESS', h5out.keys())
				continue
		cmd = 'tar xzf ' + projectDataDir+CONUSDataDir+'/'+ state + '/' + infile
		res = os.system(cmd)
		h5out = h5py.File(stn + '-'+dataResDir+'.h5','r+')
		try:
			peakdir = list(h5out['asos_one_min']['peakdir'][:])
			peakspd = list(h5out['asos_one_min']['peakspd'][:])
		except:
			try:
				peakdir = list(h5out['peakdir'][:])
				peakspd = list(h5out['peakspd'][:])				
			except:
				print ('PROBLEM WITH FINAL HD5 FILE FROM JESS', h5out.keys())
				continue
			
		peakspd_stdev = np.array(peakspd)[np.array(peakspd)>0].std()
		
		dyr = 0
		pkdir_dict = {}
		pdir_list = []
		pspd_list = []
		spd2_list = []
		dir2_list = []
		date_list = []
		cnt = 0
		for val in range(len(dates)):

			if cnt ==0:
				start_date = datetime.datetime(int(dates[val][0:4]),int(dates[val][4:6]),int(dates[val][6:8]),int(dates[val][8:10]),int(dates[val][10:12]))
		
			cnt = cnt + 1
			new_date = datetime.datetime(int(dates[val][0:4]),int(dates[val][4:6]),int(dates[val][6:8]),int(dates[val][8:10]),int(dates[val][10:12]))

			if new_date < sonic_date:
				atype = 'B'
			else:
				atype = 'V'
				
			if dyr != str(int(dates[val][0:4]))+atype:
				dyr = str(int(dates[val][0:4]))+atype

				
			if new_date == sonic_date:
				spd2_list = []
				dir2_list = []
				pspd_list = []
				pdir_list = []
				cnt = 0
				continue

			if peakdir[val]<0:
				cnt = 0	
				spd2_list = []
				dir2_list = []
				pspd_list = []
				pdir_list = []
				pdir_list = []
				continue
			if peakspd[val] < 0 or peakspd[val] > 5*peakspd_stdev:
				cnt = 0
				spd2_list = []
				dir2_list = []
				pspd_list = []
				pdir_list = []
				continue
			if windspd[val] <0 or windspd[val]>200:
				cnt = 0
				spd2_list = []
				dir2_list = []
				pspd_list = []
				pdir_list = []
				continue
			if winddir[val] <0:
				cnt = 0
				spd2_list = []
				dir2_list = []
				pspd_list = []
				pdir_list = []
				continue
		
		
			if (new_date-start_date).total_seconds()<600:	
				pdir_list.append(peakdir[val])
				pspd_list.append(peakspd[val])
				spd2_list.append(windspd[val])
				dir2_list.append(winddir[val])
				date_list.append(new_date)
				
				continue
			elif (new_date-start_date).total_seconds()>600:
				#print (new_date-start_date,'>600',dates[val],pdir_list,start_date)
				pspd_list = [peakspd[val]]
				pdir_list = [peakdir[val]]
				spd2_list=[windspd[val]]
				dir2_list=[winddir[val]]
				date_list = [new_date]
				start_date = new_date
				
				continue
			elif (new_date-start_date).total_seconds()==600:
				
				if len(pdir_list)<10:
					print ('SHORT', date_list,len(pspd_list),new_date,start_date)
					pspd_list = [peakspd[val]]
					pdir_list = [peakdir[val]]
					spd2_list=[windspd[val]]
					dir2_list=[winddir[val]]
					start_date = new_date
					continue
		
				if len(pdir_list)!=10:
					print (pdir_list,new_date,start_date)
					duh
		
				### only want very other 2 minute wind in the to avoid double-counting winds in the 1-minuta data.  Use np.mean(art[::2])
		
		
				pspd = max(pspd_list)
		
				spd2_stdev = np.std(spd2_list)
				
				pspd_list = np.asarray(pspd_list)[::2]
				pdir_list = np.asarray(pdir_list)[::2]
				spd2_list = np.asarray(spd2_list)[::2]
				dir2_list = np.asarray(dir2_list)[::2]
		
		
				spd2 = spd2_list.mean()
		
				too_variable = 0
				std = yamartino(dir2_list)
		
				vect,scale,vect_dir = vector_average_spd(spd2_list,dir2_list)
		
				p_vect,p_scale,p_vect_dir = vector_average_spd(pspd_list,pdir_list)
		
				pdir = p_vect_dir
				dir2 = vect_dir
		
				if std>10:
					too_variable = 1
		
				if abs(vect-scale) > 1: ## one kt
					too_variable = 1
		
		
				RAtest = RA_test(spd2_list)
		
				meanRemovedDir = dir2_list-vect_dir
		
				peak_meanRemovedDir = pdir_list-p_vect_dir
		
				alongCompSpd = spd2_list*np.cos(np.radians(meanRemovedDir))
				acrossCompSpd = spd2_list*np.sin(np.radians(meanRemovedDir))		
		
				p_alongCompSpd = pspd_list*np.cos(np.radians(peak_meanRemovedDir))
				p_acrossCompSpd = pspd_list*np.sin(np.radians(peak_meanRemovedDir))
		
				spd2AlongComp = alongCompSpd.mean()
				pspdAlongComp = max(p_alongCompSpd)
				
				if abs(np.mean(acrossCompSpd))>0.01:print ('ERROR ACROSS COMP SPEED DOES NOT EQUAL ZERO', acrossCompSpd)
				if abs(np.mean(p_acrossCompSpd))>0.01:print ('ERROR ACROSS COMP PEAK SPEED DOES NOT EQUAL ZERO', p_acrossCompSpd)
		
				spd2_list = []
				dir2_list = []
				pspd_list = []
				pdir_list = []
				pdir_list.append(peakdir[val])
				pspd_list.append(peakspd[val])
				spd2_list.append(windspd[val])
				dir2_list.append(winddir[val])
				date_list = [new_date]
		
				start_date = new_date
		
			if spd2*(0.51444)< 5 or too_variable == 1:         #5 is windspeed in m/s * 0.51444 converts kt to m/s
				continue  
		
			pdir =  p_vect_dir

		
			if pdir<=0:continue
		
			if pdir <= 11.25: 
				pcom = 'N'
			elif pdir >11.25 and pdir <=33.75:
				pcom = 'NNE'
			elif pdir >33.75 and pdir <=56.25:
				pcom = 'NE'
			elif pdir >56.25 and pdir <=78.75:
				pcom = 'ENE'
			elif pdir > 78.75 and pdir <=101.25:
				pcom = 'E'
			elif pdir >101.25 and pdir <=123.75:
				pcom = 'ESE'
			elif pdir >123.75 and pdir <=146.25:
				pcom = 'SE'
			elif pdir >146.25 and pdir <=168.75:
				pcom = 'SSE'
			elif pdir >168.75 and pdir <=191.25:
				pcom = 'S'
			elif pdir >191.25 and pdir <=213.75:
				pcom = 'SSW'
			elif pdir >213.75 and pdir <=236.25:
				pcom = 'SW'
			elif pdir >236.25 and pdir <=258.75:
				pcom = 'WSW'
			elif pdir >258.75 and pdir <=281.25:
				pcom = 'W'
			elif pdir >281.25 and pdir <=303.75:
				pcom = 'WNW'
			elif pdir >303.75 and pdir <=326.25:
				pcom = 'NW'
			elif pdir >326.25 and pdir <=348.75:
				pcom = 'NNW'
			elif pdir>348.47 and pdir <=360:
				pcom = 'N'	
		
		
			
			if dyr not in pkdir_dict:
				pkdir_dict[dyr]={}
			if pcom not in pkdir_dict[dyr]:
				pkdir_dict[dyr][pcom]=[[],[],[]]
		
			pkdir_dict[dyr][pcom][0].append(pspd)
		
			pkdir_dict[dyr][pcom][1].append(spd2AlongComp)
		
			pkdir_dict[dyr][pcom][2].append(pspd/spd2AlongComp)

		avail_yrs = list(pkdir_dict.keys())
		avail_yrs.sort
				
		for val in range(len(avail_yrs)):
			yr = avail_yrs[val]
			yr_end = yr[0:4]+'-12-31'
			yr_start = yr[0:4]+'-01-01'
			
			if avail_yrs[val][-1]=='B':
				instrum = 'Belfort'
			else:
				instrum = 'Vaisala'
				
			if val+1 != len(avail_yrs):
				if yr[0:4] == avail_yrs[val+1][0:4]:
					yr_start = yr[0:4]+'-01-01'
					yr_end = str(sonic_date.year) + '-'+str(sonic_date.month).zfill(2)+'-'+str(sonic_date.day).zfill(2)
					instrum = 'Belfort'
			if val !=0:
				if yr[0:4] == avail_yrs[val-1][0:4]:
					yr_start = str(sonic_date.year) + '-'+str(sonic_date.month).zfill(2)+'-'+str(sonic_date.day).zfill(2)
					yr_end = yr[0:4]+'-12-31'
					instrum = 'Vaisala'
					
		# ======================================
		# Inputs variables
		# =================================
			for dir in dir_list:
				U_T = np.nan
				u_t = np.nan
				U_T_num = 0
				u_t_num = 0
		
				if dir in pkdir_dict[yr].keys():
					U_T = np.array(pkdir_dict[yr][dir][1])*0.51444   ### convert kt to ms-1
					u_t = np.array(pkdir_dict[yr][dir][0])*0.51444   ### convert kt to ms-1
					U_T_num = len(pkdir_dict[yr][dir][1])
					u_t_num = len(pkdir_dict[yr][dir][0])
				
					GF = u_t/U_T
		
					gust_factor_dict[dir]=(np.mean(GF),len(GF))
				else:
					gust_factor_dict[dir]=(np.nan,0)

						
			if yr[0:4] == '2000':
				ofile_GF.write('ICAO\tStartDate\tEndDate\tAnemometer\tHeight\tHeightVer\t')
				for dir in dir_list:
					ofile_GF.write(dir+'_GF\t')
				for dir in dir_list:
					ofile_GF.write(dir+'_N\t')
				ofile_GF.write('\n')
		
			for dir in dir_list:
				if dir == dir_list[0]:
					ofile_GF.write('%4s\t%s\t%s\t%s\t%2d\tverif\t%4.2f'%(stn,yr_start,yr_end,instrum,aheight,gust_factor_dict[dir][0]))
				else:
					ofile_GF.write('\t%4.2f'%(gust_factor_dict[dir][0]))
		
			for dir in dir_list:
				if dir==dir_list[-1]:
					ofile_GF.write('\t%4d\n'%(int(gust_factor_dict[dir][1])))
				else:
					ofile_GF.write('\t%4d'%(int(gust_factor_dict[dir][1])))
