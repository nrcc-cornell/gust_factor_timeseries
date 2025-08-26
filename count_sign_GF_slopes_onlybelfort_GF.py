''''
Perform piecewise regression on the annual gust factor time series.
Output is a large json file with adjusted gust factor time series for each station-wind direction combination. 

'''
import os
import json
import copy 
import h5py
import datetime,random,string
import numpy as np
import scipy.stats as ss
import piecewise_regression
import pwlf

def getStnsToAnalyze(p):
	'''Get available states/stations to analyze.

		p: path to HDF5 files of data
		output: key-value object, keyed by state with station list values
	'''
	if p[-1]!='/': p+'/'
	states = os.listdir(p)
	states.sort()
	output = {}
	for state in states:
		if state.count('.')>0:continue
		stnFiles = os.listdir(p+state)
		stnFiles.sort()
		output[state] = [ f.split('-')[0] for f in stnFiles ]
	return output
	

def readGustFactorData_art(f,miss):
	'''Read gust factor and anemometer metadata.

		Gust factors are provided by time period and direction.
		Resulting data (result) is keyed by:
			StartDate: starting dates of intervals (YYYY-MM-DD)
			EndDate: ending dates of intervals (YYYY-MM-DD)
			Anemometer: type of anemometer during intervals ('Belfort','Vaisala')
			Height: height of anemometer during intervals (m)
			gust_factor: gust factors calculated for each interval, keyed by direction
	'''

	directions = ['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW']
	colsToGet = ['StartDate','EndDate','Anemometer','Height']+[direction+'_GF' for direction in directions]

	# The reference data in excel files is 1900-01-01. We'll need this info to calculate the actual dates later.
	refDateForRaw_dt = datetime.datetime.strptime('1900-01-01','%Y-%m-%d')

	# read all data from excel file into list

	if f[-13]!='s':
		f_orig = f[0:25]+f[-12:]
	else:
		f_orig = f
	
	ifile = open(f_orig,'r')
	data = []

	header_line = ['ICAO','StartDate','EndDate','Anemometer','Height','HeightVer','N_GF','NNE_GF','NE_GF','ENE_GF','E_GF','ESE_GF','SE_GF','SSE_GF','S_GF','SSW_GF','SW_GF','WSW_GF','W_GF','WNW_GF','NW_GF','NNW_GF','N_N','NNE_N','NE_N','ENE_N','E_N','ESE_N','SE_N','SSE_N','S_N','SSW_N','SW_N','WSW_N','W_N','WNW_N','NW_N','NNW_N']

	while 1:
		line = ifile.readline()
		if len(line)==0:break
		line = line.strip()
		line = line.split('\t')

		data.append(line)
### this is a hack since some GF files do not have a header line
		if data[0][0]!='ICAO':
			print('station with no header')
			data = [header_line]
			data.append(line)

	if f_orig!=f:
		ifile = open(f,'r')
		while 1:
			line = ifile.readline()
			if len(line)==0:break
			line = line.strip()
			line = line.split('\t')
			
			if line[0]!='ICAO':
				if line[1][0:4]=='2022':
					del data[-1]
				data.append(line)				
				

	for irow in range(1,len(data)):
		data[irow][6:] = map(float,data[irow][6:])
		data[irow][4]= float(data[irow][4])
		if data[irow][4]<-99:data[irow][4]=10   #### set missing anemometer heights to 10. by default

	header = data[0]
	# For each time period, check sample sizes by direction.
	# If more than half are zero, assign empty string to all GF and N during that time period.
	# Empty strings are used to indicate empty cells in the excel files, and interpreted as missing later.
	for irow in range(1,len(data)):
		if data[irow][-16:].count(0)>8: data[irow] = data[irow][:-32] + 32*['']

	# Organize data into object.
	result = {}
	result['gust_factor'] = {}
	for c in colsToGet:
		#print (c,header)

		icol = header.index(c)
		#print (icol)
		if 'Date' in c:
			result[c] = [(data[irow][icol]) for irow in range(1,len(data)) ]
		elif 'GF' in c:
			direction = c.split('_')[0]
			icol_N = header.index(direction+'_N')
			result['gust_factor'][direction] = [ data[irow][icol] if data[irow][icol]!='' and data[irow][icol_N]>=30 else miss for irow in range(1,len(data)) ]
		else:
			result[c] = [ data[irow][icol] if data[irow][icol]!='' else miss for irow in range(1,len(data)) ]
	#s
	return result
	
def smoothGustFactorData(d,refDateForRegression,miss):

	
	'''Smooth gust factor data.
		Determine significant linear relationships, and calculate reference equations dependent on days since refDateForRegression.

		d: gust factor data read from excel file (readGustFactorData)
	'''

	# The reference date for gust factor regression equations is refDateForRegression.
	refDateForRegression_dt = datetime.datetime.strptime(refDateForRegression,'%Y-%m-%d')

	result = {}

	# calculate 'smoothed' date (keyed by 'SmoothDate')
	# SmoothDate is the middle date of the period (between 'StartDate' and 'EndDate').
	# SmoothDate has units of 'days since refDateForRegression'.
	smoothDate = []
	for sd,ed in zip(d['StartDate'],d['EndDate']):
		sd_daysSinceRef = (datetime.datetime.strptime(sd,'%Y-%m-%d') - refDateForRegression_dt).days
		ed_daysSinceRef = (datetime.datetime.strptime(ed,'%Y-%m-%d') - refDateForRegression_dt).days
		# Average of sd and ed will represent middle of the time period.
		# Using mid-period date is best when date is associated with mean gust factor over same period
		smoothDate.append( (sd_daysSinceRef + ed_daysSinceRef)/2. )
	result['SmoothDate'] = smoothDate

	# Calculate smoothed gust factors by linear regression,
	# for each wind direction (16-pt) and anemometer type ('Belfort','Vaisala').
	# y = m*x + b, where
	#	x = days since refDateForRegression
	#	y = gust factor
	#	m = slope of regression fit
	#	b = y-intercept
	# if N<10: m=0, b=mean (no significance test performed due to sample size)
	# if N>=10:
	#	if slope significantly nonzero: m,b used from regression
	#	if slope not significantly diff from zero: m=0, b=mean
	result['gust_factor_smooth'] = {}
	base_data = {}
	for anemometerType in ['Belfort','Vaisala']:
		anem_mask = np.logical_not( np.array(d['Anemometer'])==anemometerType )
		start_date_mask = np.array([elem[-6:] for elem in d['StartDate']])!='-01-01'
		end_date_mask = np.array([elem[-6:] for elem in d['EndDate']])!='-12-31'
		for direction in d['gust_factor']:
			base_data[direction] = {}
			gf_miss_mask = np.array(d['gust_factor'][direction])==miss
			gf_extreme_mask = np.array(d['gust_factor'][direction])>=3.0
			mask_full = anem_mask | gf_miss_mask | gf_extreme_mask | start_date_mask | end_date_mask
			x_ma = np.ma.array(smoothDate,mask=mask_full)
			y_ma = np.ma.array(d['gust_factor'][direction],mask=mask_full)
			x = x_ma[x_ma.mask==False].data
			y = y_ma[y_ma.mask==False].data
			
			base_data[direction]['x']=x
			base_data[direction]['y']=y
				
			if len(y)>=5:
				res = ss.linregress(x,y)
				# gf slope of -8.219e-6/day is -0.003/yr, or -0.03/decade
				# If slope smaller than this magnitude, it is acceptable to use mean value (per Forrest).
				# If slope is larger, and significant, use the regression.
				if res.pvalue<0.025 and abs(res.slope)>8.219e-6:    
					# slope is significantly nonzero and over magnitude threshold
					m = res.slope
					b = res.intercept
				else:
					m = 0
					b = np.mean(y)
			elif len(y)<5 and len(y)>0:
				m = 0
				b = np.mean(y)
			else:
				m = miss
				b = miss
			if direction not in result['gust_factor_smooth']: result['gust_factor_smooth'][direction] = {}
			result['gust_factor_smooth'][direction][anemometerType] = {'m':m, 'b':b}

			#print(m,b)
	return result,base_data

def readWeatherData(f,usingMetricUnits,miss):
	'''Read one-minute weather observations from HDF5 file.
	'''

	### temporary directory to work in
	tmpRandom = ''.join(random.choices(string.ascii_uppercase,k=6))
	tmpDir = './tmp-'+tmpRandom+'/'
	if not os.path.isdir(tmpDir): res = os.system('mkdir -p '+tmpDir)

	### untar file
	#cmd = 'tar xzf ' + inDirBase + stn+'-asos-onemin.h5.tgz --directory ' + inDirBase
	cmd = 'tar xzf ' + f + ' --directory ' + tmpDir
	res = os.system(cmd)
	fUntarred = tmpDir + f.split('/')[-1][:-4]

	### input hdf5 file : get all dates, wind, prcp data in array
	h5file = h5py.File(fUntarred,'r')
	dates = np.array(h5file['asos_one_min/dates'])
	winddir_np = np.array(h5file['asos_one_min/peakdir'])
	wind_np = np.array(h5file['asos_one_min/peakspd'])
	stn_lon = h5file['asos_one_min'].attrs['longitude']
	stn_lat = h5file['asos_one_min'].attrs['latitude']
	stn_name = h5file['asos_one_min'].attrs['name']
	stn_state = h5file['asos_one_min'].attrs['state']
	stn_id = h5file['asos_one_min'].attrs['ID']
	h5file.close()

	if usingMetricUnits:
		### convert wind speed from kts to m/s
		kts_to_mps_factor = 1./1.944
		wind_miss_mask = wind_np==miss
		wind_ma = np.ma.array(wind_np,mask=wind_miss_mask)*kts_to_mps_factor
		wind_ma.set_fill_value(miss)
		wind_np = wind_ma.data
		

	### remove untarred file and tmp directory
	res = os.remove(fUntarred)
	res = os.rmdir(tmpDir)

	results = {
		'meta': {'id':stn_id, 'state':stn_state, 'name':stn_name, 'lon':stn_lon, 'lat':stn_lat},
		'data': {'dates':dates, 'wdir':winddir_np, 'wspd':wind_np}
	}

	return results
	


def find_indices_within_range(lst, range_value):
	indices = []
    
	for i in range(len(lst)):
		for j in range(i + 1, len(lst)):
			if abs(int(lst[i]) - int(lst[j])) <= range_value:
				indices.append((i, j))
    
	return indices




#stnsToRun = ['KHSE','KSAV','KMIA','KPNS','KBTR','KCRP','KELP','KOKC','KBNA']
forceQuadraticStnList = ['KHSE','KELP','KFOK']
#stnsToRun = ['KFOK',]
skipStns = ['KDHN','KSTS','KCEC','KPSP','KDEN','KSHV','KCMA']

CONUS2_states = ['AL','AZ','AR','CA','CO','CT','DE','FL','GA','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MO','MS','MT','NC','ND','NE','NH','NJ','NM','NV','NY','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VA','VT','WA','WI','WY','WV']

#CONUS2_states = ['AL']

#readGustFactorData_artCONUS2_states = ['AL','AZ','AR','CA','CO','CT']

refDateForRegression = '2000-01-01'
miss = -99.
minutesInWindow = 4.*24.*60.
minutesInOneYear = 365.25*24.*60.
# 'wspd','wspd_adj','prcp'
wspdKeyToAnalyze = 'wspd_adj'


gfSrc = './data/CONUS_gust_factors/'


#wxSrc = './data/asos-onemin-hdf-files-with-error-check-and-manual-corrections/'
wxSrc = './data/asos-onemin-hdf-files-with-error-check-and-manual-corrections/'

useMetric = False
includeTitleAndLegend = True

kts_to_mps_factor = 1./1.944
in_to_mm_factor = 25.4

min_break_duration = 2

wspdKeyToAnalyze=='wspd_adj'


stations = getStnsToAnalyze(wxSrc)
gf_data = {}
wxMeta = {}

num_stn_trend=0
num_stn_flat = 0
trend_slope_list= []
GF_diff_list = []
dir_cnts = {}
direct_trends = {}

total_breaks = 4 
break_summary = {}

for tbrk in range(total_breaks+1):
	break_summary[tbrk]={}
	
for state,stns in stations.items():
	if state not in CONUS2_states:continue
	for stn in stns:
		#if stn!='KSAC':continue
		if stn[0:3]=='.DS':continue
		#if stn=='KGZH':break
		print(gfSrc+stn+'_GF.txt')
		if stn in skipStns:continue
		if 'stnsToRun' in globals():
			if stn not in stnsToRun: continue
		#print(state,stn,datetime.datetime.now())
		
		direct_trends[stn]={}
		# read weather data
		wxData = readWeatherData(wxSrc+state+'/'+stn+'-asos-onemin.h5.tgz',useMetric,miss)
		wxMeta[stn] = wxData['meta']
		lat=wxMeta[stn]['lat']
		lat = str(5*round(float(lat)/5))
		
		if wspdKeyToAnalyze=='wspd_adj':
			# read gust factors and create smoothing equations
			try:
				gf_data_raw = readGustFactorData_art(gfSrc+stn+'_GF.txt',miss)
			except:
				print('NO GUST FACTOR FILE',stn)
				continue
				
			gf_data_smooth,base_data = smoothGustFactorData(gf_data_raw,refDateForRegression,miss)

			counter_trnd = 0
			for key in gf_data_smooth['gust_factor_smooth']:
				counter_trnd = 0  #### NEED TO SEE IF THIS BELONGS HERE OR CAN BE OMITTED	
				if stn==stns[0]:
					if key not in dir_cnts.keys():
						dir_cnts[key]=1        ### made this 1 so that it always calculates slopes for pre Vaisala period
				direct_trends[stn][key]=1      ### made this 1 so that it always calculates slopes for pre Vaisala period
				if gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['m']!=0 and gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['m']>-99:
					#direct_trends[stn][key]=direct_trends[stn][key]+1
					#dir_cnts[key]=dir_cnts[key]+1
					counter_trnd = 1
					trend_slope_list.append(gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['m']*365.25)   ## convert from trend per day to trend per year
					GF_diff_list.append((gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['m']*base_data[key]['x'][-1]+gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['b'])-(gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['m']*base_data[key]['x'][0]+gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['b']))
					#print (gf_data_smooth['gust_factor_smooth'][key]['Vaisala'])
					#print(gf_data_raw,refDateForRegression,miss)
					#input('look')		
					
				#	if GF_diff_list[-1]>0.25:print(key,stn,GF_diff_list[-1],base_data[key]['x'],base_data[key]['y'])
		
			#	if stn=='KEET':
			#		print(key,counter_trnd,'counter',gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['m'],direct_trends[stn][key])	


			Vaisala_start = gf_data_raw['Anemometer'].index('Vaisala')+1
			
			#if Vaisala_start==1:Vaisala_start=0
			
			#print(Vaisala_start,'start')

			if counter_trnd == 1: 
				num_stn_trend= num_stn_trend+1
			
			for key in gf_data_raw['gust_factor']:
				
				if direct_trends[stn][key]==1:
					#print(gf_data_raw['gust_factor'][key])
					if  dir_cnts[key]==0:
						#print(stn,key,gf_data_raw['gust_factor'])
						#input('no key')	
						continue
					
			#		print(gf_data_smooth['gust_factor_smooth'][key]['Vaisala']['m']*365.25)   ## ## convert from trend per day to trend per year
					#print(gf_data_raw)
					#print(refDateForRegression)
					
					
					y = gf_data_raw['gust_factor'][key][0:Vaisala_start-2]
					X = list(range(1,len(y)+1))
					y_hold = copy.deepcopy(y)
					X_hold = copy.deepcopy(X)
		
					old_y = np.array(y_hold)
					valid_mask = old_y != -99
					cleaned = old_y[valid_mask]
					original_indices = np.where(valid_mask)[0]

					#print ('first orginal',original_indices)
					#print (valid_mask)
					#print (cleaned)	

					if y.count(-99)!=0:
						indices_to_delete = [i for i, val in enumerate(y) if val == -99]
						for index in sorted(indices_to_delete, reverse=True):
							del y[index]
							del X[index]
				
					y = np.array(y)
					X = np.array(X)
			
						
					if len(y)<=5:
						#print (stn,key,'y is small')
						if stn not in break_summary[0]:
							break_summary[0][stn]={}
				#		#print(stn)
						
						break_summary[0][stn][key]={}
						break_summary[0][stn][key]['B_1']={}
						break_summary[0][stn][key]['B_1']['start'] = (gf_data_raw['StartDate'][0],0)	
						break_summary[0][stn][key]['B_1']['end'] = (gf_data_raw['EndDate'][gf_data_raw['Anemometer'].index('Vaisala')-1],0)
						break_summary[0][stn][key]['B_1']['Anemometer'] = gf_data_raw['Anemometer'][gf_data_raw['Anemometer'].index('Vaisala')-1]
						break_summary[0][stn][key]['B_1']['slope']= 0
						break_summary[0][stn][key]['B_1']['inter']= 0
						break_summary[0][stn][key]['B_1']['GF_diff']= 0
						break_summary[0][stn][key]['B_1']['davies']=np.nan
						break_summary[0][stn][key]['B_1']['significant']=False	
						break_summary[0][stn][key]['B_1']['GF_series']=[]
						Belfort_no_miss = [x for x in gf_data_raw['gust_factor'][key][0:Vaisala_start-2] if x != -99]
						
						for GF_yr in range (1,Vaisala_start):
							break_summary[0][stn][key]['B_1']['GF_series'].append(np.mean(Belfort_no_miss))
				
						break_summary[0][stn][key]['SeriesType']='noTrendInRawSeries'
						break_summary[0][stn][key]['AlteredBreaks']='none'
						break_summary[0][stn][key]['Modified']='none'
						break_summary[0][stn][key]['Raw_GF_series']=gf_data_raw['gust_factor'][key]
						
						#input('look small Y')
						
						continue

					
					ms = piecewise_regression.ModelSelection(X, y, max_breakpoints=total_breaks, verbose = False)
					#ms3 = piecewise_regression.ModelSelection(X, y, max_breakpoints=3, verbose = False)
						
					
					#let piecewise_regression fit models with multiple break points nand choose the best number based on BIC	
					min_bic = 99999
					for bic_val in range(len(ms.model_summaries)):
						if ms.model_summaries[bic_val]['converged']:
							if ms.model_summaries[bic_val]['bic']<min_bic:
								min_bic = ms.model_summaries[bic_val]['bic']
								num_breaks = bic_val
					
					# now fit the model for the selected number of breaks using the pwlf package 			
					my_pwlf = pwlf.PiecewiseLinFit(X,y)
					breaks = my_pwlf.fit(num_breaks+1)
					
					
					short_break = 'none'
					#print (breaks,num_breaks,'FIRST BREAKS',my_pwlf.slopes,len(breaks))
					#print(y_hold)

					#if the len of breaks is 2 it is just the start and end of series

					if len(breaks)>2:

						close_breaks = find_indices_within_range(breaks,min_break_duration)
						#print(close_breaks)

						breaks = breaks.tolist()
						breaks_new = copy.deepcopy(breaks)
						
						del_list = []	
						

						for index in sorted(close_breaks, reverse=True):
							#print(index,num_breaks,close_breaks)	
							if index[0]==0:
								if not(breaks[index[1]] in del_list or breaks[index[0]] in del_list):
									del_list.append(breaks[index[1]])
									#### need to delete the value not the index so something like
									###del breaks_new[breaks_new.index(breaks[2])]
									###del breaks_new[breaks_new.index(breaks[3])]
									### del breaks_new[breaks[index[1]].index]
									del breaks_new[breaks_new.index(breaks[index[1]])]
									num_breaks = num_breaks -1
								#	print('deleting 1',num_breaks)
									if short_break != 'none':
										short_break = short_break+'_and_start'
									else:
										short_break = 'at_start'
								#	print('del start',index,index[1])
								
							else:
								y0_idx = int(breaks[index[0]])-1
								y1_idx = int(breaks[index[1]])-1
								y1_next_idx = int(breaks[index[1]])   #(+1-1)
								y0_next_idx = int(breaks[index[0]])   #(+1-1)
	
								if int(index[1])-int(index[0])==1:
								
									if y1_next_idx<len(X)-1:
										#print(y[y0_idx],y[y1_idx],y[y1_next_idx],y[y1_idx]-y[y0_idx],(y[y1_next_idx]-y[y1_idx]))
										if y[y1_idx]-y[y0_idx]>(y[y1_next_idx]-y[y1_idx]):
										
										#	print (breaks[index[1]],breaks[index[0]],del_list)	
											if not(breaks[index[1]] in del_list or breaks[index[0]] in del_list):
												del_list.append(breaks[index[1]])
												del breaks_new[breaks_new.index(breaks[index[1]])]
												
												num_breaks = num_breaks -1
										#		print('deleting 2',num_breaks)
										#		print(short_break)
												if short_break == 'none': 
													short_break = 'at_middle'
												elif short_break != 'at_middle':
													short_break = short_break +'_and_middle'
													
											#	print('del one',index[1],breaks[index[1]])
										else:
										#	print (breaks[index[1]],breaks[index[0]],del_list)
											if not (breaks[index[1]] in del_list or breaks[index[0]] in del_list):
												del_list.append(breaks[index[0]])
												del breaks_new[breaks_new.index(breaks[index[0]])]
											#	print ('del zero',index[0],short_break)
												num_breaks = num_breaks -1
										#		print('deleting 3',num_breaks,breaks_new)
												if short_break == 'none': 
													short_break = 'at_middle'
												elif short_break != 'at_middle':
													short_break = short_break +'_and_middle'
									else:
									#	print (breaks[index[1]],breaks[index[0]],del_list)	
										if not (breaks[index[1]] in del_list or breaks[index[0]] in del_list):
											del_list.append(breaks[index[0]])
											del breaks_new[breaks_new.index(breaks[index[0]])]
											num_breaks = num_breaks -1
									#		print('deleting 4',num_breaks)
											if short_break == 'none':
												short_break = 'at_end'
												
									#		print('del else',index,close_breaks,index[0],del_list,breaks,breaks_new)

								
								else:
								#	print(y[y0_next_idx],y[y0_idx],y[y1_idx],y[y0_next_idx],y0_next_idx,y0_idx,y1_idx)
									if y[y0_next_idx]-y[y0_idx]>y[y1_idx]-y[y0_next_idx]:
								#		print (breaks[index[1]],breaks[index[0]],del_list)
										if not (breaks[index[1]] in del_list or breaks[index[0]] in del_list):
											del_list.append(breaks[index[1]])
											del breaks_new[breaks_new.index(breaks[index[1]])]
											num_breaks = num_breaks -1
								#			print('deleting 5',num_breaks)
								#			print(short_break)
											if short_break == 'none':
												short_break = 'at_middle'
											elif short_break != 'at_middle':
												short_break = short_break +'_and_middle'
									#		print ('del one two',index[0],index,close_breaks)
								#			print (y[y0_next_idx],y[y0_idx],y[y1_idx],y[y0_next_idx],y0_next_idx,y0_idx)
								#			print (y)
											#input('look del one two')
									else:
								#		print (breaks[index[1]],breaks[index[0]],del_list,num_breaks)
										if not(breaks[index[1]] in del_list or breaks[index[0]] in del_list):
											del_list.append(breaks[index[0]])
											del breaks_new[breaks_new.index(breaks[index[0]])]
											num_breaks = num_breaks -1
								#			print('deleting 6',num_breaks)
								#			print(breaks_new)
								#			print(short_break,num_breaks)
											if short_break == 'none': 
												short_break = 'at_middle'
											elif short_break != 'at_middle':
												short_break = short_break +'_and_middle'
								#			print ('del zero two',index[0])
											#input('look del zero two')
									

						breaks = [int(val) for val in breaks]
						breaks_new = [int(val) for val in breaks_new]
						breaks = np.array(breaks)
						breaks_new = np.array(breaks_new)
				
				#		print (breaks_new,num_breaks,'SECOND BREAKS',breaks)
				#		print (y)
				#		print (y_hold)
						#input('look')
						#continue
			
	
					if len(breaks)>2: breaks = breaks_new  ### added 4/8
					
					my_pwlf = pwlf.PiecewiseLinFit(X, y)
					my_pwlf.fit_with_breaks(breaks)

					added_end = []
					segment_end = []
					for brks in range (len(breaks)):
						if brks==0:
							segment_end.append((breaks[brks],breaks[brks+1]))
							added_end.append(breaks[brks])
							added_end.append(breaks[brks])    ## was breaks[brks]+1
						elif brks < len(breaks) and brks >1:
							segment_end.append((breaks[brks-1],breaks[brks]))   ## was breaks[brks-1]+1
							added_end.append(breaks[brks-1])  ## was breaks[brks-1]+1
							added_end.append(breaks[brks])
				#	print(added_end)
				#	print(segment_end)
					my_pwlf = pwlf.PiecewiseLinFit(X, y)
					added_end = np.array(added_end)
					my_pwlf.fit_with_breaks(added_end)
					
					slp_list = []
					inter_list = []
					for slopes in range(len(my_pwlf.slopes)):
						if slopes == 0:
							slp_list.append(my_pwlf.slopes[0])
							inter_list.append(my_pwlf.intercepts[0])
						if len(my_pwlf.slopes)==2:
							if slopes == 1:
								slp_list.append(my_pwlf.slopes[slopes])
								inter_list.append(my_pwlf.intercepts[slopes])
						elif slopes%2==1:
							continue
						else:
							if slopes!=0:
								slp_list.append(my_pwlf.slopes[slopes])
								inter_list.append(my_pwlf.intercepts[slopes])
	
				#	if ms.model_summaries[bic_val]['davies']<0.05:
						
					if stn not in break_summary[num_breaks]:
						break_summary[num_breaks][stn]={}
				#		#print(stn)
						
					break_summary[num_breaks][stn][key]={}
					#print(key)
					#print(break_summary[num_breaks].keys())
					#print(break_summary[num_breaks][stn].keys())
						
						
					add1 = 0	
					add2 = 0


					break_summary[num_breaks][stn][key]['Raw_GF_series']=gf_data_raw['gust_factor'][key]

					#print('BEFORE THIS LOOP')

					#print(key,add1,stn,breaks)
					#input('look')
					cnt = 0
					#for brk in range(1,len(breaks)):
					
					start_check = []
					end_check = []
					for brk in range(0,len(segment_end)):
					#	if brk%2==1:continue
						cnt = cnt +1
					#	print(brk,add1)	
						modified_break = 'none'


						#print(stn,key,brk,add1,add2,breaks,gf_data_raw['Anemometer'],round(breaks[brk-1])+add1,round(breaks[brk]+add2))
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]= {}
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['start'] = (gf_data_raw['StartDate'][0+round(segment_end[brk][0])],added_end[brk])
				#		print('setting start 1',gf_data_raw['StartDate'],gf_data_raw['Anemometer'].index('Vaisala'),round(segment_end[brk][0]),gf_data_raw['StartDate'][0+round(segment_end[brk][0])],cnt,brk)
				#		print(break_summary[num_breaks][stn][key]['B_'+str(cnt)]['start'],'SUMMARY for ',cnt)
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['end'] = (gf_data_raw['EndDate'][0+round(segment_end[brk][1])],added_end[brk])
						
						
				#		print(key,cnt,break_summary[num_breaks][stn][key]['B_'+str(cnt)]['start'],break_summary[num_breaks][stn][key]['B_'+str(cnt)]['end'],round(segment_end[brk][0]))
						#input('look')
						
							
						#if brk==1 or len(segment_end)==1:   #### added len(segment_end)==1 4/17
						
						#	break_summary[num_breaks][stn][key][cnt]['start'] = (gf_data_raw['StartDate'][gf_data_raw['Anemometer'].index('Vaisala')],added_end[brk])
						#	print('setting start 2',cnt,brk,gf_data_raw['StartDate'][gf_data_raw['Anemometer'].index('Vaisala')])
							#print(key,cnt,break_summary[num_breaks][stn][key][cnt]['start'],break_summary[num_breaks][stn][key][cnt]['end'],round(segment_end[brk][0]))
							#input('I am here')
													
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['Anemometer'] = gf_data_raw['Anemometer'][0+round(segment_end[brk][0])]
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['Height'] = gf_data_raw['Height'][0+round(segment_end[brk][0])]
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['slope']= slp_list[brk]
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['inter']= inter_list[brk]
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_diff']= ((slp_list[brk]*segment_end[brk][1])+ inter_list[brk])-((slp_list[brk]*segment_end[brk][0])+ inter_list[brk])
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['davies']=ms.model_summaries[bic_val]['davies']
							
						#if brk == 1:add1 = 0
							
						if abs(break_summary[num_breaks][stn][key]['B_'+str(cnt)]['slope']) > 0.003 and break_summary[num_breaks][stn][key]['B_'+str(cnt)]['davies']<=0.05 and abs(break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_diff'])>0.02:    ### changed 5/1
							break_summary[num_breaks][stn][key]['B_'+str(cnt)]['significant']= True
						else:
							break_summary[num_breaks][stn][key]['B_'+str(cnt)]['significant']= False
			
						
#old_y = gf_data_raw['gust_factor'][key][Vaisala_start:]
#>>> old_y = np.array(old_y)
#>>> valid_mask = old_y != -99
#>>> cleaned = old_y[valid_mask]
#>>> original_indices = np.where(valid_mask)[0]
#>>> original_indices
###array([ 0,  1,  2,  3,  4,  5,  7,  8,  9, 10, 11, 12, 14])

					#	print (original_indices, 'ORIGINAL')
						
						if round(segment_end[brk][0]) < round(segment_end[brk][1]):   #was segment_end[brk][0])+1
						
							if brk ==0:
								new_index = original_indices[(original_indices>=round(segment_end[brk][0])-1) & (original_indices<=round(segment_end[brk][1]))]  #was (original_indices>=round(segment_end[brk][0])-1) & (original_indices<=round(segment_end[brk][1])-1)]
							else:
								new_index = original_indices[(original_indices>=round(segment_end[brk][0])+1) & (original_indices<=round(segment_end[brk][1]))]
								
							list_check =  np.array(y_hold)[new_index]
							
							if brk == 0:
								start_check.append(gf_data_raw['StartDate'][0])
							else:
								start_check.append(gf_data_raw['StartDate'][0+round(segment_end[brk][0])+1])
								#print('start_check in teh else', gf_data_raw['StartDate'],gf_data_raw['Anemometer'].index('Vaisala'),round(segment_end[brk][0]),brk,cnt)
							end_check.append(gf_data_raw['EndDate'][0+round(segment_end[brk][1])])
							
							#print (brk,list_check,'WEIRD',round(segment_end[brk][0])+1,round(segment_end[brk][1]),segment_end[brk],original_indices)
							#print(start_check, end_check)
							#print (new_index)
							#print(y_hold)
								
							if len(list_check)>1:
								list_max = max(list_check)
								#list_min = min(list_check)
								list_std = np.std(list_check)
								
								#list_truncated = list_check.tolist()
								#list_truncated.remove(max(list_check))
								#list_std = np.std(np.array(list_truncated))
								
								list_avg_truncated = (sum(list_check)-max(list_check))/(len(list_check)-1)

								#print(max(list_check),np.std(list_check),(sum(list_check)-max(list_check))/(len(list_check)-1))

								if list_max> list_avg_truncated+2*list_std and break_summary[num_breaks][stn][key]['B_'+str(cnt)]['significant']: 
									break_summary[num_breaks][stn][key]['B_'+str(cnt)]['significant']= False
									modified_break = 'segment with outlier'
								#	print ('MAKING NOT SIGNIFICANT TOO HIGH',stn,key)				
								#	input('look')

								#list_truncated = list_check.tolist()
								#list_truncated.remove(min(list_check))
								#list_std = np.std(np.array(list_truncated))
								
								#list_avg_truncated = (sum(list_check)-min(list_check))/(len(list_check)-1)

								#print(min(list_check),np.std(list_check),(sum(list_check)-min(list_check))/(len(list_check)-1),modified_break)

								#if list_min < list_avg_truncated-2*list_std and break_summary[num_breaks][stn][key][cnt]['significant']: 
								#	break_summary[num_breaks][stn][key][cnt]['significant']= False
								#	modified_break = 'segment with outlier'
								#	print ('MAKING NOT SIGNIFICANT TOO LOW',stn,key)				
								#	input('look')		
								
							
							else:
								break_summary[num_breaks][stn][key]['B_'+str(cnt)]['significant']= False
								modified_break = 'segment with weird missing data'
							#	print ('WEIRD SEGMENTS MISSING DATA',stn,key,gf_data_raw['gust_factor'][key],y_hold)
								#input('look weird segment')	
			
						break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_series']=[]
							
						if break_summary[num_breaks][stn][key]['B_'+str(cnt)]['significant']:
							
							if len(list_check)>1:
							#	print(brk,round(segment_end[brk][0]),round(segment_end[brk][1]),list_check)
							#	print(list_check,'LIST CHECK')
								break_summary[num_breaks][stn][key]['B_'+str(cnt)]['start'] = (start_check[cnt-1],added_end[brk])
							#	print(break_summary[num_breaks][stn][key]['B_'+str(cnt)]['start'],'SUMMARY for ',cnt)
							#	print('setting start 3a',start_check,cnt,brk)
								break_summary[num_breaks][stn][key]['B_'+str(cnt)]['end'] = (end_check[cnt-1],added_end[brk])
							
							
							if cnt ==1:
								for GF_yr in range (round(segment_end[brk][0]),round(segment_end[brk][1])+1):
									break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_series'].append((slp_list[brk]*GF_yr)+ inter_list[brk])
								#print ('Art sign',round(breaks[brk-1]+add1),round(breaks[brk]+add2+1),GF_yr,y,y[round(breaks[brk-1]+add1):round(breaks[brk]+add2)])
								#if stn == 'KMSL' and key =='NNE': 
								#	print(my_pwlf.slopes[brk-1]*GF_yr+ my_pwlf.intercepts[brk-1])
								#	print(X[GF_yr-1],y[GF_yr-1])
							else:
								for GF_yr in range (round(segment_end[brk][0]),round(segment_end[brk][1])):
									break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_series'].append((slp_list[brk]*GF_yr)+ inter_list[brk])
						else:
							#list_check = y[round(segment_end[brk][0]):round(segment_end[brk][1])+1]
							if len(list_check)>1:
						
								if cnt==1:
									for GF_yr in range (round(segment_end[brk][0]),round(segment_end[brk][1])+1):
									#print ('Art1',round(segment_end[brk][0]),round(segment_end[brk][1]),y[round(segment_end[brk][0]-1):round(segment_end[brk][1])])
									#break_summary[num_breaks][stn][key][cnt]['GF_series'].append(np.mean(y[round(segment_end[brk][0]-1):round(segment_end[brk][1])]))
										break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_series'].append(np.mean(list_check))
								
								else:
									for GF_yr in range (round(segment_end[brk][0]),round(segment_end[brk][1])):
										break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_series'].append(np.mean(list_check))
										
											
							#	print(brk,round(segment_end[brk][0]),round(segment_end[brk][1]),list_check)
							#	print(list_check,'LIST CHECK')
								break_summary[num_breaks][stn][key]['B_'+str(cnt)]['start'] = (start_check[cnt-1],added_end[brk])
							#	print(break_summary[num_breaks][stn][key]['B_'+str(cnt)]['start'],'SUMMARY for ',cnt)
							#	print('setting start 3',start_check,cnt,brk)
								break_summary[num_breaks][stn][key]['B_'+str(cnt)]['end'] = (end_check[cnt-1],added_end[brk])
									#break_summary[num_breaks][stn][key][cnt]['end'] = (gf_data_raw['EndDate'][gf_data_raw['Anemometer'].index('Vaisala')+round(segment_end[brk][1])],added_end[brk])
									#break_summary[num_breaks][stn][key][cnt]['Anemometer'] = gf_data_raw['Anemometer'][gf_data_raw['Anemometer'].index('Vaisala')+round(segment_end[brk][0])]
									#break_summary[num_breaks][stn][key][cnt]['Height'] = gf_data_raw['Height'][gf_data_raw['Anemometer'].index('Vaisala')+round(segment_end[brk][0])]
									
									
									

							else:
								y_hold  = [np.nan if x == -99.0 else x for x in y_hold]
								
								if cnt==1:
									for GF_yr in range (round(segment_end[brk][0]),round(segment_end[brk][1])+1): 
										#print ('Art2',round(segment_end[brk][0]),round(segment_end[brk][1]),y[round(segment_end[brk][0]-1):round(segment_end[brk][1])])
										break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_series'].append(np.nanmean(y_hold[round(segment_end[brk][0]-1):round(segment_end[brk][1])]))
								else:
									for GF_yr in range (round(segment_end[brk][0]),round(segment_end[brk][1])):
										break_summary[num_breaks][stn][key]['B_'+str(cnt)]['GF_series'].append(np.nanmean(y_hold[round(segment_end[brk][0]-1):round(segment_end[brk][1])])) 
								
						#if add2 == 1:
						#	add1 = 2    ##changed from 2
						#else:
						#	add1 = 1
						
					#print(break_summary[num_breaks][stn][key])
					#input('look')	
					#break_summary[num_breaks][stn][key][0]= {}
					
					#if 	gf_data_raw['StartDate'][0]!= gf_data_raw['StartDate'][Vaisala_start-1]:
						
						
					#	break_summary[num_breaks][stn][key][0]['start'] = (gf_data_raw['StartDate'][0],0)	
					#	break_summary[num_breaks][stn][key][0]['end'] = (gf_data_raw['EndDate'][gf_data_raw['Anemometer'].index('Vaisala')-1],0)
					#	break_summary[num_breaks][stn][key][0]['Anemometer'] = gf_data_raw['Anemometer'][gf_data_raw['Anemometer'].index('Vaisala')-1]
					##	break_summary[num_breaks][stn][key][0]['slope']= 0
					#	break_summary[num_breaks][stn][key][0]['inter']= 0
					#	break_summary[num_breaks][stn][key][0]['GF_diff']= 0
					#	break_summary[num_breaks][stn][key][0]['davies']=np.nan
					#	break_summary[num_breaks][stn][key][0]['significant']=False	
					#	break_summary[num_breaks][stn][key][0]['GF_series']=[]
					#else:
					#	break_summary[num_breaks][stn][key][0]['start'] = ('none',0)
					#	break_summary[num_breaks][stn][key][0]['end'] = ('none',0)
					#	break_summary[num_breaks][stn][key][0]['Anemometer'] = 'Belfort'
					#	break_summary[num_breaks][stn][key][0]['Height'] = 'none'
					#	break_summary[num_breaks][stn][key][0]['slope']= 0
					#	break_summary[num_breaks][stn][key][0]['inter']= 0
					#	break_summary[num_breaks][stn][key][0]['GF_diff']= 0
					#	break_summary[num_breaks][stn][key][0]['davies']=np.nan
					#	break_summary[num_breaks][stn][key][0]['significant']=False
					#	break_summary[num_breaks][stn][key][0]['GF_series']=[]
						
					
					#Belfort_no_miss = [x for x in gf_data_raw['gust_factor'][key][0:Vaisala_start-2] if x != -99]
					#for GF_yr in range (1,Vaisala_start):
					#		break_summary[num_breaks][stn][key][0]['GF_series'].append(np.mean(Belfort_no_miss))
					
				#	print(break_summary[num_breaks][stn][key][1]['GF_series'])
					if break_summary[num_breaks][stn][key]['B_1']['significant']:
						break_summary[num_breaks][stn][key]['B_1']['GF_series'].insert(0,my_pwlf.intercepts[0])
					#	print ('inserted ',Vaisala_start-1,break_summary[num_breaks][stn][key][1]['GF_series'])
					#	input('look')
					else:
						if num_breaks !=0:
							break_summary[num_breaks][stn][key]['B_1']['GF_series'].insert(0,break_summary[num_breaks][stn][key]['B_1']['GF_series'][0])
							#print('inserted brk-1 GF[0]')
						else:
							break_summary[num_breaks][stn][key]['B_1']['GF_series'].insert(0,break_summary[num_breaks][stn][key]['B_1']['GF_series'][0])
							#print('inserted brk GF[0]')

						#if num_breaks == 1 and breaks[-1]-breaks[-2]>2:
						#	print ('THIS IS A SIGNIFICANT CHANGE')
						#elif num_breaks ==0:
						#	print ('SINGLE SLOPE THROUGH TIME')
						#elif num_breaks >1:
						#	print('MORE THAN ONE BREAK')
				
				#	print('RAW',gf_data_raw['gust_factor'][key])
				#	print (break_summary[num_breaks][stn][key])
				#	if short_breaks !='none': input('got a slope, but altered breaks')	
				
					break_summary[num_breaks][stn][key]['AlteredBreaks']=short_break
					break_summary[num_breaks][stn][key]['SeriesType']='trendInRawSeries'
					break_summary[num_breaks][stn][key]['Modified']=modified_break
				else:
				#	print('no sig at start',num_breaks)
				
					num_breaks = 0
			
					if stn not in break_summary[0]:
						break_summary[num_breaks][stn]={}
					break_summary[num_breaks][stn][key]={}
					
					#break_summary[num_breaks][stn][key][0]= {}
					break_summary[num_breaks][stn][key]['B_1']= {}
					
					if 	gf_data_raw['StartDate'][0]	!= gf_data_raw['StartDate'][Vaisala_start-1]:
						
						#break_summary[num_breaks][stn][key][0]['start'] = (gf_data_raw['StartDate'][0],0)
						#break_summary[num_breaks][stn][key][0]['end'] = (gf_data_raw['EndDate'][gf_data_raw['Anemometer'].index('Vaisala')-1],0)
						#break_summary[num_breaks][stn][key][0]['Anemometer'] = gf_data_raw['Anemometer'][gf_data_raw['Anemometer'].index('Vaisala')-1]
						#break_summary[num_breaks][stn][key][0]['Height'] = gf_data_raw['Height'][gf_data_raw['Anemometer'].index('Vaisala')-1]
						#break_summary[num_breaks][stn][key][0]['slope']= 0
						#break_summary[num_breaks][stn][key][0]['inter']= 0
						#break_summary[num_breaks][stn][key][0]['GF_diff']= 0
						#break_summary[num_breaks][stn][key][0]['davies']=np.nan
						#break_summary[num_breaks][stn][key][0]['significant']=False
						#break_summary[num_breaks][stn][key][0]['GF_series']=[]
						continue
					else:
						break_summary[num_breaks][stn][key]['B_1']['start'] = ('none',0)
						break_summary[num_breaks][stn][key]['B_1']['end'] = ('none',0)
						break_summary[num_breaks][stn][key]['B_1']['Anemometer'] = 'Belfort'
						break_summary[num_breaks][stn][key]['B_1']['Height'] = 'none'
						break_summary[num_breaks][stn][key]['B_1']['slope']= 0
						break_summary[num_breaks][stn][key]['B_1']['inter']= 0
						break_summary[num_breaks][stn][key]['B_1']['GF_diff']= 0
						break_summary[num_breaks][stn][key]['B_1']['davies']=np.nan
						break_summary[num_breaks][stn][key]['B_1']['significant']=False
						break_summary[num_breaks][stn][key]['B_1']['GF_series']=[]
						
						
					#break_summary[num_breaks][stn][key][1]= {}
					#print(stn,key,'BBEFORE YEAR 1 START')
					#break_summary[num_breaks][stn][key][1]['start'] = (gf_data_raw['StartDate'][gf_data_raw['Anemometer'].index('Vaisala')],1)
					
					#break_summary[num_breaks][stn][key][1]['end'] = (gf_data_raw['EndDate'][-1],1)
						
					#break_summary[num_breaks][stn][key][1]['Anemometer'] = gf_data_raw['Anemometer'][gf_data_raw['Anemometer'].index('Vaisala')]
						
					#break_summary[num_breaks][stn][key][1]['Height'] = gf_data_raw['Height'][gf_data_raw['Anemometer'].index('Vaisala')]
						
					#break_summary[num_breaks][stn][key][1]['slope']= 0
						
					#break_summary[num_breaks][stn][key][1]['inter']= 0
						
					#break_summary[num_breaks][stn][key][1]['GF_diff']= 0

					#break_summary[num_breaks][stn][key][1]['davies']=np.nan
						
					#break_summary[num_breaks][stn][key][1]['significant']=False	
					
					#break_summary[num_breaks][stn][key][1]['GF_series']=[]

					
					Belfort_no_miss = [x for x in gf_data_raw['gust_factor'][key][0:Vaisala_start-2] if x != -99]
					Vaisala_no_miss = [x for x in gf_data_raw['gust_factor'][key][Vaisala_start:] if x != -99]
						
					for GF_yr in range (1,Vaisala_start):
						break_summary[num_breaks][stn][key]['B_1']['GF_series'].append(np.mean(Belfort_no_miss))
					
				#	for GF_yr in range (Vaisala_start-1,len(gf_data_raw['EndDate'])):
				#		break_summary[num_breaks][stn][key][1]['GF_series'].append(np.mean(Vaisala_no_miss))	
					
				#	print (break_summary[num_breaks][stn][key])
				#	input('never a trend')	
				
					break_summary[num_breaks][stn][key]['SeriesType']='noTrendInRawSeries'
					break_summary[num_breaks][stn][key]['AlteredBreaks']='none'
					break_summary[num_breaks][stn][key]['Modified']='none'
					break_summary[num_breaks][stn][key]['Raw_GF_series']=gf_data_raw['gust_factor'][key]

					
					num_stn_flat = num_stn_flat + 1

			
			
				
					#input('no trend at start')
				#print(stn,key,Vaisala_start)
				#print (break_summary[num_breaks][stn][key])
		#		input('look')
			

with open("break_summary_belfort.json", "w") as outfile:
    json.dump(break_summary, outfile, default=str)
    
print (num_stn_trend,' TREND ',num_stn_flat, 'FLAT')
print (dir_cnts)
