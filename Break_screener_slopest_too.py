import json
import numpy as np
import copy
import piecewise_regression
import pwlf
import copy


def find_maximum_numeric_value(string_list):
	"""
	Finds the maximum numeric value in a list of strings.

	Args:
	string_list: A list of strings, some of which may represent numbers.

	Returns:
	The maximum numeric value found in the list, or None if no numeric values are present.
	"""
	numeric_values = []
	for s in string_list:
		try:
			numeric_values.append(int(s[-1]))
		except ValueError:
			pass  # Ignore non-numeric strings

	if numeric_values:
		return max(numeric_values)
	else:
		return None

def find_indices_within_range(lst, range_value):
	indices = []
    
	for i in range(len(lst)):
		for j in range(i + 1, len(lst)):
			if abs(int(lst[i]) - int(lst[j])) <= range_value:
				indices.append((i, j))
    
	return indices


with open('break_summary_thru24.json', 'r') as file:
	break_data = json.load(file)	

for iteration in range (3):
	for key in break_data.keys():
		for stn in break_data[key].keys():
			#if stn !='KFTW':continue
			for dir in break_data[key][stn].keys():
				#if break_data[key][stn][dir]['Modified']!='none':
	
				GF_hold = copy.deepcopy(break_data[key][stn][dir]["Raw_GF_series"])
				GF_hold = [np.nan if x == -99.0 else x for x in GF_hold]
						
				for brk in range(0,6):
					if 'V_'+str(brk) not in break_data[key][stn][dir].keys():continue
					#if 'V_'+str(brk+1) not in break_data[key][stn][dir].keys():continue
	
						
					if brk!=0 and 'V_'+str(brk+1) in break_data[key][stn][dir].keys():
						if int(break_data[key][stn][dir]['V_'+str(brk)]['end'][0][0:4])+1 != int(break_data[key][stn][dir]['V_'+str(brk+1)]['start'][0][0:4]):
							print ('DATE MISMATCH', key,dir)
							print()
							print()
					else:
						if break_data[key][stn][dir]['V_'+str(brk)]['end'][0][0:4]!='none':
							if 'V_'+str(brk+1) in break_data[key][stn][dir].keys():
								if int(break_data[key][stn][dir]['V_'+str(brk)]['end'][0][0:4]) != int(break_data[key][stn][dir]['V_'+str(brk+1)]['start'][0][0:4]):
									print ('DATE MISMATCH Belfort switch', stn,key,dir)
									print()
									print()
								#print (break_data[key][stn][dir][str(brk)])
								#print (break_data[key][stn][dir][str(brk+1)])
								#print (break_data[key][stn][dir]["Raw_GF_series"])
						
						if break_data[key][stn][dir]['V_'+str(brk)]['end'][0][0:4]=='none':
							if not break_data[key][stn][dir]['V_1']['significant']: 
							
								if np.nanmean(break_data[key][stn][dir]['V_1']["GF_series"])!=np.nanmean(GF_hold[0:len(break_data[key][stn][dir]['V_1'])]):
									#print ('MEAN MISMATCH NO Belfort switch', stn,key,dir,break_data[key][stn][dir]['V_1']['significant'])
									
									num_yrs = len(break_data[key][stn][dir]['V_1']["GF_series"])
									break_data[key][stn][dir]['V_1']["GF_series"]=[]
									for val in range (num_yrs):
										break_data[key][stn][dir]['V_1']["GF_series"].append(np.nanmean(GF_hold[0:num_yrs]))
									
									print (break_data[key][stn][dir],"CHANGED SERIES ASSOCIATED WITH BREAK 1")
	
	
					if brk>1:
	
						print(break_data[key][stn][dir],'BEFORE',stn,dir)
						if abs(break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"][0]- break_data[key][stn][dir]['V_'+str(brk)]["GF_series"][-1])<0.02:  #  and not break_data[key][stn][dir]['V_'+str(brk-1)]['significant'] and not break_data[key][stn][dir]['V_'+str(brk)]['significant']:
							print('NO REAL DIFFERENCE IN GF',break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"],break_data[key][stn][dir]['V_'+str(brk)]["GF_series"],brk)
					
							skip_yrs = 1
							for val in range(0,brk-1):
								print('SUMMING SKIPS',val,skip_yrs,len(break_data[key][stn][dir]['V_'+str(val)]["GF_series"]))
								skip_yrs = skip_yrs + len(break_data[key][stn][dir]['V_'+str(val)]["GF_series"])
								
							print(break_data[key][stn][dir])
							
							new_GF_list = []
	
							print ('SKIPS',skip_yrs,len(break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"]),len(break_data[key][stn][dir]['V_'+str(brk)]["GF_series"]),skip_yrs+len(break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"])+len(break_data[key][stn][dir]['V_'+str(brk)]["GF_series"]))
							for val in GF_hold[skip_yrs:skip_yrs+len(break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"])+len(break_data[key][stn][dir]['V_'+str(brk)]["GF_series"])-1]:
								new_GF_list.append(val)
							
							
							num_yrs = len(break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"]) + len(break_data[key][stn][dir]['V_'+str(brk)]["GF_series"])
							
							break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"]=[]
							
							for val in range (num_yrs):
								break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"].append(np.nanmean(new_GF_list))
							
							print ('I AM GF LIST NEW', new_GF_list,len(break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"]),len(break_data[key][stn][dir]['V_'+str(brk)]["GF_series"]),num_yrs)
							print ("intermed ", break_data[key][stn][dir]['V_'+str(brk-1)]["GF_series"])
							break_data[key][stn][dir]['V_'+str(brk-1)]['end']=break_data[key][stn][dir]['V_'+str(brk)]['end']
							
							check_sign = 0
							for val in range (brk,find_maximum_numeric_value(break_data[key][stn][dir].keys())+1):
								if val < find_maximum_numeric_value(break_data[key][stn][dir].keys()):
									if break_data[key][stn][dir]['V_'+str(val)]['significant']or break_data[key][stn][dir]['V_'+str(val+1)]['significant']:
										check_sign = 1
									break_data[key][stn][dir]['V_'+str(val)]=break_data[key][stn][dir]['V_'+str(val+1)]
									

									if break_data[key][stn][dir]['V_'+str(val)]["GF_series"][0]==break_data[key][stn][dir]['V_'+str(val)]["GF_series"][-1]:
										break_data[key][stn][dir]['V_'+str(val)]['significant']=False
									
								else:
									del break_data[key][stn][dir]['V_'+str(val)]
							
							if check_sign == 1:
								print()
								print (break_data[key][stn][dir],'AFTER',stn,dir)
								#input('look')		
	
						#print(break_data[key][stn][dir],'After2',stn,dir)
	

with open("break_summary_smoothed_thru24.json", "w") as outfile:
	json.dump(break_data, outfile, default=str)
    
#print (num_stn_trend,' TREND ',num_stn_flat, 'FLAT')
#print (dir_cnts)

				
