import json
import numpy as np
import copy
import piecewise_regression
import pwlf


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


with open('break_summary_belfort.json', 'r') as file:
	break_data = json.load(file)	


total_breaks = 4
min_break_duration = 2

for key in break_data.keys():
	for stn in break_data[key].keys():
		#if stn !='KFTW':continue
		for dir in break_data[key][stn].keys():
			#if break_data[key][stn][dir]['Modified']!='none':
			
			for brk in range (0,6): #(1,int(key)+1):
				if 'B_'+str(brk) not in break_data[key][stn][dir].keys():continue
				#if 'B_'+str(brk+1) not in break_data[key][stn][dir].keys():continue

				if brk!=0 and 'B_'+str(brk+1) in break_data[key][stn][dir].keys():
					if int(break_data[key][stn][dir]['B_'+str(brk)]['end'][0][0:4])+1 != int(break_data[key][stn][dir]['B_'+str(brk+1)]['start'][0][0:4]):
						print ('DATE MISMATCH', key,dir)
						print()
						print()
				else:
					if break_data[key][stn][dir]['B_'+str(brk)]['end'][0][0:4]!='none':
						if 'B_'+str(brk+1) in break_data[key][stn][dir].keys():
							if int(break_data[key][stn][dir]['B_'+str(brk)]['end'][0][0:4]) != int(break_data[key][stn][dir]['B_'+str(brk+1)]['start'][0][0:4]):
								print ('DATE MISMATCH Belfort switch', stn,key,dir)
								print()
								print()
							#print (break_data[key][stn][dir][str(brk)])
							#print (break_data[key][stn][dir][str(brk+1)])
							#print (break_data[key][stn][dir]["Raw_GF_series"])
					
					if break_data[key][stn][dir]['B_'+str(brk)]['end'][0][0:4]=='none':
						if not break_data[key][stn][dir]['B_1']['significant']: 
							if np.mean(break_data[key][stn][dir]['B_1']["GF_series"])!=np.mean(break_data[key][stn][dir]["Raw_GF_series"][0:len(break_data[key][stn][dir]['B_1'])]):
								print ('MEAN MISMATCH NO Belfort switch', stn,key,dir,break_data[key][stn][dir]['B_1']['significant'])
								
								num_yrs = len(break_data[key][stn][dir]['B_1']["GF_series"])
								break_data[key][stn][dir]['B_1']["GF_series"]=[]
								for val in range (num_yrs):
									break_data[key][stn][dir]['B_1']["GF_series"].append(np.mean(break_data[key][stn][dir]["Raw_GF_series"][0:num_yrs]))
								
								print (break_data[key][stn][dir],"Changed series associated with break1")
						else:
							print (' Belfort switch WITH TRENDs', stn,key,dir)
							print (break_data[key][stn][dir])
							

### flag segments with a big undetected jump at the end
				if len(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"])>=4:
					if abs(np.mean(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"][-1:])- np.mean(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"][:-2]))>0.01 and not break_data[key][stn][dir]['B_'+str(brk)]['significant']:
						print('JUMP AT END',break_data[key][stn][dir]['B_'+str(brk)]["GF_series"],break_data[key][stn][dir]['B_'+str(brk)]['significant'],brk)

### flag neighboring segments that are flat and don't change by much from one segment to the next
				if brk>1:
				
					print(break_data[key][stn][dir],'BEFORE')
					
					if abs(break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"][0]- break_data[key][stn][dir]['B_'+str(brk)]["GF_series"][-1])<0.02:  #  and not break_data[key][stn][dir]['B_'+str(brk-1)]['significant'] and not break_data[key][stn][dir]['B_'+str(brk)]['significant']:
						print('NO REAL DIFFERENCE IN GF',break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"],break_data[key][stn][dir]['B_'+str(brk)]["GF_series"],brk)
				
						skip_yrs = 1
						for val in range(1,brk-1):
							print('SUMMING SKIPS',val,skip_yrs,len(break_data[key][stn][dir]['B_'+str(val)]["GF_series"]))
							skip_yrs = skip_yrs + len(break_data[key][stn][dir]['B_'+str(val)]["GF_series"])
							
						print(break_data[key][stn][dir])
						
						new_GF_list = []

						print ('SKIPS',skip_yrs,len(break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"]),len(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"]),skip_yrs+len(break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"])+len(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"]))
						for val in break_data[key][stn][dir]["Raw_GF_series"][skip_yrs:skip_yrs+len(break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"])+len(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"])-1]:
							new_GF_list.append(val)
						
						
						num_yrs = len(break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"]) + len(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"])
						
						break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"]=[]
						
						for val in range (num_yrs):
							break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"].append(np.mean(new_GF_list))
						
						print ('I AM GF LIST NEW', new_GF_list,len(break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"]),len(break_data[key][stn][dir]['B_'+str(brk)]["GF_series"]),num_yrs)
						print ("intermed ", break_data[key][stn][dir]['B_'+str(brk-1)]["GF_series"])
						break_data[key][stn][dir]['B_'+str(brk-1)]['end']=break_data[key][stn][dir]['B_'+str(brk)]['end']
						
						for val in range (brk,find_maximum_numeric_value(break_data[key][stn][dir].keys())+1):
							if val < find_maximum_numeric_value(break_data[key][stn][dir].keys()):
								break_data[key][stn][dir]['B_'+str(val)]=break_data[key][stn][dir]['B_'+str(val+1)]
								if break_data[key][stn][dir]['B_'+str(val)]["GF_series"][0]==break_data[key][stn][dir]['B_'+str(val)]["GF_series"][-1]:
										break_data[key][stn][dir]['B_'+str(val)]['significant']=False
							else:
								del break_data[key][stn][dir]['B_'+str(val)]
						
						print(stn,dir)
						print(break_data[key][stn][dir],'AFTER')
						#input('look)')

with open("break_summary_smoothed_belfort_slopes_too.json", "w") as outfile:
	json.dump(break_data, outfile, default=str)
    
#print (num_stn_trend,' TREND ',num_stn_flat, 'FLAT')
#print (dir_cnts)

				
