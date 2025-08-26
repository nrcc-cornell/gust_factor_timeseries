import json
import numpy as np

with open('break_summary_smoothed_thru24.json', 'r') as file:
	break_data_V = json.load(file)

with open('break_summary_smoothed_belfort.json', 'r') as file:
	break_data_B = json.load(file)
	
Belfort_brks = ['B_0','B_1','B_2','B_3']
for key in break_data_V.keys():
	for stn in break_data_V[key].keys():
		for dir in break_data_V[key][stn].keys():
			for bel_brk in Belfort_brks:
				for key2 in break_data_B.keys():			
					try:
						art = break_data_B[key2][stn][dir][bel_brk]

						if int(break_data_B[key2][stn][dir][bel_brk]['end'][0][0:4])>2009:
						
							if "B_1" not in break_data_V[key][stn][dir].keys():
								print ('in here',stn,dir,key,key2,bel_brk)
								break_data_V[key][stn][dir][bel_brk]={"start": ["none", 0], "end": ["none", 0], "Anemometer": "Belfort", "Height": "none", "slope": 0, "inter": 0, "GF_diff": 0, "davies": np.nan, "significant": False, "GF_series": []}
								if len(break_data_B[key2][stn][dir][bel_brk]["GF_series"])==0:print(stn,dir,'empty GF')
								if "V_0" in break_data_V[key][stn][dir].keys():del break_data_V[key][stn][dir]["V_0"]
						else:
							break_data_V[key][stn][dir][bel_brk]=break_data_B[key2][stn][dir][bel_brk]
							if "V_0" in break_data_V[key][stn][dir].keys():
								#print('deleting V0',stn,dir)
								del break_data_V[key][stn][dir]["V_0"]
								
					except:
						#if stn=='KBHM':print(key,stn,dir,'in continue')
						continue
			if "V_0" in break_data_V[key][stn][dir].keys():
				#print('deleting V0',stn,dir)
				break_data_V[key][stn][dir]["B_1"]={"start": ["none", 0], "end": ["none", 0], "Anemometer": "Belfort", "Height": "none", "slope": 0, "inter": 0, "GF_diff": 0, "davies": np.nan, "significant": False, "GF_series": []}
				del break_data_V[key][stn][dir]["V_0"]

no_key_brks = {}
for key in break_data_V.keys():
	for stn in break_data_V[key].keys():
		if stn not in no_key_brks.keys():
			no_key_brks[stn]={}
		for dir in break_data_V[key][stn].keys():
			if dir not in no_key_brks[stn].keys():
				no_key_brks[stn][dir]=break_data_V[key][stn][dir]

	
				
with open("break_summary_combined_thru24_nokey.json", "w") as outfile:
	json.dump(no_key_brks, outfile, default=str)
					
