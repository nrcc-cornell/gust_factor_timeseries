This is a multiple-step process that involves running a set of programs

1) obtain one-minute data from NCEI archive by running the program runDownloadASOS_GF.py.  The program requires an input file that lists the stations to be downloaded.
  The file all_stns.txt includes just two stations in PA as a test.  The file stn_inventory.txt includes all 1-minute ASOS stations.  The output files are written to a  ../data/asos-onemin/ directory, 
  with a separate folder created for each state that is processed. Within each state folder is a folder for each station. Depending on the server load at NCEI, this program can take a while to run.
  
2) Convert the downloaded csv files to hdf5 format by running the program runCreateHdfFileFromOneMinData_GF.py .
   The output files are written to a  ../data/asos-onemin/ directory, with a separate folder created for each state that is processed. Within each state folder is a single file for each station.

3) Perform qc on the data in the hdf5 files by running the program runCreateHdfFileFromOneMinData_withErrorCheck_GF.py .

   The program runCreateHdfFileFromOneMinData_withErrorCheck_GF.py creates two files manual_data_dates.txt and manual_data_checks.txt.  These files contain data points that the automated qc 
   deemed suspect and should be screened and manually adjusted if necessary.  An example line from the file manual_data_dates.txt file is:
   	
   	PA,KABE,wind,b'200907071342',31.0,88.0,0.0,obs,101
   	
   A sample line from the manual_data_checks.txt file is:
   
   PA,KABE,wind,b'200907071342',31.0,88.0,0.0,obs,101,0
   
   manual_data_dates.txt contains observations that have been flagged as suspect and manual_data_checks.txt contains data from surrounding minutes to assist in manually validting the suspect observations.
   
4) Incorporate any manual corrections identified in step 3 by running the program runCreateHdfFileFromOneMinData_manualCorrections_GF.py  This requires the creation of a new file (manual_data_dates_long_POR_invalidSaved.txt)  prior to running. This file is
   simply the lines from manual_data_dates.txt that the analyst determines are erroneous.  For example if the analysts determines that the 88.0 wind speed in the above line 
   (PA,KABE,wind,b'200907071342',31.0,88.0,0.0,obs,101,0) is a bad observation, this line should be copied to manual_data_dates_long_POR_invalidSaved.txt.  The program reads through this file and sets the wind speed to missing 
   for those stations dates and times that are included.  The program creates new files with the edited missing values included in the asos-onemin-hdf-files-with-error-check-and-manual-corrections directory.
   

5) Run the program GustFactor_FileWriter.py to compute annual gust factor time series based on the manually edited file produced in step 4.  Requires the file ASOS_anem_metadata.txt

6) Run count_sign_GF_slopes_v2_thru24.py to perform piecewise regression on the annual gust factor time series produced in step 5.  Output is a large json file with adjusted gust factor time series for each station-wind direction combination. 

7) Run count_sign_GF_slopes_onlybelfort.py to perform piecewise regression on the annual gust factor time series during the period with a belfort anemometer in use at the station. Output is a large json file with adjusted gust factor time series for each station-wind direction combination.

8) Run the program Break_screener_slopest_too.py to make adjustments to the raw piecewise regression results.  This primary combines adjacent periods of time in which the difference in gust factors is small from a practical standpoint.

9) Run the program Break_screener_belfort_slopest_too.py to make adjustments to the raw piecewise regression results for the belfort era.  This primary combines adjacent periods of time in which the difference in gust factors is small from a practical standpoint.

10 Run the program smoother_json_combine.py  to combine the separe results for the Belfort an Vaisala periods into a single json output file.

Final output

The final output file is in the json format.  The file is keyed on 4-letter ICAO (station) id and then wind direction.  For each wind direction there is a key for each segment of the gust factor series. The segment keys are designated V_0,V_1, V-2, etc.  
For each segment there is a dictionary that gives the segment's starting date "start", ending date "end", Anemometer type "Anemometer", Anemometer height, "Height", the slope of the segment, "slope",
 the  y intercept of the regression "inter", gust factor difference between the current segment and the previous segment "GF_diff", the Davie's statistic for the piecewise fit "davies", 
 whether a significant slope was present (true or false) "significant", and the smoothed gust factor series for the segment "GF_series".  In addition to the segment keys, for each wind direction keys also indicate the raw gust factor series "Raw_GF_series",
 whether the breaks given by the piecewise regresson needed to be manually altered "AlteredBreaks", the general type of series ", "SeriesType", and whether the data required any manual modification "Modified".
 
 AlteredBreaks can either be "none", "at_start", "at_end", "at_middle", or a combination (e.g. at_start_and_middle). These describe cases in with the piecewise regression identifies segments that are two or less years in length. Such breaks are omitted and the piecewise regression refit.  
 The AlteredBreaks key indicates if this occurred and where in the GF series the short break segment was located. "at_start" signifies the first Vaisala segment was too short and "at_end" denotes the last segment.
 
 Modified can have three values: "none", "segment with outlier" and "segment with weird missing data".  Segments with outliers have a significant slope that is influenced by a value that exceeds the average gust factor of a segment by more than two standard deviations. 
 This value is ignored and the segment is assigned the average of the remaining gust factors.  Segments in which missing data influence the slope are treated in a similar way and identified as segments with weird missing data.
 
 SeriesType can either be "trendInRawSeries" or "noTrendInRawSeries".  When identified as "noTrendInRawSeries", the piecewise regression did not identify any beaks in the Vaisala anemometer segment.

  