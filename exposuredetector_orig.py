## Script to determine "exposure to coastline" for GNAF properties.
## Requires 10.1 for features incl. arcpy.da cursors & IN_MEMORY workspace for rasters
## Not original author of initial script

import os
import itertools
import time
import arcpy
import csv
import cPickle as pickle
from math import log
from arcpy import env
from arcpy.sa import * # I hate that they ask me to do that
from arcgisscripting import ExecuteError

print "\n"+("="*79)+"\nStarted\n"
starttime = time.time()

##========== check out SA ext'n
arcpy.CheckOutExtension("spatial")

##========== utility functions
def save_object(obj,fi):
    with open(fi,"w") as f:
        pickle.dump(obj,f)
    return

def load_rows(fi):
    """Attempt to load the pickled set in fi, else return a new empty set"""
    try:
        with open(fi,"r") as f:
            return pickle.load(f)
    except IOError:
        return set()

def directionfieldnames():
    ##========== Set up the fields to take the measurements
    prefs = ["South","Southwest","West","Northwest","North","Northeast","East","Southeast"]
    suffs = ["ExposedAtAll","HighlyExposed","ExposureRate"]
    return ["ExposedAtAllInAnyDirection"]+ \
           ["{0}{1}".format(i,j) for i,j in itertools.product(prefs,suffs)]

def write_header_row(outcsv):
    try:
        with open(outcsv, "r") as f:
            pass
    except IOError:
        with open(outcsv, "w") as f:
            f.write(",".join(["GNAF_PID","lat","lon"]+directionfieldnames()) + "\n")

##========== logic functions
def getExposureResults(zonalstats_results_table):#,output_fieldnames):
    """Return whether a point is exposed at all, or majorly, from each direction.
    Relies on a particular format - the output of the Zonal Statistics as Table tool.
    Assumes there are always 8 classes numbered 1-8 where 1 is South proceeding clockwise
    to 8 being Southeast."""

    # get the values of the table
    with arcpy.da.SearchCursor(zonalstats_results_table,["MAX","MAJORITY","MEAN"]) as tblcur:
        results = (r for r in tblcur)
    # results is now a generator containing tuples, so flatten it
    flat_results = itertools.chain.from_iterable(results)
    # return a tuple for easy appending to the searchcursor results
    return tuple(flat_results)

def exposure(in_property_layer,nearby_properties_layer,ocean_layer):
    print "\t...commencing raster processing..."
    try:
        cellsize = 10
        # remap objects for the reclassifications
        combined_remap = RemapValue([["NODATA", 0],
                                     [-1, 0],
                                     [0, 0],
                                     [10, 10],
                                     [11, 10],
                                     [200, 1]])
        euc_dir_remap = RemapRange([[-1, -1, 0],
                                    [-1, 22.5, 1],
                                    [22.5, 67.5, 2],
                                    [67.5, 112.5, 3],
                                    [112.5, 157.5, 4],
                                    [157.5, 202.5, 5],
                                    [202.5, 247.5, 6],
                                    [247.5, 292.5, 7],
                                    [292.5, 337.5, 8],
                                    [337.5, 360, 1]])

        # generate land_ocean raster
        arcpy.FeatureToRaster_conversion(ocean_layer, "HY_CS_CODE", "_ocean", cellsize)
        land_ocean = Reclassify("_ocean", "VALUE", combined_remap)

        # create surface from nearby property points
        arcpy.FeatureToRaster_conversion(nearby_properties_layer, offset_fieldname,
                                         "point_heights", cellsize)
        vis_surface = Reclassify("point_heights", "VALUE", combined_remap, "DATA")

        # Apparently we have to save the input FC to disk now.
        # This adds a whole second to the processing =\
        arcpy.FeatureClassToFeatureClass_conversion(in_property_layer,
                                                    r"C:\ExposureDetection\scratch", "inprop.shp")

        # Do the viewshed & find visible ocean
        viewshed = Viewshed(vis_surface, r"C:\ExposureDetection\scratch\inprop.shp")
        visible_ocean = (land_ocean == 1) & (viewshed > 0)

        # Euclidean direction & classify into zones
        euc_dir = EucDirection(in_property_layer,analysis_distance)
        euc_zones = Reclassify(euc_dir, "VALUE", euc_dir_remap, "DATA")

        # Do the stats
        print "\t...calculating statistics..."
        ZonalStatisticsAsTable(euc_zones, "VALUE", visible_ocean, "zStats")

        # Now go off and get the deets from that table.
        return getExposureResults("zStats")

    except ExecuteError as e:
        print " Raster processing failed! ".center(79,"=")
        print "Error was:", e

        # Return null values for failed rows.
        return (None, ) * 8


##========== Input parameters & other "Constants" =============================
gnaf_fc = r"C:\ExposureDetection\input\Alpha_Phase3.gdb\GNAF_Nowra_to_Mackay_1km_From_Coast_GALambert"
#hydro_fc = r"C:\ExposureDetection\input\Alpha_Phase3.gdb\Nowra_to_Mackay_Ocean_GALambert"
hydro_fc = r"C:\ExposureDetection\input\Alpha_Phase3.gdb\Nowra_to_Mackay_Ocean_GALambert_simplifiedNoHoles"

out_csv = r"C:\ExposureDetection\output\Output.csv"
processedrowset = r"C:\ExposureDetection\output\ProcessedRowSet.pickle"

env.workspace = "IN_MEMORY"
env.overwriteOutput = 1
analysis_distance = 1000 # in meters
buff_geom = arcpy.Geometry() # to hold the 1km buffer, to get & set the extent

##temp_zonalstats_table = "IN_MEMORY/zStats"

offset_fieldname = "OFFSETA"

cellsize = '10'
env.cellsize = cellsize

##========== Stuff to keep track of ourselves =================================
maxRows = 25 # this is how often the script will restart. 250 seems overly cautious.
rowCount = 0
rowTimes = []
rowsProcessed = load_rows(processedrowset) # pickeld set object of rows we've already seen
current_row_oid = -1

##========== Prepare the cursor ===============================================
oid_field = arcpy.Describe(gnaf_fc).OIDFieldName

# This where clause is how we processed the different subsets on the different machines
# Each VM got different start & end rows to make sure the whole dataset was covered
start_row = 0
end_row = 100
where_clause = "%s BETWEEN %s AND %s" % (oid_field,start_row,end_row) # "" #
# These are the fields which will appear in the cursor
##cur_fields = ["OID@","GNAF_PID","lat","lon"] # <-- original fields
cur_fields = ["OID@","ID","LATITUDE","LONGITUDE"]

## Write the header row - will only happen if file doesn't already exist
write_header_row(out_csv)

##========== Go ===============================================================
with arcpy.da.SearchCursor(gnaf_fc,cur_fields,where_clause) as gnaf_cur:
    with open(out_csv,"a") as f: # note the file is opened in append mode
        writer = csv.writer(f)
        print "Commencing processing..."
        for row in gnaf_cur:
            if rowCount >= maxRows:
                print "...that's enough for now, exiting."
                break
            rowStartTime = time.time()

            # set the extent to be everything for the SelectByAttribute operation
            env.extent = "MAXOF" # DEBUG
            current_row_oid = row[0]

            print "\nProcessing Row %s..."%current_row_oid

            # Check if we've seen this before (this may or may not be doing anything)
            if current_row_oid in rowsProcessed:
                print "\t...definitely seen this one before, moving along."
                continue

            print "\t...haven't done this one yet..."

            print "\t...extracting vectors..."
            # Select the point we're looking at to its own temp dataset
            current_row_selector = '"%s" = %s'%(oid_field,current_row_oid)
            arcpy.Select_analysis(gnaf_fc,"tempTarget",current_row_selector)
            print "\t finished select"
            arcpy.FeatureClassToFeatureClass_conversion("tempTarget",r"C:\ExposureDetection\scratch", "Target.shp")# DEBUG
            print "Dumped to C:/temp"# DEBUG

            # Set the processing extent to an %analysis_distance buffer from current point
            buff_geom_list = arcpy.Buffer_analysis("tempTarget",
                                                   buff_geom,
                                                   "%s Meters"%analysis_distance)
            env.extent = buff_geom_list[0].extent

            # Extract the applicable hydro layer for this feature, only containing ocean
            arcpy.Clip_analysis(hydro_fc,buff_geom_list,"tempHydro")
            print "\t finished clip oceans" # DEBUG
            # Try using intersect instead of clip # DEBUG
##            arcpy.Intersect_analysis([hydro_fc,buff_geom_list],"tempHydro") # DEBUG
##            print "\t finished intersecting oceans" # DEBUG
            arcpy.FeatureClassToFeatureClass_conversion("tempHydro",r"C:\ExposureDetection\scratch", "Hydro.shp") # DEBUG
            print "\t finished dumping to temp" # DEBUG
            # Extract the nearby properties
            arcpy.Clip_analysis(gnaf_fc,buff_geom_list,"tempNearby")
            print "\t finished clip properties"
            # Remove intersecting features
            arcpy.SymDiff_analysis("tempTarget", "tempNearby","tempNearbyNoIntersect")
            print "\t finished SymDiff"
            # Set their elev to 10
            arcpy.CalculateField_management("tempNearbyNoIntersect",offset_fieldname,10,"PYTHON_9.3")
            print "\t finished Calc Field"

            # Do the processing
            stats = exposure("tempTarget","tempNearbyNoIntersect","tempHydro")

            # A little feedback to check it's working...
            if max(stats) > 0:
                print " Exposure detected! ".center(79,"*")
            else:
                print "\t...row not exposed..."

            # Write the results to our out CSV
            # the "(max(stats(,)" part is the "exposedatall" attribute
            # as a single-element tuple for easy cat'ing to the rest of the row
            # (will return 1 if exposed anywhere and 0 otherwise)
            writer.writerow(row[1:]+(max(stats),)+stats)

            # admin for tracking
            t = time.time()-rowStartTime
            rowTimes.append(t)
            rowCount += 1
            rowsProcessed.add(current_row_oid)
            save_object(rowsProcessed,processedrowset) #just overwrite

            print "\t...done. Row %s took "%current_row_oid, t, "seconds."
            print "%s rows processed this session so far."%rowCount

endtime = time.time()
print " Session complete ".center(79,"=")
print "Finished %s rows."%rowCount
try:
    print "Mean row time:",sum(rowTimes)/len(rowTimes),"seconds"
    print "Total time taken:",(endtime-starttime)/60,"minutes"

except ZeroDivisionError:
    print "\n"
    print "*"*79
    print "Zero rows processed - this machine has probably finished!".center(79,"*")
    print "*"*79
