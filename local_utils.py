#Define auxiliary functions for beacon list filtering. 

def f2m(item_in_feet):
    "converts feet to meters"
    try:
        return item_in_feet / 3.28084
    
    except TypeError:
        return float(item_in_feet) / 3.28084
    
    return converted


def find_near_transmitters(station_lat, station_lon, radius):
    import pandas as pd
    from geographiclib.geodesic import Geodesic

    beacon_types = ['VOR','VOR-DME','TACAN']
    path_to_navaids_file = "http://ourairports.com/data/navaids.csv"

    pd.set_option('display.max_rows', 1000)
    navaidsData = pd.read_csv(path_to_navaids_file, index_col='id')    
    
    transmitters = navaidsData[navaidsData['type'].isin(beacon_types)]
    near_transmitters = pd.DataFrame()
    
    for transmitter in transmitters.index:
        # calculate distance between transmitter and signal reception station. 
        distance = Geodesic.WGS84.Inverse(station_lat, station_lon, transmitters.ix[transmitter].latitude_deg, transmitters.ix[transmitter].longitude_deg)['s12']
        distance = distance / 1000 # convert distance to kilometers
        if distance < radius :
            near_transmitters = near_transmitters.append(transmitters.ix[transmitter])
    return near_transmitters

def dip_angle_lat(lat2,long2,alt2,lat1,long1,alt1,lcorr_for_refraction= True):

    from math import *

    ################################################################################
    # if we are going to correct for refraction, calculate the 
    # distance so we can correct alt2
    # yoeli_with_viewshed_refraction_paper
    # awesome_paper_on_viewshed_analysis
    # horizon_calculation.pdf
    ################################################################################
    
    theta1 = lat1*pi/180
    theta2 = lat2*pi/180
    phi1 = long1*pi/180
    phi2 = long2*pi/180
    
    a = 6378137.0 # earth ellisoid parameters
    b = 6356752.3
    
    R1 = sqrt(((a**2 * cos(theta1))**2 + (b**2*sin(theta1))**2)/((a*cos(theta1))**2+(b*sin(theta1))**2))
    R2 = sqrt(((a**2 * cos(theta2))**2 + (b**2*sin(theta2))**2)/((a*cos(theta2))**2+(b*sin(theta2))**2))
    
    ################################################################################
    # calculate the correction for refraction
    # yoeli_with_viewshed_refraction_paper
    # awesome_paper_on_viewshed_analysis
    # distance between latitude and longitude is
    # http://williams.best.vwh.net/avform.htm#Dist
    ################################################################################

    dlat = (theta2-theta1)
    dlon = (phi2-phi1)
    a = (sin(dlat/2)**2)+(sin(dlon/2)**2)*cos(theta1)*cos(theta2)
    c = 2*atan2(sqrt(a),sqrt(1-a))
    dist = R2*c
    if lcorr_for_refraction:
        k = 0.13
        adiff = (dist**2)*(k)/(R2*2)
    else:
        adiff = 0
 
    ################################################################################
    # now calculate the dip angle
    # this is a good online calculator view-source:http://cosinekitty.com/compass.html
    # first calculate the cartesian coordinates of each of the points
    # gamma is the angle between the vectors
    ################################################################################
    x1 = (R1+alt1)*cos(phi1)*sin(pi/2-theta1)
    y1 = (R1+alt1)*sin(phi1)*sin(pi/2-theta1)
    z1 = (R1+alt1)*cos(pi/2-theta1)
    x2 = (R2+alt2)*cos(phi2)*sin(pi/2-theta2)
    y2 = (R2+alt2)*sin(phi2)*sin(pi/2-theta2)
    z2 = (R2+alt2)*cos(pi/2-theta2)
    d1 = sqrt(x1**2+y1**2+z1**2)
    d2 = sqrt(x2**2+y2**2+z2**2)
    gamma = abs(acos((x1*x2+y1*y2+z1*z2)/(d1*d2)))
    
    ################################################################################
    # using law of cosines, calculate the distance between the points, w,
    # at the ends of the vectors, given that we know the angle between them
    # we also have
    # R2+h2+v = (R1+h1)/cos(gamma)
    # (thus v = (R1+h1)/cos(gamma)-R2-h2 )
    # and
    # w/sin(90-gamma) = v/sin(delta)
    # (thus sin(delta) = v*sin(90-gamma)/w)
    ################################################################################
    lnew = True
        
    if lnew:
        w = sqrt((R2+alt2+adiff)**2+(R1+alt1)**2-2*(R1+alt1)*(R2+alt2+adiff)*cos(gamma))
        v = ((alt1+R1)/cos(gamma)-R2-alt2-adiff)
        dip_angle = -asin(v*sin(pi/2-gamma)/w)
        dip_angle = dip_angle*180/pi
            
    else:
        l = sqrt((R2+alt2)**2+(R1+alt1)**2-2*(R1+alt1)*(R2+alt2+adiff)*cos(gamma))
        t = ((alt1+R1)/cos(gamma))
        dip_angle = asin((R2+alt2+adiff-t)/(l*sin(pi/2+gamma)))
        dip_angle = dip_angle*180/pi
    
    a = sin(phi2-phi1)*cos(theta2)
    b = cos(theta1)*sin(theta2) - sin(theta1)*cos(theta2)*cos(phi2-phi1)
    bearing = atan2(a,b)
    bearing = bearing*180/pi
    
    return(dip_angle,bearing,(w/1000.0))
