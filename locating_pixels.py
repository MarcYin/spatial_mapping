import numpy as np
import gdal
from osgeo import osr
import itertools


#a function transform between different projections....
def transform(a = 1):
    # from prof. lewis
    wgs84 = osr.SpatialReference( ) # Define a SpatialReference object
    wgs84.ImportFromEPSG( 4326 ) # And set it to WGS84 using the EPSG code
    modis_sinu = osr.SpatialReference() # define the SpatialReference object
    # In this case, we get the projection from a Proj4 string
    modis_sinu.ImportFromProj4 ( \
                    "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
    
    utm = osr.SpatialReference( )
    utm.ImportFromProj4 ('+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs')
    # add more ....
    
    if a ==1:
        # from modis to wgs 84
        tx = osr.CoordinateTransformation(modis_sinu, wgs84)
    elif a == 2:
        # from wgs 84 to modis
        tx = osr.CoordinateTransformation(wgs84,modis_sinu)
    elif a == 3:
        # from wgs 84 to modis
        tx = osr.CoordinateTransformation(wgs84, utm)
    elif a == 4:
        tx = osr.CoordinateTransformation(utm, wgs84)
    # the following two transforms shold be tested before the usage....
    # generally WGS84 is reccommonded as the base to convert between
    # different projections
    elif a == 5:
        tx = osr.CoordinateTransformation(modis_sinu, utm)
    elif a == 6:
        tx = osr.CoordinateTransformation(utm, modis_sinu)
    # elif: 
    #   even more .....
    else:
        tx = 0
        
        print 'please define your own transformation....'
    
    return tx
# since we have the origins and the steps, we can write a function
# to calculate modis tile number with the lat and lon as the inputs




x_step = -463.31271653
y_step = 463.31271653
m_y0, m_x0 = -20015109.354, 10007554.677
def mtile_cal(lat, lon):
    # a function calculate the tile number for MODIS, based on the lat and lon
    tx = transform( a = 2)# from wgs84 to modis 
    ho,vo,z = tx.TransformPoint(lon, lat)# still use the function instead of using the equation....
    h = int((ho-m_y0)/(2400*y_step))
    v = int((vo-m_x0)/(2400*x_step))
    return h,v

def get_m_corners(h,v):
    tx = transform(a=1) # from modis to wgs84
    # Work out the boundaries of the new dataset in the target projection
    x_size, y_size = 2400, 2400
    geo_t0, geo_t3 = m_y0 + h*2400*y_step , m_x0 + v*2400*x_step

    (ul_lon, ul_lat, ulz ) = tx.TransformPoint( geo_t0, geo_t3)

    (lr_lon, lr_lat, lrz ) = tx.TransformPoint( geo_t0 + y_step*y_size, \
                                          geo_t3 + x_step*x_size )

    (ll_lon, ll_lat, llz )  = tx.TransformPoint( geo_t0 , \
                                          geo_t3 + x_step*x_size )

    (ur_lon, ur_lat, urz ) = tx.TransformPoint( geo_t0 + y_step*y_size, \
                                          geo_t3  )
    
    return ul_lon, ul_lat, lr_lon, lr_lat, ll_lon, ll_lat, ur_lon, ur_lat


def bilineanr(coords, dic, cors):
    
    '''
    basically a bilinear interpolation in matrix form (https://en.wikipedia.org/wiki/Bilinear_interpolation) 
    --------------------------------------------------------------------------------------------------------------
    coords is the coordinates [(x1, y1), (x2, y2)...] needed to transfer
    
    dic is a dictionary of the Upper lfet (UL), UR, LL, LR 's lat and lons
    
    example: dic ={'LL_LAT': 36.35288,
                   'LL_LON': 113.00651,
                   'LR_LAT': 36.41186,
                   'LR_LON': 115.6326,
                   'UL_LAT': 38.51077,
                   'UL_LON': 112.88999,
                   'UR_LAT': 38.57451,
                   'UR_LON': 115.59258}
    cors = {'ulx':ulx, 'llx': llx, 'urx': urx, 'lrx': lrx, 'uly': uly, 'lly': lly, 'ury': ury, 'lry': lry}
    
    corners is the (x, y) corresponding to the shape of the area (array)            
    '''
    
    a = np.matrix([[1, dic['UL_LAT'], dic['UL_LON'], dic['UL_LAT']*dic['UL_LON']],
                   [1, dic['LL_LAT'], dic['LL_LON'], dic['LL_LAT']*dic['LL_LON']],
                   [1, dic['UR_LAT'], dic['UR_LON'], dic['UR_LAT']*dic['UR_LON']],
                   [1, dic['LR_LAT'], dic['LR_LON'], dic['LR_LAT']*dic['LR_LON']],
                   ])
    
    convs = np.ones((4,len(coords[0])))
    convs[1] = coords[0]
    convs[2] = coords[1] 
    convs[3] = (coords[0]* coords[1])
    convs = np.matrix(convs)

    x = np.matrix([cors['ulx'], cors['llx'], cors['urx'], cors['lrx']])*((a**-1).T)*convs
    y = np.matrix([cors['uly'], cors['lly'], cors['ury'], cors['lry']])*((a**-1).T)*convs
    
    return np.array(x).ravel(), np.array(y).ravel()

def Find_corresponding_pixels(H_res_fname, destination_res=500):
    
    '''
    A function for the finding of corresponding pixels indexes 
    in H (sentinel 2) and L (MODIS) resolution image.
    
    args:
        H_res_fname -- the high resolution image filename, need to have geoinformation enbeded in the file
        destination_res -- for the calculation of pixel number in one MODIS tile
    return:
        index: a dictionary contain both the MODIS tile name and pixels indexes
    
    '''
    
    if destination_res % 250 == 0:
        pass
    else:
        print 'destination resolution can only be 250, 500 and 1000 !!!'
        raise IOError

    g = gdal.Open(H_res_fname)
    geo_t = g.GetGeoTransform()
    x_size, y_size = g.RasterXSize, g.RasterYSize

    wgs84 = osr.SpatialReference( ) # Define a SpatialReference object
    wgs84.ImportFromEPSG( 4326 ) # And set it to WGS84 using the EPSG code
    H_res_geo = osr.SpatialReference( )
    raster_wkt = g.GetProjection()
    H_res_geo.ImportFromWkt(raster_wkt)
    tx = osr.CoordinateTransformation(H_res_geo, wgs84)
    # so we need the four corners coordiates to check whether they are within the same modis tile
    (ul_lon, ul_lat, ulz ) = tx.TransformPoint( geo_t[0], geo_t[3])

    (lr_lon, lr_lat, lrz ) = tx.TransformPoint( geo_t[0] + geo_t[1]*x_size, \
                                          geo_t[3] + geo_t[5]*y_size )

    (ll_lon, ll_lat, llz ) = tx.TransformPoint( geo_t[0] , \
                                          geo_t[3] + geo_t[5]*y_size )

    (ur_lon, ur_lat, urz ) = tx.TransformPoint( geo_t[0] + geo_t[1]*x_size, \
                                          geo_t[3]  )

    #print (ul_lon, ul_lat), (lr_lon, lr_lat), (ll_lon, ll_lat), (ur_lon, ur_lat)
    # now its the s2 corners latitudes and longtitudes
    s_dic ={'UR_LAT': ur_lat,
           'UR_LON': ur_lon,
           'LR_LAT': lr_lat,
           'LR_LON': lr_lon,
           'UL_LAT': ul_lat,
           'UL_LON': ul_lon,
           'LL_LAT': ll_lat,
           'LL_LON': ll_lon}
    s_cors = ulx, uly, lrx, lry, llx, lly, urx, ury = 0,0,x_size, y_size, x_size,0, 0, y_size 
    s_corners = dict(zip(['ulx', 'uly', 'lrx', 'lry', 'llx', 'lly', 'urx', 'ury'], s_cors))
    
    
    
    a0, b0 = None, None
    corners = [(ul_lon, ul_lat), (lr_lon, lr_lat), (ll_lon, ll_lat), (ur_lon, ur_lat)]
    tiles = []
    for i,j  in enumerate(corners):
        h,v = mtile_cal(j[1], j[0])
        if (h==a0) &(v==b0):
            pass
        else:
            tiles.append([i,h,v]) # 0--ul;1--lr;2--ll;3--ur
            a0, b0 = h,v
    
    # The modis defults
    pix_num = 4800/(destination_res/250)
    xs = np.arange(0,pix_num).repeat(pix_num)
    ys = np.tile(np.arange(0,pix_num), pix_num)
    m_coords = [xs,ys]

    # corners coordinate in image
    m_dic ={'UR_LAT': 0,
           'UR_LON': pix_num-1,
           'LR_LAT': pix_num-1,
           'LR_LON': pix_num-1,
           'UL_LAT': 0,
           'UL_LON': 0,
           'LL_LAT': pix_num-1,
           'LL_LON': 0}

    inds = {}
    for i in tiles:
        m_cors = get_m_corners(i[1],i[2])
        m_corners = dict(zip(['uly', 'ulx', 'lry', 'lrx', 'lly', 'llx', 'ury', 'urx'], m_cors))
        latitudes, longtitudes = bilineanr(m_coords, m_dic, m_corners)

        cors = bilineanr([latitudes, longtitudes], s_dic, s_corners)

        hinds =np.array([cors[0][(cors[0]>=0)&(cors[0]<x_size)&(cors[1]>=0)&(cors[1]<y_size)],
                     cors[1][(cors[0]>=0)&(cors[0]<x_size)&(cors[1]>=0)&(cors[1]<y_size)]]).astype(int)
        minds = np.where(((cors[0]>=0)&(cors[0]<x_size)&(cors[1]>=0)&(cors[1]<y_size)).reshape((pix_num,pix_num)))
        inds['h%02dv%02d'%(i[1],i[2])] = [hinds, np.array(minds)]

    return inds
