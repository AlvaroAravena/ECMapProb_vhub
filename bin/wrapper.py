import Rappture
import elevation
import tifffile
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, atan2, radians, log, factorial, tan
import sys
import os
from PIL import Image, ImageDraw

# Auxiliary functions

def distance_two_points(lat1, lat2, lon1, lon2):

	R = 6373.0

	lat1 = radians(lat1)
	lon1 = radians(lon1)
	lat2 = radians(lat2)
	lon2 = radians(lon2)

	dlon = lon2 - lon1
	dlat = lat2 - lat1

	a = sin(dlat / 2.0) ** 2.0 + cos(lat1) * cos(lat2) * sin(dlon / 2.0) ** 2.0
	c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
	return ( R * c ) * 1000.0

def interpol_pos(lon1, lat1, step_lon_deg, step_lat_deg, lon_cen, lat_cen, cells_lon, cells_lat, Topography):

	dlon = int(np.floor( (lon_cen - lon1 )/ (step_lon_deg) ))
	dlat = (cells_lat - 2) - int(np.floor( (lat_cen - lat1) / (step_lat_deg) ))

	if(dlon >= ( cells_lon - 1.0 ) or dlat >= ( cells_lat - 1.0 ) or dlon < 0.0 or dlat < 0.0):
		return 99999

	aux_lon = 2.0 * ( lon_cen - ( dlon * step_lon_deg + lon1 ) - step_lon_deg / 2.0 ) / step_lon_deg
	aux_lat = 2.0 *( - lat_cen + ( (cells_lat - 1.0 - dlat) * step_lat_deg + lat1 ) - step_lat_deg / 2.0 ) / step_lat_deg

	dc = ( Topography[dlat][dlon] + Topography[dlat][dlon+1] + Topography[dlat+1][dlon] + Topography[dlat+1][dlon+1] ) / 4
	[x3, y3, z3] = [0.0, 0.0, dc]

	if( aux_lon >= 0.0 and abs(aux_lon) >= abs(aux_lat) ):
		[x1,y1,z1] = [1.0, 1.0, Topography[dlat+1][dlon+1]] 
		[x2,y2,z2] = [1.0, -1.0, Topography[dlat][dlon+1]] 
	elif( aux_lat >= 0.0 and abs(aux_lon) < abs(aux_lat) ):
		[x1,y1,z1] = [-1.0, 1.0, Topography[dlat+1][dlon]] 
		[x2,y2,z2] = [1.0, 1.0, Topography[dlat+1][dlon+1]] 
	elif( aux_lon < 0.0 and abs(aux_lon) >= abs(aux_lat) ):
		[x1,y1,z1] = [-1.0, 1.0, Topography[dlat+1][dlon]] 
		[x2,y2,z2] = [-1.0, -1.0, Topography[dlat][dlon]] 
	else:
		[x1,y1,z1] = [-1.0, -1.0, Topography[dlat][dlon]] 
		[x2,y2,z2] = [1.0, -1.0, Topography[dlat][dlon+1]]
 
	f1 = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1)
	f2 = (z2-z1)*(x3-x1) - (z3-z1)*(x2-x1)
	f3 = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

	return ((- aux_lon * f1 - aux_lat * f2) / f3 + dc)

##########################
###### MAIN PROGRAM ######
##########################

# INPUT PARAMETERS

# READING INPUT PARAMETERS
print "Reading input parameters \n"

inDeckName = sys.argv[1]
lib = Rappture.library(inDeckName)

run_name = lib.get('input.string(run_name).current')
topo_source = lib.get('input.choice(topography).current')

if(topo_source == 'Upload_UTM'):
	DEM = lib.get('input.string(uploadedFile_UTM).current')
	east_cen = float(lib.get('input.number(east).current'))
	north_cen = float(lib.get('input.number(north).current'))
elif(topo_source == 'Upload_deg'):
	DEM = lib.get('input.string(uploadedFile_deg).current')
	lon_cen = float(lib.get('input.number(lonc).current'))
	lat_cen = float(lib.get('input.number(latc).current'))
else:
	lon1 = float(lib.get('input.number(lon1).current'))
	lon2 = float(lib.get('input.number(lon2).current'))
	lat1 = float(lib.get('input.number(lat1).current'))
	lat2 = float(lib.get('input.number(lat2).current'))
	lon_cen = float(lib.get('input.number(lonc).current'))
	lat_cen = float(lib.get('input.number(latc).current'))

source_type = lib.get('input.choice(source_type).current')

if(source_type == 'Linear'):
	azimuth_lin = float(lib.get('input.number(azimuth_lin).current'))
	length_lin = float(lib.get('input.number(length_lin).current'))

if(source_type == 'Radial'):
	radius_rad = float(lib.get('input.number(radius_rad).current'))
	ang1_rad = float(lib.get('input.number(ang1_rad).current'))
	ang2_rad = float(lib.get('input.number(ang2_rad).current'))

var_cen = float(lib.get('input.number(varc).current'))
height = float(lib.get('input.number(h).current'))
var_height = float(lib.get('input.number(varh).current'))
hl = float(lib.get('input.number(hl).current'))
var_hl = float(lib.get('input.number(varhl).current'))
N = int(lib.get('input.integer(N).current'))
cone_levels = int(lib.get('input.integer(c_levels).current'))

input_dist = lib.get('input.choice(input_dist).current')
resolution = lib.get('input.choice(resolution).current')
if(resolution == 'Ten'):
	angstep = 10
elif(resolution == 'Twenty'):
	angstep = 20

log = ''
current_path = os.getcwd()

# IMPORT MAP
if(topo_source == 'SRTM'):

	aux_lon = np.array([lon1, lon2])
	aux_lat = np.array([lat1, lat2])
	lon1 = min(aux_lon)
	lon2 = max(aux_lon)
	lat1 = min(aux_lat)
	lat2 = max(aux_lat)
	elevation.clip(bounds=(lon1, lat1, lon2, lat2), output = current_path + '/' + run_name + '.tif')

# READ MAP
if(topo_source == 'SRTM'):
	fp = run_name + '.tif'
	with tifffile.TIFFfile(fp) as tif:
		data = tif.asarray()
		for page in tif:
			for tag in page.tags.values():
				t = tag.name, tag.value
			image = page.asarray()
	elevation.clean()
	Topography = np.array(image)
	Topography_Save = Topography

	Topography_Sea = Topography + 0.0
	Topography_Sea[ Topography_Sea[:,:] <= 0] = -1.0 * np.sqrt(-1.0 * Topography_Sea[ Topography_Sea[:,:] <= 0])
	Topography_Sea[ Topography_Sea[:,:] > 0] =  np.nan
	Topography_Sea = Topography_Sea * -1.0

	Topography  = (Topography  + abs(Topography)) / 2.0
	cells_lon = Topography.shape[1]
	cells_lat = Topography.shape[0]

if(topo_source == 'Upload_UTM'):
	n_north = -1
	n_east = -1
	cellsize = -1
	indexini = -1
	nodata = -1
	lines = DEM.splitlines()

	for i in range(0,10):
		aux = lines[i].split()

		if(aux[0] == 'ncols'):
			n_north = int(aux[1])
		if(aux[0] == 'nrows'):
			n_east = int(aux[1])
		if(aux[0] == 'cellsize'):
			cellsize = float(aux[1])
		if(aux[0] == 'xllcorner'):
			east_cor = float(aux[1])
		if(aux[0] == 'yllcorner'):
			north_cor = float(aux[1])
		if(aux[0] == 'NODATA_value'):
			nodata = float(aux[1])
		if(len(aux) >= 10):
			indexini = i
			break

	Topography = np.zeros((n_north,n_east))
	for i in range(indexini, indexini + n_north):
		aux = lines[i].split()
		for j in range(0, n_east):
			Topography[i-indexini,j] = float(aux[j])
	Topography_Save = Topography

	Topography_Sea = Topography + 0.0
	Topography_Sea[ Topography_Sea[:,:] <= 0] = -1.0 * np.sqrt(-1.0 * Topography_Sea[ Topography_Sea[:,:] <= 0])
	Topography_Sea[ Topography_Sea[:,:] > 0] =  np.nan
	Topography_Sea = Topography_Sea * -1.0
	Topography  = (Topography  + abs(Topography)) / 2.0

if(topo_source == 'Upload_deg'):

	lon1 = -1
	lon2 = -1
	lat1 = -1
	lat2 = -1
	cells_lon = -1
	cells_lat = -1
	lines = DEM.splitlines()

	for i in range(0,10):
		aux = lines[i].split()
		if(aux[0] == 'lon1'):
			lon1 = float(aux[1])
		if(aux[0] == 'lon2'):
			lon2 = float(aux[1])
		if(aux[0] == 'lat1'):
			lat1 = float(aux[1])
		if(aux[0] == 'lat2'):
			lat2 = float(aux[1])
		if(aux[0] == 'cells_lon'):
			cells_lon = int(aux[1])
		if(aux[0] == 'cells_lat'):
			cells_lat = int(aux[1])
		if(len(aux) >= 10):
			indexini = i
			break

	Topography = np.zeros((cells_lat,cells_lon))
	for i in range(indexini, indexini + cells_lat):
		aux = lines[i].split()
		for j in range(0, cells_lon):
			Topography[i-indexini,j] = float(aux[j])
	Topography_Save = Topography

	Topography_Sea = Topography + 0.0
	Topography_Sea[ Topography_Sea[:,:] <= 0] = -1.0 * np.sqrt(-1.0 * Topography_Sea[ Topography_Sea[:,:] <= 0])
	Topography_Sea[ Topography_Sea[:,:] > 0] =  np.nan
	Topography_Sea = Topography_Sea * -1.0
	Topography  = (Topography  + abs(Topography)) / 2.0

# DEFINE THE MATRIX OF COORDINATES

if(topo_source == 'SRTM' or topo_source == 'Upload_deg' ):
	distance_lon = distance_two_points(lat1,lat1,lon1,lon2)
	distance_lat = distance_two_points(lat1,lat2,lon1,lon1)

	step_lon_m = distance_lon / (cells_lon-1)
	step_lat_m = distance_lat / (cells_lat-1)

	matrix_lon = np.zeros((cells_lat,cells_lon))
	matrix_lat = np.zeros((cells_lat,cells_lon))

	for i in range(0,cells_lon):
		matrix_lon[:,i] = lon1 + (lon2 - lon1) * (i) / (cells_lon-1)
	for j in range(0,cells_lat):
		matrix_lat[j,:] = lat1 + (lat2 - lat1) * (cells_lat-1-j) / (cells_lat-1)

	step_lon_deg = (lon2 - lon1)/(cells_lon - 1)
	step_lat_deg = (lat2 - lat1)/(cells_lat - 1)

if(topo_source == 'Upload_UTM' ):
	matrix_north = np.zeros((n_north,n_east))
	matrix_east = np.zeros((n_north,n_east))

	for i in range(0,n_east):
		matrix_east[:,i] = (east_cor + cellsize * i)
	for j in range(0,n_north):
		matrix_north[j,:] = (north_cor + cellsize * j)
	matrix_north = matrix_north[ range(len(matrix_north[:,0]) -1 , -1 , -1 ) , : ]

# CREATE VECTORS OF INPUT PARAMETERS AND DELETE NEGATIVE DATA
print 'Creating input vectors'

if(var_height > 0.0):
	if(input_dist == 'Gaussian'):
		height_vector = np.random.normal(height,var_height,N)
	if(input_dist == 'Uniform'):
		height_vector = np.random.uniform(height - var_height, height + var_height, N)
else:
	height_vector = np.ones(N) * height

if(var_hl > 0.0):
	if(input_dist == 'Gaussian'):
		hl_vector = np.random.normal(hl,var_hl,N)
	if(input_dist == 'Uniform'):
		hl_vector = np.random.uniform(hl - var_hl, hl + var_hl, N)
else:
	hl_vector = np.ones(N) * hl

if(var_height > 0.0):
	while( 1 == 1 ):
		aux_boolean = 0
		for i in range(0,N):
			if(height_vector[i] < 0):
				if(input_dist == 'Gaussian'):
					height_vector[i] = np.random.normal(height,var_height,1)
				if(input_dist == 'Uniform'):
					height_vector[i] = np.random.uniform(height - var_height, height + var_height,1)			
				aux_boolean = 1
		if(aux_boolean == 0):
			break

if(var_hl > 0.0):
	while( 1 == 1 ):
		aux_boolean = 0
		for i in range(0,N):
			if(hl_vector[i] < 0.05):
				if(input_dist == 'Gaussian'):
					hl_vector[i] = np.random.normal(hl,var_hl,1)
				if(input_dist == 'Uniform'):
					hl_vector[i] = np.random.uniform(hl - var_hl, hl + var_hl,1)
				aux_boolean = 1
		if(aux_boolean == 0):
			break

if(topo_source == 'SRTM'  or topo_source == 'Upload_deg' ):
	if( var_cen > 0.0):
		if(input_dist == 'Gaussian'):
			lon_cen_vector = np.random.normal(lon_cen, var_cen * step_lon_deg / step_lon_m , N)
			lat_cen_vector = np.random.normal(lat_cen, var_cen * step_lat_deg / step_lat_m , N)
		if(input_dist == 'Uniform'):
			lon_cen_vector = np.random.uniform(lon_cen - var_cen * step_lon_deg / step_lon_m, lon_cen + var_cen * step_lon_deg / step_lon_m, N)
			lat_cen_vector = np.random.uniform(lat_cen - var_cen * step_lat_deg / step_lat_m, lat_cen + var_cen * step_lat_deg / step_lat_m, N)
			while( 1 == 1 ):
				aux_boolean = 0
				for i in range(0,N):
					if(np.power((lon_cen_vector[i] - lon_cen) * step_lon_m / step_lon_deg ,2) + np.power((lat_cen_vector[i] - lat_cen) * step_lat_m / step_lat_deg , 2) > np.power(var_cen,2)):
						lon_cen_vector[i]  = np.random.uniform(lon_cen - var_cen * step_lon_deg / step_lon_m, lon_cen + var_cen * step_lon_deg / step_lon_m, 1)
						lat_cen_vector[i]  = np.random.uniform(lat_cen - var_cen * step_lat_deg / step_lat_m, lat_cen + var_cen * step_lat_deg / step_lat_m, 1)
						aux_boolean = 1
				if(aux_boolean == 0):
					break

	else:
		lon_cen_vector = np.ones(N) * lon_cen
		lat_cen_vector = np.ones(N) * lat_cen

	if(source_type == 'Linear'):
		pos_structure = np.random.uniform(-1,1,N)
		lon_cen_vector = lon_cen_vector + pos_structure * np.sin(azimuth_lin * np.pi/180) * length_lin * step_lon_deg / step_lon_m
		lat_cen_vector = lat_cen_vector + pos_structure * np.cos(azimuth_lin * np.pi/180) * length_lin * step_lat_deg / step_lat_m

	if(source_type == 'Radial'):
		pos_structure = ang1_rad + np.random.uniform(0,1,N)*(ang2_rad - ang1_rad)
		lon_cen_vector = lon_cen_vector + np.cos(pos_structure * np.pi/180) * radius_rad * step_lon_deg / step_lon_m
		lat_cen_vector = lat_cen_vector + np.sin(pos_structure * np.pi/180) * radius_rad * step_lat_deg / step_lat_m

if(topo_source == 'Upload_UTM'):

	if( var_cen > 0.0):
		if(input_dist == 'Gaussian'):
			east_cen_vector = np.random.normal(east_cen,var_cen,N)
			north_cen_vector = np.random.normal(north_cen,var_cen,N)
		if(input_dist == 'Uniform'):
			east_cen_vector = np.random.uniform(east_cen - var_cen, east_cen + var_cen, N)
			north_cen_vector = np.random.uniform(north_cen - var_cen, north_cen + var_cen,N)
			while( 1 == 1 ):
				aux_boolean = 0
				for i in range(0,N):
					if(np.power((east_cen_vector[i] - east_cen) ,2) + np.power((north_cen_vector[i] - north_cen) , 2) > np.power(var_cen,2)):
						east_cen_vector[i]  = np.random.uniform(east_cen - var_cen , east_cen + var_cen , 1)
						north_cen_vector[i]  = np.random.uniform(north_cen - var_cen, north_cen + var_cen , 1)
						aux_boolean = 1
				if(aux_boolean == 0):
					break

	else:
		east_cen_vector = np.ones(N) * east_cen
		north_cen_vector = np.ones(N) * north_cen

	if(source_type == 'Linear'):
		pos_structure = np.random.uniform(-1,1,N)
		east_cen_vector = east_cen_vector + pos_structure * np.sin(azimuth_lin * np.pi/180) * length_lin
		north_cen_vector = north_cen_vector + pos_structure * np.cos(azimuth_lin * np.pi/180) * length_lin

	if(source_type == 'Radial'):
		pos_structure = ang1_rad + np.random.uniform( 0 , 1 , N ) * ( ang2_rad - ang1_rad )
		east_cen_vector = east_cen_vector  + np.cos(pos_structure * np.pi/180 ) * radius_rad 
		north_cen_vector = north_cen_vector + np.sin(pos_structure * np.pi/180 ) * radius_rad

# ENERGY CONES
anglen = 360 / angstep
pix_min = 0.0
distep = 10

factor_mult = 50.0
center_elim = 0.5
aux_backward = 1 / (1 + np.exp(factor_mult * (np.linspace(0.0, 1.0, anglen/2 + 1) - center_elim) ) )
vector_backward_1 = np.zeros(anglen)
vector_backward_1[0:anglen/2 - 1] = aux_backward[anglen/2-1:0:-1]
vector_backward_1[anglen/2-1:] = aux_backward[:]
vector_backward_1[vector_backward_1 < 1e-3] = 0
vector_backward_1[vector_backward_1 > 1.0 - 1e-3] = 1.0
aux_backward = 1 / (1 + np.exp(factor_mult * (np.linspace(1.0/(anglen/2), 1.0 - 1.0/(anglen/2), anglen/2 ) - center_elim) ) )
vector_backward_2 = np.zeros(anglen)
vector_backward_2[0:anglen/2] = aux_backward[::-1]
vector_backward_2[anglen/2:] = aux_backward[:]
vector_backward_2[vector_backward_2 < 1e-3] = 0
vector_backward_2[vector_backward_2 > 1.0 - 1e-3] = 1.0
index_max = anglen/2 - 1
vector_correc = np.zeros(anglen)

if(topo_source == 'SRTM'  or topo_source == 'Upload_deg'):
	data_cones = np.zeros((cells_lat,cells_lon))
	data_cones_save = np.zeros((cells_lat,cells_lon))
	data_aux_t = np.ones((cells_lat,cells_lon))
	data_aux_b = np.zeros((cells_lat,cells_lon))
	vec_ang = range(0, 360, angstep)

	for i in range(0,N):

		current_level = 0
		data_step = np.zeros((cells_lat,cells_lon))
		polygon = []
		height_eff = height_vector[i] + interpol_pos(lon1, lat1, step_lon_deg, step_lat_deg, lon_cen_vector[i], lat_cen_vector[i], cells_lon, cells_lat, Topography)
		polygon.append((lon_cen_vector[i], lat_cen_vector[i],  height_eff, 1.0, -1, height_vector[i] ))
		sum_pixels = 0
		hl_current = hl_vector[i]		

		for j in range(10000): 

			if(j == len(polygon)):
				if( N == 1 ):			
					data_cones = data_cones + data_step
				break
			if( cone_levels < polygon[j][3] ):
				if( N == 1 ):
					data_cones = data_cones + data_step
				break
			elif(current_level < polygon[j][3]):
				current_level = polygon[j][3]
				if( N == 1 ):
					data_cones = data_cones + data_step

			polygon_xy = []
			polygons_new = []

			for angle_deg in vec_ang:
				angle_rad = angle_deg * np.pi /180
				for distance in range(0, 100000, distep):
					h = interpol_pos(lon1, lat1, step_lon_deg, step_lat_deg, polygon[j][0] + distance * cos(angle_rad) * step_lon_deg / step_lon_m , polygon[j][1] + distance*sin(angle_rad)*step_lat_deg/step_lat_m , cells_lon, cells_lat, Topography)
					if( h >= polygon[j][2] - hl_current * distance ):
						polygon_xy.append((int((polygon[j][0] + (distance - distep)*cos(angle_rad)*step_lon_deg/step_lon_m - lon1) * cells_lon / (lon2 - lon1)),int((polygon[j][1] + (distance - distep)*sin(angle_rad)*step_lat_deg/step_lat_m - lat1) * cells_lat / (lat2 - lat1))))
						polygons_new.append(distance - distep)
						break

			if( polygon[j][4] > -1 ):
				lim = np.int(polygon[j][4])
				if( polygon[j][4] == np.int(polygon[j][4]) ):
					for ii in range(anglen):
						vector_correc[ii] = vector_backward_1[int((ii - polygon[j][4] + index_max) % anglen)]

				else:
					for ii in range(anglen):
						vector_correc[ii] = vector_backward_2[int((ii - polygon[j][4] + index_max) % anglen)]

				polygons_new = polygons_new * vector_correc

			img = Image.new('L', (cells_lon, cells_lat), 0)
			if( len(polygon_xy) > 0 ):
				draw = ImageDraw.Draw(img).polygon(polygon_xy, outline = 1 , fill = 1)
				data_step = np.maximum( np.minimum(data_aux_t, data_step + np.array(img)), data_aux_b)

			if( cone_levels > polygon[j][3] and sum(sum(data_step)) > sum_pixels + pix_min ):

				aux = np.zeros(len(polygons_new)+2) 
				aux[1:len(polygons_new)+1] = np.array(polygons_new) 
				aux[0] = polygons_new[len(polygons_new)-1]
				aux[len(polygons_new)+1] = polygons_new[0]
				der1 = (aux[1:len(aux)-1] - aux[2:len(aux)])
				der2 = (aux[1:len(aux)-1] - aux[0:len(aux)-2])
				wh1 = np.where(der1 >= 0)
				wh2 = np.where(der2 >= 0)
				wh_max = np.intersect1d(wh1[0], wh2[0])
				wh_grouped = np.split(wh_max, np.where(np.diff(wh_max) > 1)[0] + 1 )
				wh3 = np.where( abs(der1) > 0)
				wh4 = np.where( abs(der2) > 0)
				wh5 = np.intersect1d(wh_max, wh3[0])
				wh6 = np.intersect1d(wh_max, wh4[0])
				grouped_filter = np.zeros(len(wh_grouped))

				for x_grouped in range(len(wh_grouped)):
					if( len(np.intersect1d(wh_grouped[x_grouped],wh5)) > 0 and len(np.intersect1d(wh_grouped[x_grouped],wh6)) > 0):
						grouped_filter[x_grouped] = 1

				if( np.min(wh_grouped[0]) == 0 and np.max(wh_grouped[len(wh_grouped)-1]) == anglen - 1):

					if( len(np.intersect1d(wh_grouped[0],wh5)) > 0 and len(np.intersect1d(wh_grouped[len(wh_grouped)-1],wh6)) > 0):
						grouped_filter[len(wh_grouped) - 1] = 1

					aux_grouped = np.concatenate((wh_grouped[len(wh_grouped)-1], wh_grouped[0] + len(polygons_new)))
					aux_filter = grouped_filter[len(wh_grouped)-1] + grouped_filter[0]
					wh_grouped = wh_grouped[1:-1]
					wh_grouped.append(aux_grouped)
					grouped_filter = np.append(grouped_filter[1:-1],aux_filter)

				wh_max = []
				for k in range(len(grouped_filter)):
					if(grouped_filter[k] > 0 ):
						if(np.mean(wh_grouped[k]) < len(polygons_new) and np.mean(wh_grouped[k]) >= 0.0):
							wh_max.append(np.mean(wh_grouped[k]))
						elif( np.mean(wh_grouped[k]) < len(polygons_new) ):
							wh_max.append(len(polygons_new) + np.mean(wh_grouped[k]))
						else:
							wh_max.append(- len(polygons_new) + np.mean(wh_grouped[k]))

				wh1 = np.where(der1 <= 0)
				wh2 = np.where(der2 <= 0)
				wh_min = np.intersect1d(wh1[0], wh2[0])
				wh_grouped = np.split(wh_min, np.where(np.diff(wh_min) > 1)[0] + 1 )
				wh3 = np.where( abs(der1) > 0)
				wh4 = np.where( abs(der2) > 0)
				wh5 = np.intersect1d(wh_min, wh3[0])
				wh6 = np.intersect1d(wh_min, wh4[0])
				grouped_filter = np.zeros(len(wh_grouped))

				for x_grouped in range(len(wh_grouped)):
					if( len(np.intersect1d(wh_grouped[x_grouped],wh5)) > 0 and len(np.intersect1d(wh_grouped[x_grouped],wh6)) > 0):
						grouped_filter[x_grouped] = 1
					
				if( np.min(wh_grouped[0]) == 0 and np.max(wh_grouped[len(wh_grouped)-1]) == anglen - 1):
					if( len(np.intersect1d(wh_grouped[0],wh5)) > 0 and len(np.intersect1d(wh_grouped[len(wh_grouped)-1],wh6)) > 0):
						grouped_filter[len(wh_grouped) - 1] = 1

					aux_grouped = np.concatenate((wh_grouped[len(wh_grouped)-1], wh_grouped[0] + len(polygons_new)))
					aux_filter = grouped_filter[len(wh_grouped)-1] + grouped_filter[0]
					wh_grouped = wh_grouped[1:-1]
					wh_grouped.append(aux_grouped)
					grouped_filter = np.append(grouped_filter[1:-1],aux_filter)

				wh_min = []

				for k in range(len(grouped_filter)):
					if(grouped_filter[k] > 0 ):
						if(np.mean(wh_grouped[k]) < len(polygons_new) and np.mean(wh_grouped[k]) >= 0.0):
							wh_min.append(np.mean(wh_grouped[k]))
						elif(np.mean(wh_grouped[k]) < len(polygons_new) ):
							wh_min.append(len(polygons_new) + np.mean(wh_grouped[k]))
						else:
							wh_min.append(- len(polygons_new) + np.mean(wh_grouped[k]) )

				wh_sum = np.zeros(len(polygons_new)) 

				if(len(wh_max) > 0):
					
					if( len(wh_max) == 1):
						for l_max_real in wh_max:
							lmax = np.int(l_max_real)
							l_it = 	len(polygons_new) - 1		
							for l in range(1,len(polygons_new)):
								l_index = lmax + l
								if(l_index >= len(polygons_new)):
									l_index = l_index - len(polygons_new)
								if( polygons_new[lmax] < polygons_new[l_index] ):
									l_it = l
									break
								wh_sum[lmax] = wh_sum[lmax] + (polygons_new[lmax] - polygons_new[l_index])							
							for l in range(1,len(polygons_new) - l_it):
								l_index = lmax - l
								if(l_index < 0):
									l_index = l_index + len(polygons_new)
								if( polygons_new[lmax] < polygons_new[l_index] ):
									break							
								wh_sum[lmax] = wh_sum[lmax] + (polygons_new[lmax] - polygons_new[l_index])

					else:

						wh_max = np.sort(wh_max)
						wh_min = np.sort(wh_min)

						if(wh_min[0] > wh_max[0]):

							for l_ind in range(len(wh_max)):
								l_max_real = wh_max[l_ind]	
								l_max_int = np.int(l_max_real)
								step_right = wh_min[l_ind] - l_max_int
								l_right_real = wh_min[l_ind]
								l_right_int = np.int(l_right_real)

								if(l_ind == 0):
									step_left = anglen + l_max_int - wh_min[len(wh_min)-1]
									l_left_real = wh_min[len(wh_min) - 1]
									left_index = len(wh_min) - 1
								else:
									step_left = l_max_int - wh_min[l_ind - 1]
									l_left_real = wh_min[l_ind - 1]
									left_index = l_ind - 1
								
								l_left_int = np.int(l_left_real)

								for l in range(1,int(step_right)):
									l_index = l_max_int + l
									if(l_index >= len(polygons_new)):
										l_index = l_index - len(polygons_new)
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])

								if( int(step_right) == step_right ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_right_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_right_int])

								for l in range(1,int(step_left)):
									l_index = l_max_int - l
									if( l_index < 0 ):
										l_index = len(polygons_new) + l_index
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])

								if( int(step_left) == step_left ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_left_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_left_int])

						else:

							for l_ind in range(len(wh_max)):
								l_max_real = wh_max[l_ind]	
								l_max_int = np.int(l_max_real)
								step_left = l_max_int - wh_min[l_ind]
								l_left_real = wh_min[l_ind]
								l_left_int = np.int(l_left_real)

								if(l_ind == len(wh_max) - 1 ):
									step_right = anglen - l_max_int + wh_min[0]
									l_right_real = wh_min[0]
									right_index = 0
								else:
									step_right =  wh_min[l_ind + 1] - l_max_int
									l_right_real = wh_min[l_ind + 1]
									right_index = l_ind + 1

								l_right_int = np.int(l_right_real)

								for l in range(1,int(step_right)):
									l_index = l_max_int + l
									if(l_index >= len(polygons_new)):
										l_index = l_index - len(polygons_new)
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])


								if( int(step_right) == step_right ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_right_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_right_int])

								for l in range(1,int(step_left)):
									l_index = l_max_int - l
									if( l_index < 0 ):
										l_index = len(polygons_new) + l_index
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])

								if( int(step_left) == step_left ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_left_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_left_int])


					wh_sum = wh_sum * hl_current * angstep / 360

					for l in wh_max:
						lint = np.int(l)
						if( wh_sum[lint] > 0 ):

							new_x = polygon[j][0] + polygons_new[lint] * cos((vec_ang[lint] + angstep*(l-lint) ) * np.pi / 180 ) *step_lon_deg/step_lon_m ; 
							new_y = polygon[j][1] + polygons_new[lint] * sin((vec_ang[lint] + angstep*(l-lint) ) * np.pi / 180 ) *step_lat_deg/step_lat_m ;
							height_eff = wh_sum[lint] + interpol_pos(lon1, lat1, step_lon_deg, step_lat_deg, new_x, new_y, cells_lon, cells_lat, Topography)
							if(interpol_pos(lon1, lat1, step_lon_deg, step_lat_deg, new_x, new_y, cells_lon, cells_lat, Topography) < 99999 ):
								polygon.append(( new_x, new_y, height_eff, polygon[j][3] + 1, l, wh_sum[lint] ))
			
			sum_pixels = sum(sum(data_step))	
			print((j, len(polygon), polygon[j][3], polygon[j][2], sum(sum(data_step)), polygon[j][4] ))

		if( N > 1 ):
			data_cones = data_cones + data_step

		print(' Simulation finished (N = ' + str(i+1) + ')')

if(topo_source == 'Upload_UTM'):
	
	data_cones = np.zeros((n_north,n_east))
	data_cones_save = np.zeros((n_north,n_east))
	data_aux_t = np.ones((n_north,n_east))
	data_aux_b = np.zeros((n_north,n_east))
	vec_ang = range(0, 360, angstep)

	for i in range(0,N):
		current_level = 0
		data_step = np.zeros((n_north,n_east))
		polygon = []
		height_eff = height_vector[i] + interpol_pos(east_cor, north_cor, cellsize, cellsize, east_cen_vector[i], north_cen_vector[i], n_east, n_north, Topography)
		polygon.append((east_cen_vector[i], north_cen_vector[i],  height_eff, 1.0, -1, height_vector[i] ))
		sum_pixels = 0
		hl_current = hl_vector[i]

		for j in range(10000): 
			if(j == len(polygon)):
				if( N == 1 ):			
					data_cones = data_cones + data_step
				break
			if( cone_levels < polygon[j][3] ):
				if( N == 1 ):
					data_cones = data_cones + data_step
				break
			elif(current_level < polygon[j][3]):
				current_level = polygon[j][3]
				if( N == 1 ):
					data_cones = data_cones + data_step

			polygon_xy = []
			polygons_new = []

			for angle_deg in vec_ang:
				angle_rad = angle_deg * np.pi /180
				for distance in range(0, 100000, distep):
					h = interpol_pos(east_cor, north_cor, cellsize, cellsize, polygon[j][0] + distance * cos(angle_rad) , polygon[j][1] + distance*sin(angle_rad) , n_east, n_north, Topography)
					if( h >= polygon[j][2] - hl_current * distance ):
						polygon_xy.append((int((polygon[j][0] + (distance-distep)* cos(angle_rad) - east_cor) * n_east / ( cellsize * ( n_east - 1 ) ) ), int((polygon[j][1] + (distance-distep)*sin(angle_rad) - north_cor) * n_north / ( cellsize * ( n_north - 1 ) ))))
						polygons_new.append(distance - distep)
						break						

			if( polygon[j][4] > -1 ):
				lim = np.int(polygon[j][4])
				if( polygon[j][4] == np.int(polygon[j][4]) ):
					for ii in range(anglen):
						vector_correc[ii] = vector_backward_1[int((ii - polygon[j][4] + index_max) % anglen)]

				else:
					for ii in range(anglen):
						vector_correc[ii] = vector_backward_2[int((ii - polygon[j][4] + index_max) % anglen)]

				polygons_new = polygons_new * vector_correc

			img = Image.new('L', (n_east, n_north), 0)
			if( len(polygon_xy) > 0 ):
				draw = ImageDraw.Draw(img).polygon(polygon_xy, outline = 0 , fill = 1)
				data_step = np.maximum( np.minimum(data_aux_t, data_step + np.array(img)), data_aux_b)

			if( cone_levels > polygon[j][3] and sum(sum(data_step)) > sum_pixels + pix_min ):
				aux = np.zeros(len(polygons_new) + 2) 
				aux[1:len(polygons_new)+1] = np.array(polygons_new) 
				aux[0] = polygons_new[len(polygons_new)-1]
				aux[len(polygons_new)+1] = polygons_new[0]
				der1 = (aux[1:len(aux)-1] - aux[2:len(aux)])
				der2 = (aux[1:len(aux)-1] - aux[0:len(aux)-2])
				wh1 = np.where(der1 >= 0)
				wh2 = np.where(der2 >= 0)
				wh_max = np.intersect1d(wh1[0], wh2[0])
				wh_grouped = np.split(wh_max, np.where(np.diff(wh_max) > 1)[0] + 1 )
				wh3 = np.where( abs(der1) > 0 )
				wh4 = np.where( abs(der2) > 0 )
				wh5 = np.intersect1d(wh_max, wh3[0])
				wh6 = np.intersect1d(wh_max, wh4[0])
				grouped_filter = np.zeros(len(wh_grouped))

				for x_grouped in range(len(wh_grouped)):
					if( len(np.intersect1d(wh_grouped[x_grouped],wh5)) > 0 and len(np.intersect1d(wh_grouped[x_grouped],wh6)) > 0):
						grouped_filter[x_grouped] = 1

				if( np.min(wh_grouped[0]) == 0 and np.max(wh_grouped[len(wh_grouped)-1]) == anglen - 1 ):
					if( len(np.intersect1d(wh_grouped[0],wh5)) > 0 and len(np.intersect1d(wh_grouped[len(wh_grouped)-1],wh6)) > 0):
						grouped_filter[len(wh_grouped) - 1] = 1
					aux_grouped = np.concatenate((wh_grouped[len(wh_grouped)-1], wh_grouped[0] + len(polygons_new)))
					aux_filter = grouped_filter[len(wh_grouped)-1] + grouped_filter[0]
					wh_grouped = wh_grouped[1:-1]
					wh_grouped.append(aux_grouped)
					grouped_filter = np.append(grouped_filter[1:-1],aux_filter)

				wh_max = []
				for k in range(len(grouped_filter)):
					if(grouped_filter[k] > 0 ):
						if(np.mean(wh_grouped[k]) < len(polygons_new) and np.mean(wh_grouped[k]) >= 0.0):
							wh_max.append(np.mean(wh_grouped[k]))
						elif( np.mean(wh_grouped[k]) < len(polygons_new) ):
							wh_max.append(len(polygons_new) + np.mean(wh_grouped[k]))
						else:
							wh_max.append(- len(polygons_new) + np.mean(wh_grouped[k]))

				wh1 = np.where(der1 <= 0)
				wh2 = np.where(der2 <= 0)
				wh_min = np.intersect1d(wh1[0], wh2[0])
				wh_grouped = np.split(wh_min, np.where(np.diff(wh_min) > 1)[0] + 1 )
				wh3 = np.where( abs(der1) > 0)
				wh4 = np.where( abs(der2) > 0)
				wh5 = np.intersect1d(wh_min, wh3[0])
				wh6 = np.intersect1d(wh_min, wh4[0])
				grouped_filter = np.zeros(len(wh_grouped))

				for x_grouped in range(len(wh_grouped)):
					if( len(np.intersect1d(wh_grouped[x_grouped],wh5)) > 0 and len(np.intersect1d(wh_grouped[x_grouped],wh6)) > 0):
						grouped_filter[x_grouped] = 1
					
				if( np.min(wh_grouped[0]) == 0 and np.max(wh_grouped[len(wh_grouped)-1]) == anglen - 1):
					if( len(np.intersect1d(wh_grouped[0],wh5)) > 0 and len(np.intersect1d(wh_grouped[len(wh_grouped)-1],wh6)) > 0):
						grouped_filter[len(wh_grouped) - 1] = 1
					aux_grouped = np.concatenate((wh_grouped[len(wh_grouped)-1], wh_grouped[0] + len(polygons_new)))
					aux_filter = grouped_filter[len(wh_grouped)-1] + grouped_filter[0]
					wh_grouped = wh_grouped[1:-1]
					wh_grouped.append(aux_grouped)
					grouped_filter = np.append(grouped_filter[1:-1],aux_filter)

				wh_min = []

				for k in range(len(grouped_filter)):
					if(grouped_filter[k] > 0 ):
						if(np.mean(wh_grouped[k]) < len(polygons_new) and np.mean(wh_grouped[k]) >= 0.0):
							wh_min.append(np.mean(wh_grouped[k]))
						elif(np.mean(wh_grouped[k]) < len(polygons_new) ):
							wh_min.append(len(polygons_new) + np.mean(wh_grouped[k]))
						else:
							wh_min.append(- len(polygons_new) + np.mean(wh_grouped[k]) )

				wh_sum = np.zeros(len(polygons_new)) 

				if(len(wh_max) > 0):
					
					if( len(wh_max) == 1):
						for l_max_real in wh_max:
							lmax = np.int(l_max_real)
							l_it = 	len(polygons_new) - 1		
							for l in range(1,len(polygons_new)):
								l_index = lmax + l
								if(l_index >= len(polygons_new)):
									l_index = l_index - len(polygons_new)
								if( polygons_new[lmax] < polygons_new[l_index] ):
									l_it = l
									break
								wh_sum[lmax] = wh_sum[lmax] + (polygons_new[lmax] - polygons_new[l_index])							
							for l in range(1,len(polygons_new) - l_it):
								l_index = lmax - l
								if(l_index < 0):
									l_index = l_index + len(polygons_new)
								if( polygons_new[lmax] < polygons_new[l_index] ):
									break							
								wh_sum[lmax] = wh_sum[lmax] + (polygons_new[lmax] - polygons_new[l_index])

					else:

						wh_max = np.sort(wh_max)
						wh_min = np.sort(wh_min)

						if(wh_min[0] > wh_max[0]):

							for l_ind in range(len(wh_max)):
								l_max_real = wh_max[l_ind]	
								l_max_int = np.int(l_max_real)
								step_right = wh_min[l_ind] - l_max_int
								l_right_real = wh_min[l_ind]
								l_right_int = np.int(l_right_real)

								if(l_ind == 0):
									step_left = anglen + l_max_int - wh_min[len(wh_min)-1]
									l_left_real = wh_min[len(wh_min) - 1]
									left_index = len(wh_min) - 1
								else:
									step_left = l_max_int - wh_min[l_ind - 1]
									l_left_real = wh_min[l_ind - 1]
									left_index = l_ind - 1
								
								l_left_int = np.int(l_left_real)

								for l in range(1,int(step_right)):
									l_index = l_max_int + l
									if(l_index >= len(polygons_new)):
										l_index = l_index - len(polygons_new)
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])

								if( int(step_right) == step_right ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_right_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_right_int])

								for l in range(1,int(step_left)):
									l_index = l_max_int - l
									if( l_index < 0 ):
										l_index = len(polygons_new) + l_index
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])

								if( int(step_left) == step_left ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_left_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_left_int])

						else:

							for l_ind in range(len(wh_max)):
								l_max_real = wh_max[l_ind]	
								l_max_int = np.int(l_max_real)
								step_left = l_max_int - wh_min[l_ind]
								l_left_real = wh_min[l_ind]
								l_left_int = np.int(l_left_real)

								if(l_ind == len(wh_max) - 1 ):
									step_right = anglen - l_max_int + wh_min[0]
									l_right_real = wh_min[0]
									right_index = 0
								else:
									step_right =  wh_min[l_ind + 1] - l_max_int
									l_right_real = wh_min[l_ind + 1]
									right_index = l_ind + 1

								l_right_int = np.int(l_right_real)

								for l in range(1,int(step_right)):
									l_index = l_max_int + l
									if(l_index >= len(polygons_new)):
										l_index = l_index - len(polygons_new)
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])


								if( int(step_right) == step_right ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_right_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_right_int])

								for l in range(1,int(step_left)):
									l_index = l_max_int - l
									if( l_index < 0 ):
										l_index = len(polygons_new) + l_index
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_index])

								if( int(step_left) == step_left ):
									wh_sum[l_max_int] = wh_sum[l_max_int] + 0.5 * (polygons_new[l_max_int] - polygons_new[l_left_int])
								else:
									wh_sum[l_max_int] = wh_sum[l_max_int] + (polygons_new[l_max_int] - polygons_new[l_left_int])

					wh_sum = wh_sum * hl_current * angstep / 360

					for l in wh_max:
						lint = np.int(l)
						if( wh_sum[lint] > 0 ):
							new_x = polygon[j][0] + polygons_new[lint] * cos((vec_ang[lint] + angstep*(l-lint) ) * np.pi / 180 ) ; 
							new_y = polygon[j][1] + polygons_new[lint] * sin((vec_ang[lint] + angstep*(l-lint) ) * np.pi / 180 ) ;
							height_eff = wh_sum[lint] + interpol_pos(east_cor, north_cor, cellsize, cellsize, new_x, new_y, n_east, n_north, Topography)
							if(interpol_pos(east_cor, north_cor, cellsize, cellsize, new_x, new_y, n_east, n_north, Topography) < 99999 ):
								polygon.append(( new_x, new_y, height_eff, polygon[j][3] + 1, l, wh_sum[lint] ))

			sum_pixels = sum(sum(data_step))	
			print((j, len(polygon), polygon[j][3], polygon[j][2], sum(sum(data_step)), polygon[j][4] ))

		if( N > 1 ):
			data_cones = data_cones + data_step

		print(' Simulation finished (N = ' + str(i+1) + ')')

if(topo_source == 'SRTM'  or topo_source == 'Upload_deg'):

	data_cones = data_cones[ range(len(data_cones[:,0]) -1 , -1 , -1 ) , : ] / N
	data_cones_save[:,:] = data_cones[:,:]
	line_val = data_cones.max()
	data_cones[data_cones[:,:] == 0] =  np.nan
	val_up = np.floor((line_val + 0.1 - 1.0 / N ) * 10.0) / 20.0
	val_down = np.maximum( val_up / 10.0 , 0.02 )
	plt.figure(1)
	cmapg = plt.cm.get_cmap('Greys')
	cmapr = plt.cm.get_cmap('Reds')
	cmaps = plt.cm.get_cmap('Blues') 

	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()

	if( N > 1 ):
		CS_Topo = plt.contourf(matrix_lon,matrix_lat,Topography, 100, alpha = 1.0, cmap = cmapg ,antialiased=True, lw=0.0001)
		CS_Sea = plt.contourf(matrix_lon,matrix_lat,Topography_Sea, 100, alpha = 0.5, cmap = cmaps ,antialiased=True, lw=100)
		CS = plt.contourf(matrix_lon, matrix_lat, data_cones, 100, vmin = 0.0, vmax = 1.0,  alpha= 0.3, interpolation='linear', cmap=cmapr, antialiased=True, lw=0.01)	
		fmt = '%.2f'
		plt.colorbar()
		CS_lines = plt.contour(matrix_lon,matrix_lat,data_cones, np.array([val_down, val_up]), colors='r', interpolation='linear', lw=0.01)
		plt.clabel(CS_lines, inline=0.1, fontsize = 7, colors='k', fmt=fmt)
	else:
		CS_Topo = plt.contourf(matrix_lon,matrix_lat,Topography, 100, alpha = 1.0, cmap = cmapg ,antialiased=True, lw=0.0001)
		CS_Sea = plt.contourf(matrix_lon,matrix_lat,Topography_Sea, 100, alpha = 0.5, cmap = cmaps ,antialiased=True, lw=100)
		CS = plt.contourf(matrix_lon, matrix_lat, data_cones, 100, alpha= 0.3, interpolation='nearest', cmap=cmapr, antialiased=True, lw=0.01)

	plt.axes().set_aspect(step_lat_m/step_lon_m)
	plt.xlabel('Longitude $[^\circ]$')
	plt.ylabel('Latitude $[^\circ]$')
	plt.xlim(lon1, lon2 )
	plt.ylim(lat1, lat2 )

	for i in range(0,N):
		plt.plot( lon_cen_vector[i], lat_cen_vector[i], 'r.', markersize=2)

	if( N == 1 ):
		for i in range(1,len(polygon)):
			plt.plot( polygon[i][0],polygon[i][1], 'b.', markersize=2)
	try:
		os.stat('Results')
	except:
		os.mkdir('Results')
	try:
		os.stat('Results/' + run_name)
	except:
		os.mkdir('Results/' + run_name)

	plt.savefig('Results/' + run_name + '/' + 'Map.png', bbox_inches = 'tight')

if(topo_source == 'Upload_UTM'):

	data_cones = data_cones[ range(len(data_cones[:,0]) -1 , -1 , -1 ) , : ] / N
	line_val = data_cones.max()
	data_cones_save[:,:] = data_cones[:,:]
	data_cones[data_cones[:,:] == 0] =  np.nan
	val_up = np.floor((line_val + 0.1 - 1.0 / N ) * 10.0) / 20.0
	val_down = np.maximum( val_up / 10.0 , 0.02 )
	plt.figure(1)

	cmapg = plt.cm.get_cmap('Greys')
	cmapr = plt.cm.get_cmap('Reds')
	cmaps = plt.cm.get_cmap('Blues') 

	mng = plt.get_current_fig_manager()
	mng.full_screen_toggle()

	if( N > 1 ):
		CS_Topo = plt.contourf(matrix_east / 1e3,matrix_north / 1e3,Topography, 100, alpha = 1.0, cmap = cmapg ,antialiased=True, lw=0.0001)
		CS_Sea = plt.contourf(matrix_east / 1e3,matrix_north / 1e3,Topography_Sea, 100, alpha = 0.5, cmap = cmaps ,antialiased=True, lw=100)
		CS = plt.contourf(matrix_east / 1e3,matrix_north / 1e3, data_cones, 100, vmin = 0.0, vmax = 1.0,  alpha= 0.3, interpolation='linear', cmap=cmapr, antialiased=True, lw=0.01)	
		fmt = '%.2f'
		plt.colorbar()
		CS_lines = plt.contour(matrix_east / 1e3 ,matrix_north / 1e3, data_cones, np.array([val_down, val_up]), colors='r', interpolation='linear', lw=0.01)
		plt.clabel(CS_lines, inline=0.1, fontsize = 7, colors='k', fmt=fmt)
	else:
		CS_Topo = plt.contourf(matrix_east / 1e3,matrix_north / 1e3,Topography, 100, alpha = 1.0, cmap = cmapg ,antialiased=True, lw=0.0001)
		CS_Sea = plt.contourf(matrix_east / 1e3,matrix_north / 1e3,Topography_Sea, 100, alpha = 0.5, cmap = cmaps ,antialiased=True, lw=100)
		CS = plt.contourf(matrix_east / 1e3 ,matrix_north / 1e3 ,data_cones, 100, alpha= 0.3, interpolation='nearest', cmap=cmapr, antialiased=True, lw=0.01)

	plt.axes().set_aspect(1.0)
	plt.xlabel('East [km]')
	plt.ylabel('North [km]')
	plt.xlim((east_cor)/1e3, (east_cor + cellsize * (n_east - 1))/1e3 )
	plt.ylim((north_cor)/1e3,(north_cor +cellsize * (n_north - 1))/1e3 )

	for i in range(0,N):
		plt.plot( east_cen_vector[i] / 1e3, north_cen_vector[i] / 1e3, 'r.', markersize = 2 )

	if( N == 1 ):
		for i in range(1,len(polygon)):
			plt.plot( polygon[i][0], polygon[i][1], 'b.', markersize = 2 )

	try:
		os.stat('Results')
	except:
		os.mkdir('Results')
	try:
		os.stat('Results/' + run_name)
	except:
		os.mkdir('Results/' + run_name)

	plt.savefig('Results/' + run_name + '/' + 'Map.png', bbox_inches = 'tight')

log =  log + '\n'   + 'The probability map was saved in Results/' + run_name + '/Map.png'

if(N > 1 and var_height >= 1E-5):
	log =  log + '\n' + 'The histogram of the distribution of collapse height is presented in Histogram Height'

if(N > 1 and var_hl >= 1E-5):
	log =  log + '\n'  + 'The histogram of the distribution of H/L is presented in Histogram H/L'
	log =  log + '\n' + 'Topography data was saved in Results/' + run_name + '/Topography.txt'
	log =  log + '\n' + 'Probability data was saved in Results/' + run_name + '/Probability.txt'
	log =  log + '\n' + 'Mesh data was saved in Results/' + run_name + '/Mesh_East.txt and Results/' + run_name + '/Mesh_North.txt'

j = 0

lib.put('output.log(multi'+str(j)+').about.label','Readme')
lib.put('output.log(multi'+str(j)+')', log)

if(var_height >= 1E-5):
	j = j+1
	aux_h = np.histogram(height_vector,min(10,N))
	xy1 = []
	for i in range(len(aux_h[0])):
		str1 = str(0.5*(aux_h[1][i]+aux_h[1][i+1]))+' '+str(aux_h[0][i]/float(N))
		str1 = str1.replace(",", "")
		str1 = str1.replace("[", "")
		str1 = str1.replace("]", "")
 		xy1.append(str1)
	xy1 =  "\n".join(xy1)

	lib.put('output.histogram(multi'+str(j)+').about.group','Histogram Height')
	lib.put('output.histogram(multi'+str(j)+').xaxis.label','Collapse Height [m]')
	lib.put('output.histogram(multi'+str(j)+').yaxis.label','Frequency')
	lib.put('output.histogram(multi'+str(j)+').component.xy',xy1)

if(var_hl >= 1E-5):
	j = j+1
	aux_h = np.histogram(hl_vector,min(10,N))
	xy1 = []
	for i in range(len(aux_h[0])):
		str1 = str(0.5*(aux_h[1][i]+aux_h[1][i+1]))+' '+str(aux_h[0][i]/float(N))
		str1 = str1.replace(",", "")
		str1 = str1.replace("[", "")
		str1 = str1.replace("]", "")
 		xy1.append(str1)
	xy1 =  "\n".join(xy1)

	lib.put('output.histogram(multi'+str(j)+').about.group','Histogram H/L')
	lib.put('output.histogram(multi'+str(j)+').xaxis.label','H/L')
	lib.put('output.histogram(multi'+str(j)+').yaxis.label','Frequency')
	lib.put('output.histogram(multi'+str(j)+').component.xy',xy1)


np.savetxt('Results/' + run_name + '/Topography.txt', Topography_Save, fmt='%.2e')
np.savetxt('Results/' + run_name + '/Probability.txt', data_cones_save, fmt='%.2e')

if(topo_source == 'SRTM' or topo_source == 'Upload_deg'):
	np.savetxt('Results/' + run_name + '/Mesh_East.txt', matrix_lon, fmt='%.2e')
	np.savetxt('Results/' + run_name + '/Mesh_North.txt', matrix_lat, fmt='%.2e')
if(topo_source == 'Upload_UTM'):
	np.savetxt('Results/' + run_name + '/Mesh_East.txt', matrix_east, fmt='%.2e')
	np.savetxt('Results/' + run_name + '/Mesh_North.txt', matrix_north, fmt='%.2e')

Rappture.result(lib)

sys.exit(0)
