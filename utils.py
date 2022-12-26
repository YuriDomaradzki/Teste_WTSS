import shapely
import geopandas as gpd
from datetime import datetime


def sub_list(ln_pnts, step):
    for i in range(0, len(ln_pnts), step):
        yield ln_pnts[i:i+step]
    
def generating_points(df_area: gpd.GeoDataFrame, grid: tuple):
    bounds = df_area['geometry'].bounds.iloc[0]
    res_x, res_y = grid
    points = []
        
    bounds = (float(coord) for coord in bounds)
    x_min, y_min, x_max, y_max = bounds

    resolution_x = (x_max - x_min) / res_x
    resolution_y = (y_max - y_min) / res_y

    for py in range(0, res_y, 1):
        for px in range(0, res_x, 1):
            x1 = x_min + px * resolution_x
            y1 = y_min + py * resolution_y

            point = shapely.Point((x1, y1))
                
            if df_area.contains(point).bool():
                points.append(point)

    points = {'ID': list(range(0, len(points))), 'geometry': points}
    df_samples = gpd.GeoDataFrame(points, crs='EPSG:4326')

    return df_samples 

def to_datetime(date):
    return datetime.strptime(date, "%Y-%m-%d")
