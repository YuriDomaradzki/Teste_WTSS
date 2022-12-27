import glob
import folium
import webbrowser
import numpy as np
import pandas as pd
import multiprocessing
import geopandas as gpd
import concurrent.futures
from wtss import *
from datetime import datetime, timedelta
from utils import sub_list, generating_points, to_datetime
from wtss.utils import to_datetime as bdc_to_datetime

class WTSS_Harmonize:
    
    def __init__(self, token):
        WTSS_Harmonize._service = WTSS('https://brazildatacube.dpi.inpe.br/', access_token=token)

    @classmethod
    def _get_coverage(cls):
        return cls._coverage

    @classmethod
    def _set_coverage(cls, name):
        cls._coverage = cls._service[name]

    @classmethod
    def _get_region_info(cls):
        return cls._region_info

    @classmethod
    def _set_region_info(cls, info: tuple):
        cls._region_info = info
        
    @property
    def list_coverages(cls):
        return cls._service._list_coverages()
    
    @classmethod
    def list_coverage_attributes(cls, name):
        return [coverage['name'] for coverage in cls._service._describe_coverage(name)['attributes']]

    @classmethod
    def _query_wtss(cls, items_process: dict):
        df_ts = gpd.GeoDataFrame()
        for point in items_process['points']:
            ts = cls._get_coverage().ts(attributes=items_process['attributes'],
                            latitude=float(point.y), longitude=float(point.x),
                            start_date=items_process['start_date'], end_date=items_process['end_date'])
            timeline = bdc_to_datetime([date[0:10] for date in ts.timeline])
            geometry = np.full(len(timeline), point)
            points = np.full(len(timeline), '{:.5f} {:.5f}'.format(point.x, point.y))
            lat = np.full(len(timeline), '{:.5f}'.format(point.x))
            long = np.full(len(timeline), '{:.5f}'.format(point.y))

            tl = {'timeline': timeline}
            tl.update({chave: ts.values(chave) for chave in ts.attributes})
            tl.update({'point': points, 'lat': lat, 'long': long, 'geometry':geometry})

            df = gpd.GeoDataFrame(tl, crs="EPSG:4326")
            df_ts = pd.concat([df_ts, df])
        
        df_ts.reset_index(inplace=True, drop=True)
        return df_ts

    @classmethod      
    def get_timeseries(cls, municipio: str, uf: str, start_date: str, end_date: str, 
                       coverage: str, attributes: tuple, grid: tuple):
        
        cls._set_coverage(coverage)
        cls._set_region_info((uf.upper(), municipio.title()))

        # open shp according to federative unit and filter by municipality
        file = glob.glob(f'shapefile/{uf.upper()}/*.shp')[0]
        df_municipio = gpd.read_file(file)
        df_municipio = df_municipio.loc[df_municipio['NM_MUN'] == municipio.title()].set_crs("EPSG:4326", allow_override=True)   
        points_interior = generating_points(df_municipio, grid)['geometry']
        
        items_process = []
        step = len(points_interior) // multiprocessing.cpu_count()
        for points in list(sub_list(points_interior, step)):
            items_process.append({'points': points,
                                'start_date': start_date,
                                'end_date': end_date,
                                'coverage': coverage,
                                'attributes': attributes})
            
        df_ts = gpd.GeoDataFrame()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            mapped = executor.map(cls._query_wtss, items_process)
            
            for df in mapped:
                df_ts = pd.concat([df_ts, df])
        
        df_ts.reset_index(inplace=True, drop=True)
        return df_ts

    @classmethod
    def plot(cls, df: gpd.GeoDataFrame, attributes: list, aggregation: str, start_date: str,  
             end_date=None, qtd_days=None, save_path=None):

        uf, municipio = cls._get_region_info()

        try:
            start_date = to_datetime(start_date).date()
        except ValueError:
            return 'Insert a start_date in the YYYY-mm-dd format'
        
        if not end_date and not qtd_days:
            return "It's necessary to define an end_date or qtd_days" 
        
        if not end_date:
            end_date = (start_date + timedelta(days=qtd_days))
        
        if start_date < df['timeline'][0]:
            return 'Start_date is outside the search range done using WTSS'
        
        if end_date > df['timeline'][len(df)-1]:
            return 'End_date is outside the search range done using WTSS'
        
        file = glob.glob(f'shapefile/{uf}/*.shp')[0]
        df_municipio = gpd.read_file(file)
        df_municipio = df_municipio.loc[df_municipio['NM_MUN'] == 
                                        municipio].set_crs('EPSG:4326',
                                                                   allow_override=True)
        centroid_x, centroid_y = df_municipio['geometry'].iloc[0].centroid.xy
        
        for attribute in attributes:

            m = folium.Map([centroid_y[0], centroid_x[0]], tiles='Cartodb Positron',     
                            zoom_start=10, width='%100', height='%90', overlay=True, 
                            control=True, show=True)
            
            # Create municipality Shape Layer
            for _, shape in df_municipio.iterrows():
                sim_geo = gpd.GeoSeries(shape['geometry']).simplify(tolerance=0.001)
                geo_j = sim_geo.to_json()
                geo_j = folium.GeoJson(data=geo_j,
                                    name=f'Limits of {municipio}',
                                    style_function=lambda x: {'color': 'black',
                                                            'opacity': 1,
                                                            'fillColor': 'white',
                                                            'fillOpacity': 1})
                geo_j.add_to(m)
            
            df_time = df.loc[(start_date <= df.timeline) & (df.timeline <= end_date)]
            df_aggregate = df_time.dissolve(by='point', aggfunc=aggregation, as_index=False)

            df_aggregate.explore(column=attribute, 
                                    k=10, 
                                    m=m,
                                    legend_kwds=dict(colorbar=True), 
                                    tooltip=['point', attribute],
                                    tooltip_kwds=dict(labels=True),
                                    name=attribute)
            folium.LayerControl(collapsed=False).add_to(m)
            
            name_coverage = cls._get_coverage().name
            loc = f"{name_coverage} {aggregation.title()} {attribute.upper()} between {start_date}/{end_date}".upper()
            title_html = '''
                        <h3 align="center" style="font-size:16px"><b>{}</b></h3>
                        '''.format(loc)
            m.get_root().html.add_child(folium.Element(title_html))
            
            path = f'{save_path}/{name_coverage}_{attribute}.html' if save_path else f'Maps/{name_coverage}_{attribute}.html' 
            m.save(path)
            webbrowser.open(path)


if __name__ == "__main__":
    service = WTSS_Harmonize(token='change_me')

    start=datetime.now()
    x, y = 10, 10
    df_ts = service.get_timeseries(municipio='cametá', uf='pa', start_date='2020-01-01', end_date='2020-12-31', 
                               coverage='MOD13Q1-6', attributes=('EVI', 'NDVI'), grid=(x, y))
    print(datetime.now()-start)
    
    service.plot(df=df_ts, attributes=['EVI', 'NDVI'], aggregation='mean', start_date='2020-01-01', qtd_days=60)
