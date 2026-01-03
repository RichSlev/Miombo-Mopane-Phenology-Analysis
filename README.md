# Miombo-Mopane-Phenology-Analysis

## Data overview

This repository uses a combination of land-cover, vegetation index, and climate datasets to analyse phenological dynamics across Miombo and Mopane woodlands in north-western Zimbabwe, with a focus on spatial heterogeneity across land-use and protection gradients.

---

## Raster datasets

### Land use / land cover (LULC)

#### `LULC_MCC.tif`
Land use / land cover map derived from Sentinel-2 imagery using the IOTA² processing chain.  
Spatial resolution: 10 m.

Coverage includes:
- Mucheni Community Conservancy
- Chizarira
- Chete
- Sijaria

This product distinguishes multiple woodland structural types.

**Class legend**
- 1: Water  
- 2: Cropland  
- 3: Bare soil  
- 41: Miombo shrub (<3 m height)  
- 42: Miombo open forest (>3 m height, <70% canopy cover)  
- 43: Miombo closed forest (>3 m height, >70% canopy cover)  
- 51: Mopane shrub  
- 52: Mopane open forest  
- 53: Mopane closed forest  
- 6: Grassland  
- 7: Riparian woodland  
- 10: Settlements  

---

#### `Detailed_LULC_Kariba_area.tif`
Land use / land cover map for the Kariba region, derived from Sentinel-2 imagery using the IOTA² processing chain.  
Spatial resolution: 10 m.

This product provides a coarser thematic classification than `LULC_MCC.tif` and is used for regional-scale stratification.

**Class legend**
- 1: Water  
- 2: Cropland  
- 3: Bare soil  
- 4: Miombo woodland  
- 5: Mopane woodland  
- 6: Grassland  
- 7: Riparian woodland  
- 10: Settlements  

---

### Climate data

#### `CHIRPS_pentad.tif`
Rainfall dataset derived from CHIRPS via Google Earth Engine.  
- Temporal resolution: pentadal (5-day)  
- Spatial resolution: 0.05°  
- Units: millimetres per pentad  
- Structure: one raster band per pentad date  

Used to characterise seasonal rainfall dynamics and interannual variability.

---

### Vegetation indices

#### `NDVI.tif`
NDVI time series derived from Google Earth Engine.  
- Temporal resolution: 16-day composites  
- Spatial resolution: 250 m  
- Value range: −2000 to 10000 (scaled integer NDVI)  
- Structure: one raster band per acquisition date  

Used for phenological metric extraction (e.g. season onset, amplitude, duration).

---

## Vector datasets

### Shapefiles

- `Study_zone_Kariba_lake.shp`  
  Study extent encompassing national parks, state forests, safari areas, and communal lands surrounding Lake Kariba.

- `MCC.shp`  
  Boundary of Mucheni Community Conservancy.

- `Protected_areas.shp`  
  Protected area boundaries sourced from the World Database on Protected Areas (WDPA).

---

## Notes on workflow

For exploratory and method development phases, it is recommended to:
- begin analyses within the combined MCC + Chizarira + Chete + Sijaria extent  
- validate phenological patterns at local scale before extending to the full Kariba region  

This approach limits confounding effects of coarse spatial resolution and heterogeneous land-use mosaics during early analysis stages.
