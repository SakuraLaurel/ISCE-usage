本节包括以下内容：

1. Sentinel-1、TerraSAR-X、ALOS-2影像的批量处理。
2. ISCE的干涉、多视、相干性计算三个环节的复现。
3. SAR影像的间接几何定位，不含误差校正内容。

其中

- xmlGenerate文件夹下是**内容1**的实现。给定主影像和其他影像的路径后，它会搜索文件夹下的所有影像数据，然后为每个辅影像单独生成一个bash脚本，最后令python将脚本全部执行即可。
- `src/interferometry.py`是对**内容2**的实现。
- `src/geolocate.py`是间接几何定位的实现。