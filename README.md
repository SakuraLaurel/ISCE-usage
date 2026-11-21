这是星载SAR影像几何定位精度研究的相关代码。论文地址：[High-Resolution Spaceborne SAR Geolocation Accuracy Analysis and Error Correction](https://www.mdpi.com/2072-4292/16/22/4210)

# 1 介绍

本仓库包含对流层校正、距离-多普勒方程解算、TerraSAR-X和天绘二号数据读取的实现，也有电离层延迟等的少许相关代码。因为isce2的环境要求比较苛刻，所以本项目中与isce2处理相关的内容放在了[另一个仓库](https://github.com/SakuraLaurel/ISCE-usage)。

# 2 结构

- `functions`: 包含了算法逻辑代码
  - `dem.py`: 与DEM读取和应用有关的代码
  - `image.py`: 与TerraSAR-X、天绘二号、Sentinel-1读取有关的代码
  - `trans.py`: 与时间、坐标系、方向换算有关的代码
  - `wrf.py`: 与WRF数据处理结果的应用有关的代码
- `src`： 包含了代码示例
  - `1_read_tx.ipynb`: TerraSAR-X数据的几何定位、对流层校正有关的代码示例
  - `2_read_th.ipynb`: 天绘二号数据的几何定位、对流层校正有关的代码示例
  - `3_coregistration.ipynb`: 使用isce2将多景TerraSAR-X配准，对配准结果的统计
  - `4_troposphere.ipynb`: 对Ray tracing、Zenith delay mapping方法和gacos效果的比较
  - `5_ionosphere.ipynb`: 电离层校正相关内容，因为对X波段的SAR影像影响过于不显著，所以未加入`functions`中
  - `6_dem.ipynb`: 对DEM高程精度与几何定位误差间的关系的探索
  - `7_others.ipynb`: Sentinel-1影像的读取和定位

本项目中不包含`data`文件夹，即用于实验的数据，因为数据实在太大。