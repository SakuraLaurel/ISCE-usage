# stripmapApp.py解析

执行`stripmapApp.py --steps --help`观察其标准流程：

```bash
...
['startup', 'preprocess', 'cropraw', 'formslc', 'cropslc']
['verifyDEM', 'topo', 'geo2rdr', 'coarse_resample', 'misregistration']
['refined_resample', 'dense_offsets', 'rubber_sheet_range', 'rubber_sheet_azimuth', 'fine_resample']
['split_range_spectrum', 'sub_band_resample', 'interferogram', 'sub_band_interferogram', 'filter']
['filter_low_band', 'filter_high_band', 'unwrap', 'unwrap_low_band', 'unwrap_high_band']
['ionosphere', 'geocode', 'geocodeoffsets', 'endup']
...
```

再观察`stripmapApp.py`文件内容，找到标准流程对应的命令：

```python
...
self.timeStart = time.time()
self.help()
# Run a preprocessor for the two sets of frames
self.runPreprocessor()
#Crop raw data if desired
self.runCrop(True)
self.runFormSLC()
self.runCrop(False)
#Verify whether user defined  a dem component.  If not, then download
# SRTM DEM.
self.verifyDEM()
# run topo (mapping from radar to geo coordinates)
self.runTopo()
# run geo2rdr (mapping from geo to radar coordinates)
self.runGeo2rdr()
# resampling using only geometry offsets
self.runResampleSlc('coarse')
# refine geometry offsets using offsets computed by cross correlation
self.runRefineSecondaryTiming()
# resampling using refined offsets
self.runResampleSlc('refined')
# run dense offsets
self.runDenseOffsets()
...
```

虽然标准流程中没有出现`self.startup()`这一命令，但找到`startup`函数

```python
def startup(self):
    self.help()
    self._insar.timeStart = time.time()
```

可以发现它其实就是标准流程中的前两行。其余步骤的功能大致如下：

- preprocess：通过`applications/make_raw.py`中的函数真正解析SAR影像的内容。主程序`stripmapApp.py`其实只是一个整合包，preprocess的实现位于`$ISCE_HOME/components/isceobj/StripmapProc/runPreprocessor.py`中，且之后大多数步骤的实现都位于这个路径下，只是文件名可能和步骤名不一致；少数步骤的实现在其它模块内。
- cropraw：对原始影像进行裁剪。如果原始影像不是raw而是slc，并且没有regionOfInterest，那么这个步骤会被直接跳过。
- formslc：将raw数据聚焦为slc数据，函数的实现在`runROI.py`中。同样由于现在大多数影像都是slc影像而常常被跳过。
- cropslc：调用的和cropraw其实是一个函数，主要处理regionOfInterest。
- verifyDEM：检查DEM数据是否有效，如果不符合要求就自己下载。下载使用的是`components/contrib/demUtils/DemStitcher.py`中的代码。注意如果一直连接不到官方网站，最后下载下来的DEM会是空白的，最好用mdx检查一下。
- topo：计算主影像上每个像元投影到WGS84坐标系下的位置，最后调用的是`components/zerodop/topozero/Topozero.py`中的代码。
- geo2rdr：将主影像WGS84坐标投影到辅影像的SAR坐标系下，调用的是`components/zerodop/geo2rdr/Geo2rdr.py`。最后的结果是offsets的形式，方便辅影像直接采样到主影像坐标系下。由于这里仅仅通过几何定位来实现配准，因此配准结果可能是不精确的，称下一步为coarse采样。
- coarse_resample：将辅影像采样到主影像的SAR坐标系下。
- misregistration：这一步骤将主影像和采样后的辅影像进行配准。可能是因为代码直接来源于其他程序的DEM强度模拟影像与SAR影像配准步骤，相关代码注释说是simulated影像与SAR影像的配准。这是不准确的，所有影像都是真实的SAR影像，而且是复数影像。这一步骤实际调用的函数来自`components/isceobj/StripmapProc/runRefineSecondaryTiming.py`。

对于当前案例，`stripmapApp.py`的步骤研究到这里就足够了，已经找到了两幅影像的强度匹配程序在哪里。

# 配准处理

如果不依赖`stripmapApp.py`而是自行利用ISCE的配准功能对两幅SAR影像进行配准，可以参考<a href="./TestAmpcor.py">TestAmpcor.py</a>，程序参照了`components/isceobj/StripmapProc/runRefineSecondaryTiming.py`的写法，对其中一些参数进行了固定。只运行这一个函数，即可将采样后的辅影像与主影像进行配准，配准结果与`stripmapApp.py`计算结果一致。

作为对照，还将原始的辅影像与主影像进行了配准，但是程序显示强度匹配的效果很差，最后计算结果也远远偏离offsets的计算值。