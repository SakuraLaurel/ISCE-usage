# 程序安装

推荐使用conda安装。安装时可能遇到依赖解析耗时极久、甚至失败的问题，推荐使用libmamba来解析依赖。

```
conda install isce2 -c conda-forge --solver libmamba
```

开发者说GitHub版的isce2缺少部分组件，无法使用较老的insarApp和isceApp接口。由于stripmapApp和topsApp两个新版接口已经覆盖了insar处理的绝大部分功能，因此也没有必要再使用老式接口。

