# Closure

本项目旨在研究分布式散射体（Distributed Scatters, DS）研究中的闭合环残差问题。Sentinel-1数据已经过ISCE软件处理。

## 结构

源码介绍：

- `functions.py`: 功能性的函数，包括map、reduce的封装接口，文件路径获取函数，样本区域获取函数。
- `s1.py`：读取Sentinel-1数据，对影像进行burst合并、裁切等基本处理。
- `landcover.py`: 展示样本区域地表覆盖情况，以及裁切出来的样本区域状况。
- `closure.py`: 实验结果。

## 问题

Jupyter Notebook与Python的multiprocessing存在冲突。我的程序在Linux上编写运行，表现良好；但在Windows和Macos上测试都无法正常运行。临时解决办法是将函数复制到`.py`文件中运行。