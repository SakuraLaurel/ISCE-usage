我是Arch linux用户，因此后文都按照linux的使用习惯来展开介绍。理解本项目的内容需要一定的linux基础。

# 程序安装

推荐使用conda安装。安装时可能遇到依赖解析耗时极久、甚至失败的问题，推荐使用libmamba来解析依赖。

```bash
conda install isce2 -c conda-forge --solver libmamba
```

开发者说GitHub版的isce2缺少部分组件，无法使用较老的insarApp和isceApp接口。由于stripmapApp和topsApp两个新版接口已经覆盖了insar处理的绝大部分功能，因此也没有必要再使用老式接口。

因为stripmapApp等高级接口在applications下，显示组件mdx在bin下，为了方便使用，可以将isce主目录及相关路径加到环境变量中。

```bash
# 我的环境名为isce，python版本是3.11
export ENV_PATH=~/.conda/envs/isce
export ISCE_HOME=$ENV_PATH/lib/python3.11/site-packages/isce
export PATH=$PATH:$ISCE_HOME/applications:$ISCE_HOME/bin
```

# 基本结构

进入到`$ISCE_HOME`路径下，可以看到applications、components等文件夹，其中最重要的是这两个文件夹。
- applications是对components组件的高级封装，包含了相对完整的处理流程。如果需要用ISCE进行完整的干涉处理，那就要用到这个文件夹下的程序。

  InSAR的处理流程可以参考<a href="./isce_tutorial.pdf">isce_tutorial.pdf</a>这个文件。本项目内的3份pdf文件都比较老了，主要记录了insarApp.py时期的使用方法，但是ISCE的使用逻辑一直没有改变，用来理解ISCE的工作流程是足够的。如果对insarApp.py每一步骤的作用感兴趣，那么可以参考<a href="./ISCE.pdf">ISCE.pdf</a>。

- components内有许多功能模块，如果需要利用ISCE完成相对单一的功能，就需要自己调用这里的模块。理解模块的功能可以参考<a href="./ISCE_modules_Piyush.pdf">ISCE_modules_Piyush.pdf</a>。