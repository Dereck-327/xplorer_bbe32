# BBE32_EP Simulate

## 快速导航 (Quick Links)

* [**版本发布说明 (Release Notes)**](release.md)

## clangd 跳转使用说明

`.clangd` 文件中关于项目的路径需要替换为本机的路径. 编译器的路径如果和文件中一样则不需要修改
`.clangd-tie-shim.h`是clangd需要忽略的警告和检查错误,与编译无关

```.clangd
    ***
    - D:/WorkProgram/21.DSP_Simulate/xplorer_bbe32/.clangd-tie-shim.h
    - -ID:/WorkProgram/21.DSP_Simulate/xplorer_bbe32/bbe32ep_library/include
    - -ID:/WorkProgram/21.DSP_Simulate/xplorer_bbe32/bbe32ep_library/include_private
    ***
```

## `bbe32ep_library` 的编译

在`xplorer`中需要设置`bbe32ep_library` 项目去忽略掉一些源文件, 需要忽略的如下
![exclude](tools\xplorer_setting\exclude.png)