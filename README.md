# SIMO-Package

## Introduction

SIMO-Package is a basic NURBS-based IsoGeometric Analysis package written in MATLAB.

## Modules

* Rotation-free beam theory (Kirchhoff-Love beam theory)
* Rotation-free Kirchhoff-Love plate
* Linear Static Elasticity (including 1D, 2D, 3D elements for single patch and muti-patches problems)
* Gradient Elasticity (we are still working on it)
* IGLA (IsoGeometric Limit Analysis)
* Linear Static Heat Tranfer

## Dependencies
* Paraview (http://www.paraview.org/) for visualization

## Active Developers

* Hung Nguyen-Xuan
* Khanh Chau-Nguyen

## Contributors

* Tuan Nguyen-Ngoc
* Son Thai-Bao
* Hien Van Do
* Hoang Xuan Nguyen

## HSwang，Frank Wang，Max Zhang
-  We add annotations to code together
-  我们一起给代码添加中文注释

## 中文注释乱码的问题

由于Matlab采用的是ANSI编码，一般编辑器中文编码都是UTF-8。我们尽量使用ANSI，但难免会不小心使用UTF-8编码，如果你使用Matlab打开脚本文件发现中文注释是乱码的，只需要使用Notepad++或者其他文本编辑器将其编码转成ANSI就可以正常的在Matlab里阅读注释了。