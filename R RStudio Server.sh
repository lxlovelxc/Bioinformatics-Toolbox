#!/bin/bash
# R + RStudio Server 安装与管理脚本（适用于 Ubuntu 20.04 / 22.04 / 24.04）

# ========================
# 1. 安装 R 语言环境
# ========================
sudo apt update
sudo apt install -y r-base

# ========================
# 2. 下载 RStudio Server 安装包（适用于 Ubuntu 系统）
# ========================
wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2025.05.1-513-amd64.deb

# ========================
# 3. 安装 gdebi（推荐的 .deb 安装工具，可自动处理依赖）
# ========================
sudo apt install -y gdebi-core

# ========================
# 4. 使用 gdebi 安装 RStudio Server
# ========================
sudo gdebi rstudio-server-2025.05.1-513-amd64.deb

# ========================
# 5. 设置系统用户以登录 RStudio Web 界面
# ========================
# （可选）新建一个名为 luoxi 的用户，用于登录 RStudio Server
sudo adduser luoxi

# ========================
# 6. 控制 RStudio Server 服务的启动与停止
# ========================

# 停止 RStudio Server（关闭 Web 访问）
sudo systemctl stop rstudio-server

# 启动 RStudio Server（开启 Web 访问）
sudo systemctl start rstudio-server

# 再次停止（例如调试或关闭服务）
sudo systemctl stop rstudio-server

# ========================
# 7. 登录 RStudio Server 的地址
# ========================
# 打开浏览器并访问：http://127.0.0.1:8787
# 使用你设置的系统用户名和密码登录

