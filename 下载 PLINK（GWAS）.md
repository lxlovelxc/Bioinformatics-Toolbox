# 下载 PLINK v1.9（适合绝大多数 GWAS）
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip
unzip plink_linux_x86_64_20220402.zip
sudo mv plink /usr/local/bin/
sudo chmod +x /usr/local/bin/plink


plink --version


# 下载 GCTA 工具包（包含 GSMR）版本 1.94.4，适用于 Linux（内核 3.x 或更高）
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.4-linux-kernel-3-x86_64.zip

# 解压下载的 ZIP 包
unzip gcta-1.94.4-linux-kernel-3-x86_64.zip

# 进入解压后的目录
cd gcta-1.94.4-linux-kernel-3-x86_64/

# 查看目录内容，确认可执行文件是否存在
ls

# 为 gcta64 设置可执行权限（如果还未设置）
chmod +x gcta64

# 测试运行 GCTA，看是否正常输出帮助信息（包括 GSMR 参数）
./gcta64

# 将 gcta64 安装到系统路径中，便于在任何目录下直接调用
sudo cp gcta64 /usr/local/bin/

# 回到主目录
cd ~

# 再次运行 gcta64，确认可以在任意路径下使用（说明安装成功）
gcta64

