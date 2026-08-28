# NetMHCpan 安装包目录

NetMHCpan **不能**通过 pip/conda 直接安装，需从 DTU 获取：

1. 注册：https://services.healthtech.dtu.dk/cgi-bin/request.cgi?tool_id=NetMHCpan  
2. 邮件中的链接下载 Linux 包，例如 `netMHCpan-4.1b.linux.tar.gz`  
3. 将文件放到本目录后执行：

```bash
cd /path/to/neoag_event_pipeline
bash scripts/install_netmhcpan.sh vendor/netMHCpan-4.1b.linux.tar.gz
```

`data.tar.gz` 会在安装时自动从 DTU 下载（若尚未解压）。
