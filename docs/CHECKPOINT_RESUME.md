# 手动保存和续算

此接口仅适用于支持`--resume`的新程序。正在运行的旧程序不会因源码更新而获得该功能。

新任务使用原来的计算参数启动即可，程序默认写`mc_checkpoint.control.json`。
可选`--checkpoint-every 100`每100sweep保存；默认只在收到信号时保存。

在仓库根目录请求20022上的**一个**子任务保存后退出：

```bash
python3 scripts/checkpoint_job.py \
  --port 20022 --job 'JOB_ID[INDEX].mgt' \
  --control '/home/sunxr/RUN_ROOT/runs/CHAIN/mc_checkpoint.control.json' \
  --stop
```

去掉`--stop`则保存后继续。20015改相应端口即可。
脚本核对作业、节点与PID身份，等待`checkpoint_saved: true`才返回成功；
不会取消PBS任务。一个sweep可能持续数分钟，可用`--timeout 1800`延长等待。
若节点SSH不允许或程序不支持接口，脚本报错且不改作业。
不要直接对Bash/Python包装器`qsig`，也不要在未确认保存前`qdel`。

有执行节点访问权限时，也可向control中记录并核实身份的**MC PID**发送
`kill -USR1 PID`（保存继续）或`kill -USR2 PID`（保存退出），随后检查receipt。
保存只能在当前完整sweep及其测量结束后发生；SIGKILL无法保存。
SIGTERM只是尽力保护，PBS随后强杀时未必来得及完成当前sweep。

用原来**同一可执行文件和库环境**在新PBS输出目录续算：

```bash
/path/to/original/runtime/kitaev_ssh_mc \
  --resume /absolute/old/path/mc_checkpoint.dat \
  --output /absolute/new/path/mc
```

原参数、warmup进度、RNG及统计自动继承，不需要重复填写。可传
`--sweeps 30000`把全链总measurement目标改为30000；这不表示再加30000。
恢复不重新热化、不重新随机抽初态，不允许改lambda、dt、更新方式或测量集合。
新目录须已存在，输出前缀必须全新，避免覆盖历史文件。

保存后退出码为75，results状态`checkpointed`；包装器需识别这个码，
不能宣称生产完成。续算正常结束才输出`complete`。

恢复后的`mc_blocks.csv`已经包含检查点以前的全部完整块，未满块也会继续补齐。
分析时用最新续算结果代表该seed，**不要再与旧blocks相加**。
诊断CSV和主运行时间只描述本次进程，sweep编号接续。
保存了RNG并不保证一般相互作用链跨缓存重建逐比特重复；详细范围见
[ADR-027](decisions/ADR-027-signal-checkpoint-resume.md)。
