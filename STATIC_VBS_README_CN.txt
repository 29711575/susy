静态VBS软模工具（quick1010，2026-10-10）

当前 note 已更正，零场模长峰任务见 VBS_AMPLITUDE_README_CN.txt 与 ADR-019。
以下工具按旧 note 完成，保留作独立静态参考；以下“note参数”均指旧版本。

从/home/asus/newsusy执行：
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR=/tmp/quick1010-static-mpl \
python3 quick1010/kitaev_ssh/static_vbs_response.py \
  --output-dir fast-1006/audit/1010/YOUR_NEW_TASK \
  --Jx 1 --Jy 1 --Jz .5 --g .3 --mass 1.7 --omega-ph .8 \
  --spring-k 1.088 --uniform-grid 800 --reference-grids 200,400 \
  --geometries '6,6,0;12,3,0;12,3,3' --finite-step .001 --plot

选择新的输出目录，不覆盖历史结果。--spring-k是m*omega^2一致性检查。
不指定上述物理选项时默认采用note.docx中的静态演示点；不是生产QMC默认参数。
800均匀背景共同用于所有几何；另计算400背景检查。有限体系扫描所有非零允许q。
JSON保存完整三键cos/sin矩阵、软子空间、三档真实差分和两档同q积分。
CSV为全q摘要；--plot只导出PNG；q_candidates文件含允许测量波矢及其共轭。
候选文件用于--vbs-q-file和selected/both等时测量，susceptibility接口仍需all。
K(q)是静态能量曲率，不能替换QMC等时S(q)或有限投影响应chi(q)。
零模保留为不能认证状态及能量扫描，不加展宽；软搜索缺点时显式标记。

专项验证：
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest discover -s quick1010/kitaev_ssh -p test_static_vbs_response.py -v

本次实际结果和复现：../fast-1006/audit/1010/static-vbs-note-20261010/REPORT.txt
定义和边界：docs/decisions/ADR-018-fixed-background-static-curvature.md
