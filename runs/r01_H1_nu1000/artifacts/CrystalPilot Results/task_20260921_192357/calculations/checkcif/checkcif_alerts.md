# checkCIF alerts - 2026-09-21 20:06

target: publication:CrystalPilot Results/task_20260921_192357/final.cif
counts: A=9  B=1  C=20  G=17
delta vs job_20260921_195857_748814_4c068ca1: new [], resolved []

## A alerts (must fix or justify)

### 020_ALERT_3_A  The Value of Rint is Greater Than 0.12 .........      0.566 Report
- meaning: Rint 明显偏大（与 RINTA01 同源：>0.10 C 级、>0.15 B 级、>0.20 A 级；实践中 >0.18 即 B 级）
- causes: 晶体质量/吸收校正不足/对称性选高了/坏的采集轮次混入
- remedy: 核对 Laue 组与吸收校正；多 run 数据按轮分组查 Rint，坏轮整体舍弃优于调 rejection 阈值；如实报告
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 023_ALERT_3_A  Resolution (too) Low [sin(theta)/Lambda < 0.6]..       0.50 Ang-1
- meaning: 数据分辨率过低 (sin θ/λ < 0.6)
- causes: 数据确实只收到低角；或 CIF 没带反射数据导致按 0 计算
- remedy: 如数据确实有限，说明采集条件（同步辐射短波长、弱衍射晶体）
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 183_ALERT_1_A  Missing _cell_measurement_reflns_used Value ....     Please Do !
- meaning: 缺 _cell_measurement_reflns_used
- causes: 未记录晶胞测定用的反射数（来自指标化软件）
- remedy: 从数据处理软件（SAINT/CrysAlis/DIALS）提取，或如实留空并说明
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 184_ALERT_1_A  Missing _cell_measurement_theta_min Value ......     Please Do !
- meaning: 缺 _cell_measurement_theta_min
- causes: 同 183
- remedy: 同 183
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 185_ALERT_1_A  Missing _cell_measurement_theta_max Value ......     Please Do !
- meaning: 缺 _cell_measurement_theta_max
- causes: 同 183
- remedy: 同 183
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 197_ALERT_1_A  Missing _cell_measurement_temperature Datum ....     Please Add
- meaning: 缺 _diffrn_ambient_temperature
- causes: 未记录采集温度
- remedy: 补 experiment.temperature_K
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 198_ALERT_1_A  Missing _diffrn_ambient_temperature   Datum ....     Please Add
- meaning: 缺 _cell_measurement_temperature
- causes: 同 197
- remedy: 同 197
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 660_ALERT_1_A  No Valid _diffrn_radiation_type Value Reported .     Please Do !
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 699_ALERT_1_A  Missing _exptl_crystal_description Value .......     Please Do !
- meaning: 缺晶体外观描述 (_exptl_crystal_description)
- causes: 未记录晶体形貌
- remedy: 补 experiment.crystal
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

## B alerts

### 196_ALERT_1_B  No TEMP record and _measurement_temperature .NE.        293 Degree
- meaning: 无 TEMP 记录且温度 ≠ 293 K
- causes: CIF 温度与内嵌 .res 的 TEMP 卡不一致/缺失
- remedy: 在 .ins 写 TEMP 卡（CrystalPilot 从 experiment.temperature_K 自动写入）
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

## C alerts

### 026_ALERT_3_C  Ratio Observed / Unique Reflections (too) Low ..        47% Check
- meaning: Check for a weak data set.
- remedy: 可以尝试用SHEL 9990.84命令或OMIT-2 50命令截去一部分高角度点，如果不行就只能重新收数据了。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 052_ALERT_1_C  Info on Absorption Correction Method   Not Given     Please Do !
- meaning: Test for specification absorption correction method [0,1].
- remedy: 在CIF文件中"_exptl_absorpt_correction_type"项给出吸收校正方法即可。若未做吸收校正，此项填写'none'或'?'，并将"_exptl_absorpt_correction_T_max"、"_exptl_absorpt_correction_T_min"和"_exptl_absorpt_process_details"项内容改为'?'即可。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 053_ALERT_1_C  Minimum Crystal Dimension Missing (or Error) ...     Please Check
- meaning: Test for specification xtal_dimension_min[0,1].
- remedy: 在_exptl_crystal_size_min项中给出晶体最小尺寸，单位为毫米。或者在ins中添加SIZE指令后精修。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 054_ALERT_1_C  Medium  Crystal Dimension Missing (or Error) ...     Please Check
- meaning: Test for specification xtal_dimension_mid [0,1].
- remedy: 在_exptl_crystal_size_mid项中给出晶体中间尺寸，单位为毫米。或者在ins中添加SIZE指令后精修。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 055_ALERT_1_C  Maximum Crystal Dimension Missing (or Error) ...     Please Check
- meaning: Test for specification xtal_dimension_max [0,1].
- remedy: 在_exptl_crystal_size_max项中给出晶体最大尺寸，单位为毫米。或者在ins中添加SIZE指令后精修。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 141_ALERT_4_C  s.u. on a - Axis Small or Missing ..............    0.00000 Ang.
- meaning: a 轴标准不确定度缺失或为 0
- causes: ZERR 未带真实晶胞 esd（粗解/中间模型常见）
- remedy: 从原始 .ins 传递 ZERR esd，或由指标化软件提供
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 143_ALERT_4_C  s.u. on c - Axis Small or Missing ..............    0.00000 Ang.
- meaning: c 轴标准不确定度缺失或为 0
- causes: 同 141
- remedy: 同 141
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 151_ALERT_1_C  No s.u. (esd) Given on Volume ..................     Please Do !
- meaning: 晶胞体积未给出 s.u.
- causes: 晶胞 esd 全零，体积 esd 无法传播
- remedy: 同 141
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 241_ALERT_2_C  High   'MainMol' Ueq as Compared to Neighbors of         O3 Check
- meaning: 原子 Ueq 明显高于相邻原子
- causes: 占有率<1 的位置按全占精修；轻元素被判成重元素
- remedy: 检查该位点元素/占有率是否有密度依据
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 241_ALERT_2_C  High   'MainMol' Ueq as Compared to Neighbors of        C11 Check
- meaning: 原子 Ueq 明显高于相邻原子
- causes: 占有率<1 的位置按全占精修；轻元素被判成重元素
- remedy: 检查该位点元素/占有率是否有密度依据
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 242_ALERT_2_C  Low    'MainMol' Ueq as Compared to Neighbors of        Zr1 Check
- meaning: 原子 Ueq 明显低于相邻原子
- causes: 重元素被判成轻元素；或该位点确为更重原子
- remedy: 同 241，反向核查
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 242_ALERT_2_C  Low    'MainMol' Ueq as Compared to Neighbors of        Zr2 Check
- meaning: 原子 Ueq 明显低于相邻原子
- causes: 重元素被判成轻元素；或该位点确为更重原子
- remedy: 同 241，反向核查
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 260_ALERT_2_C  Large Average Ueq of Residue Including       Zr1      0.124 Check
- meaning: 某原子 Ueq 异常大
- causes: 鬼原子/部分占有/无序
- remedy: 检查差值密度支持，必要时删除或设占有率
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 334_ALERT_2_C  Small <C-C> Benzene Dist.   C6       -C9       .       1.37 Ang.
- meaning: Check Average in Multiple Substituted Benzene Type C-C.
- remedy: 检查是否需要无序处理，用DFIX限制C23-C28距离后精修。 / 将此苯环加FLAT限制后精修。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 342_ALERT_3_C  Low Bond Precision on  C-C Bonds ...............    0.01514 Ang.
- meaning: Check Bond Precision for C-C in Structures (Z(max) > 39).
- remedy: 找出有问题的C-C键，限制一下试试。检查精修模型是否正确，是否有末端需要无序处理，晶胞参数误差是否过大。如果均无法解决可尝试收集低温数据。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 905_ALERT_3_C  Negative K value in the Analysis of Variance ...    -26.334 Report
- meaning: Report Negative K values in the Analysis of Variance.
- remedy: 一般不会出现，有合理解释即可。 / 使用Shel 999 0.84删除高角度的点，如无法解决有合理解释即可。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 905_ALERT_3_C  Negative K value in the Analysis of Variance ...     -4.461 Report
- meaning: Report Negative K values in the Analysis of Variance.
- remedy: 一般不会出现，有合理解释即可。 / 使用Shel 999 0.84删除高角度的点，如无法解决有合理解释即可。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 911_ALERT_3_C  Missing FCF Refl Between Thmin & STh/L=    0.500          2 Report
-4 34  3,  -1 26  7,
- meaning: θ(min) 与 sinθ/λ=0.6 之间的 FCF 反射缺失偏多
- causes: 完整度不足/挡板遮挡/扫描策略缺角
- remedy: 同 029/910；说明缺失来源
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 975_ALERT_2_C  Check Calcd Resid. Dens.  0.99Ang From O4      .       0.49 eA-3
- meaning: Test for positive density near N or O.
- remedy: Q峰游离检查是否有原子未指认，Q峰如果在轻原子周围检查是否可以无序处理，在重原子周围需要重新吸收校正或者收集数据。使用Shel 999 0.84切去部分高角度点可缓解此警告。最后检查是否数据足够好，是否存在孪晶。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 976_ALERT_2_C  Check Calcd Resid. Dens.  0.70Ang From O6      .      -0.56 eA-3
- meaning: Test for negative density near N or O.
- remedy: Q峰游离检查是否有原子未指认，Q峰如果在轻原子周围检查是否可以无序处理，在重原子周围需要重新吸收校正或者收集数据。使用Shel 999 0.84切去部分高角度点可缓解此警告。最后检查是否数据足够好，是否存在孪晶。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

## G notes (codes only)
004 x1, 009 x1, 042 x1, 045 x1, 092 x1, 606 x1, 794 x2, 802 x1, 869 x1, 910 x1, 913 x1, 969 x1, 978 x1, 984 x2, 985 x1
