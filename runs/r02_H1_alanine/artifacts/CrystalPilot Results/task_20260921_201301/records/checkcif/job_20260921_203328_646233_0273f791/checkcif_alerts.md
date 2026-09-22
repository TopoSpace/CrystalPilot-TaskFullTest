# checkCIF alerts - 2026-09-21 20:33

target: publication:CrystalPilot Results/task_20260921_201301/final.cif
counts: A=4  B=3  C=19  G=14

## A alerts (must fix or justify)

### 058_ALERT_1_A  Maximum Transmission Factor Missing ............          ?
- meaning: Test for specification Tmax [0,1].
- remedy: 查数据还原后的文件,比如sad.abs.从中寻找；或者在ins文件中加上命令 size a bc (abc是晶体的大小，测定时有记录在P4P文件中，有的老师可能不给你记录，你可以根据你的感观填上数值如size0.18 0.12 0.10)，然后refine一轮CIF中即有此数值。填入“_exptl_absorpt_correction_T_max”值即可。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 059_ALERT_1_A  Minimum Transmission Factor Missing ............          ?
- meaning: Test for specification Tmin [0,1].
- remedy: 查数据还原后的文件,比如sad.abs.从中寻找；或者在ins文件中加上命令 size a bc (abc是晶体的大小，测定时有记录在P4P文件中，有的老师可能不给你记录，你可以根据你的感观填上数值如size0.18 0.12 0.10)，然后refine一轮CIF中即有此数值。填入_exptl_absorpt_correction_T_min值即可。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 699_ALERT_1_A  Missing _exptl_crystal_description Value .......     Please Do !
- meaning: 缺晶体外观描述 (_exptl_crystal_description)
- causes: 未记录晶体形貌
- remedy: 补 experiment.crystal
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 939_ALERT_3_A  Large Value of Not (SHELXL) Weight Optimized S .     829.83 Check
- meaning: 疑存在大 |Error/esd| 离群反射
- causes: 个别反射与模型严重失配（挡板遮挡/冰环/孪晶重叠）
- remedy: Olex2 Bad reflections（或 lst 表）核查 |Err/esd|>10 的点：确有仪器成因可 OMIT 并记录；成片失配则查数据侧
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

## B alerts

### 035_ALERT_1_B  _chemical_absolute_configuration Info  Not Given     Please Do !
- meaning: Check for_chemical_absolute_configuration.
- remedy: 根据需要选择绝对构型的确定方式，非手性空间群改成句点。手性空间群此项后根据实际情况改成上述关键词中的一项。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 097_ALERT_2_B  Large Reported Max.  (Positive) Residual Density       0.89 eA-3
- meaning: 最高正残差峰较高（Max positive residual density）
- causes: 未建模原子/无序/吸收；高分辨数据的高角噪声也会抬峰（专家案例：0.48 Å 数据 Max 1.0 报 B，SHEL 截 0.70 Å 降至 0.80 转 C）
- remedy: 给出峰位置与最近原子；无化学意义时按小步（~0.05 Å）截断分辨率观察峰高单调下降，绝不为吸收残峰加假原子
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 934_ALERT_3_B  Number of (Iobs-Icalc)/Sigma(W) > 10 Outliers ..          6 Check
-1  3  1,   1  3  1,   1  0  2,   0  2  2,   0  2  3,   0  2  4,
- meaning: 报告的与计算的 Rint 不一致
- causes: 数据经过预合并或掩膜贡献
- remedy: 说明数据处理链
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

## C alerts

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

### 098_ALERT_2_C  Large Reported Min.  (Negative) Residual Density      -0.76 eA-3
- meaning: 最深负残差峰较深（Min negative residual density）
- causes: 吸收校正不足/重原子附近截断效应/占有率或元素指认过重
- remedy: 给出谷位置与最近原子，说明成因（重原子 <1 Å 的涟漪属正常并说明）
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 141_ALERT_4_C  s.u. on a - Axis Small or Missing ..............    0.00000 Ang.
- meaning: a 轴标准不确定度缺失或为 0
- causes: ZERR 未带真实晶胞 esd（粗解/中间模型常见）
- remedy: 从原始 .ins 传递 ZERR esd，或由指标化软件提供
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 142_ALERT_4_C  s.u. on b - Axis Small or Missing ..............    0.00000 Ang.
- meaning: b 轴标准不确定度缺失或为 0
- causes: 同 141
- remedy: 同 141
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

### 250_ALERT_2_C  Large U3/U1 Ratio for <U(i,j)> Tensor(Resd    1)        2.1 Note
- meaning: 平均 ADP 张量 U3/U1 比大（整体各向异性强）
- causes: 层状/链状晶体真实热振动各向异性；数据各向异性截断
- remedy: 结构合理时如实说明
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 340_ALERT_3_C  Low Bond Precision on  C-C Bonds ...............     0.0045 Ang.
- meaning: 键长 s.u. 偏大
- causes: 数据/参数比低或弱数据
- remedy: 说明数据质量限制
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 907_ALERT_2_C  Flack x > 0.5, Structure Needs to be Inverted? .       0.80 Check
- meaning: Check whether the structure needs to be inverted.
- remedy: 对结构进行翻转即可。(Olex2中可通过键入inv -f指令实现构型翻转) / 无需解决，有合理解释即可。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 915_ALERT_3_C  No Flack x Check Done: Low Friedel Pair Coverage         78 %
- meaning: Test for low Friedel Pair Coverage in non-centro structure.
- remedy: 一般不会出现，有合理解释即可。 / 数据收集不当，与精修没关系，估计是无心当有心的收了；不传fcf就没事了，如果不定绝对构型有这个错误无所谓。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 918_ALERT_3_C  Reflection(s) with I(obs) much Smaller I(calc) .          5 Check
- meaning: Test for reflections with I(obs) << I(calc).
- remedy: OMIT这些点即可。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 976_ALERT_2_C  Check Calcd Resid. Dens.  0.96Ang From O1      .      -0.79 eA-3
- meaning: Test for negative density near N or O.
- remedy: Q峰游离检查是否有原子未指认，Q峰如果在轻原子周围检查是否可以无序处理，在重原子周围需要重新吸收校正或者收集数据。使用Shel 999 0.84切去部分高角度点可缓解此警告。最后检查是否数据足够好，是否存在孪晶。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 976_ALERT_2_C  Check Calcd Resid. Dens.  0.99Ang From N1      .      -0.53 eA-3
- meaning: Test for negative density near N or O.
- remedy: Q峰游离检查是否有原子未指认，Q峰如果在轻原子周围检查是否可以无序处理，在重原子周围需要重新吸收校正或者收集数据。使用Shel 999 0.84切去部分高角度点可缓解此警告。最后检查是否数据足够好，是否存在孪晶。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 976_ALERT_2_C  Check Calcd Resid. Dens.  0.97Ang From N1      .      -0.48 eA-3
- meaning: Test for negative density near N or O.
- remedy: Q峰游离检查是否有原子未指认，Q峰如果在轻原子周围检查是否可以无序处理，在重原子周围需要重新吸收校正或者收集数据。使用Shel 999 0.84切去部分高角度点可缓解此警告。最后检查是否数据足够好，是否存在孪晶。
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 977_ALERT_2_C  Check Negative Difference Density on H1C       .      -0.47 eA-3
- meaning: reflns_number 与 hkl 数不一致
- causes: 预合并/掩膜
- remedy: 说明数据链
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 977_ALERT_2_C  Check Negative Difference Density on H3A       .      -0.31 eA-3
- meaning: reflns_number 与 hkl 数不一致
- causes: 预合并/掩膜
- remedy: 说明数据链
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

### 977_ALERT_2_C  Check Negative Difference Density on H3B       .      -0.35 eA-3
- meaning: reflns_number 与 hkl 数不一致
- causes: 预合并/掩膜
- remedy: 说明数据链
- explanation: (fill in: fixed, or the evidence-based justification for VALIDATION.md)

## G notes (codes only)
007 x1, 032 x1, 484 x3, 795 x1, 802 x1, 883 x1, 912 x1, 916 x1, 955 x1, 965 x1, 969 x1, 978 x1
