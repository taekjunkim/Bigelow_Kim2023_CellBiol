# Bigelow_Kim2023_CellBiol
This repository contains data, codes to reproduce the main results of the paper
"Dissociation in neuronal encoding of object versus surface motion in the primate brain". 

For example, cellData_sua.mat file contains spike count information of 115 single units across
stimulus conditions and trials. To read the responses of unit# 86,

```
>> cellData_sua(86)
ans =
struct with fields:
file_id: 'p190516'
unit_id: 1
rf_pix: [-72 -286]
respMtx: [7×41 double]
condition: {6×1 cell}
anova: [2×5 double]
dir_tuning: [8×5 double]
DI_base: [0.9167 0.5263 1 0.8889 0.5333]
```
