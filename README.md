### Bigelow_Kim2023_CellBiol
https://doi.org/10.1016/j.cub.2023.01.016

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

If you want to reproduce the findings:
1) Download the “Bigelow2023_CB” folder which contains ‘analysisCode’ and ‘dataFiles’
2) Open MATLAB , then move to the ‘Bigelow2023_CB > analysisCode’ folder
3) Run make_FigureXX.m code

<u>Conditions for using the data</u>
If you publish any work using the data, please cite the publication above (Bigelow et al., 2023),
and also cite the data set using the following:

Bigelow, Anthony; Kim, Taekjun; Namima, Tomoyuki; Bair, Wyeth; Pasupathy, Anitha (2022),
“Dataset: Dissociation in neuronal encoding of object versus surface motion in the primate
brain”, Mendeley Data, V1, doi: 10.17632/cs76nk38zj.1
