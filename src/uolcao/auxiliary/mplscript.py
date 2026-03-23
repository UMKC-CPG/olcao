#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv('KO_Real.dat.A', sep='\\s+')

# Arrays to hold figures, subplots, axes, and curves.
figs = []
subplots = []
curves = []
figs.append(plt.figure(1, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_s'], label='s_s', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_s'], label='x_s', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_s'], label='y_s', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_s'], label='z_s', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_s'], label='xx_s', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_s'], label='yy_s', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_s'], label='zz_s', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_s'], label='xy_s', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_s'], label='xz_s', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_s'], label='yz_s', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_s'], label='xyz_s', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_s'], label='xxy_s', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_s'], label='xxz_s', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_s'], label='yyx_s', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_s'], label='yyz_s', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_s'], label='zzx_s', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_s'], label='zzy_s', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_s'], label='xxx_s', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_s'], label='yyy_s', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_s'], label='zzz_s', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_x'], label='s_x', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_x'], label='x_x', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_x'], label='y_x', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_x'], label='z_x', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_x'], label='xx_x', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_x'], label='yy_x', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_x'], label='zz_x', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_x'], label='xy_x', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_x'], label='xz_x', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_x'], label='yz_x', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_x'], label='xyz_x', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_x'], label='xxy_x', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_x'], label='xxz_x', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_x'], label='yyx_x', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_x'], label='yyz_x', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_x'], label='zzx_x', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_x'], label='zzy_x', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_x'], label='xxx_x', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_x'], label='yyy_x', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_x'], label='zzz_x', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(2, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_y'], label='s_y', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_y'], label='x_y', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_y'], label='y_y', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_y'], label='z_y', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_y'], label='xx_y', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_y'], label='yy_y', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_y'], label='zz_y', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_y'], label='xy_y', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_y'], label='xz_y', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_y'], label='yz_y', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_y'], label='xyz_y', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_y'], label='xxy_y', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_y'], label='xxz_y', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_y'], label='yyx_y', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_y'], label='yyz_y', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_y'], label='zzx_y', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_y'], label='zzy_y', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_y'], label='xxx_y', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_y'], label='yyy_y', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_y'], label='zzz_y', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_z'], label='s_z', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_z'], label='x_z', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_z'], label='y_z', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_z'], label='z_z', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_z'], label='xx_z', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_z'], label='yy_z', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_z'], label='zz_z', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_z'], label='xy_z', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_z'], label='xz_z', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_z'], label='yz_z', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_z'], label='xyz_z', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_z'], label='xxy_z', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_z'], label='xxz_z', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_z'], label='yyx_z', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_z'], label='yyz_z', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_z'], label='zzx_z', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_z'], label='zzy_z', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_z'], label='xxx_z', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_z'], label='yyy_z', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_z'], label='zzz_z', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(3, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_xx'], label='s_xx', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_xx'], label='x_xx', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_xx'], label='y_xx', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_xx'], label='z_xx', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_xx'], label='xx_xx', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_xx'], label='yy_xx', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_xx'], label='zz_xx', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_xx'], label='xy_xx', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_xx'], label='xz_xx', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_xx'], label='yz_xx', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_xx'], label='xyz_xx', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_xx'], label='xxy_xx', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_xx'], label='xxz_xx', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_xx'], label='yyx_xx', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_xx'], label='yyz_xx', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_xx'], label='zzx_xx', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_xx'], label='zzy_xx', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_xx'], label='xxx_xx', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_xx'], label='yyy_xx', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_xx'], label='zzz_xx', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_yy'], label='s_yy', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_yy'], label='x_yy', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_yy'], label='y_yy', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_yy'], label='z_yy', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_yy'], label='xx_yy', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_yy'], label='yy_yy', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_yy'], label='zz_yy', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_yy'], label='xy_yy', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_yy'], label='xz_yy', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_yy'], label='yz_yy', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_yy'], label='xyz_yy', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_yy'], label='xxy_yy', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_yy'], label='xxz_yy', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_yy'], label='yyx_yy', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_yy'], label='yyz_yy', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_yy'], label='zzx_yy', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_yy'], label='zzy_yy', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_yy'], label='xxx_yy', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_yy'], label='yyy_yy', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_yy'], label='zzz_yy', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(4, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_zz'], label='s_zz', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_zz'], label='x_zz', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_zz'], label='y_zz', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_zz'], label='z_zz', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_zz'], label='xx_zz', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_zz'], label='yy_zz', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_zz'], label='zz_zz', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_zz'], label='xy_zz', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_zz'], label='xz_zz', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_zz'], label='yz_zz', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_zz'], label='xyz_zz', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_zz'], label='xxy_zz', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_zz'], label='xxz_zz', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_zz'], label='yyx_zz', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_zz'], label='yyz_zz', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_zz'], label='zzx_zz', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_zz'], label='zzy_zz', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_zz'], label='xxx_zz', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_zz'], label='yyy_zz', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_zz'], label='zzz_zz', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_xy'], label='s_xy', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_xy'], label='x_xy', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_xy'], label='y_xy', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_xy'], label='z_xy', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_xy'], label='xx_xy', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_xy'], label='yy_xy', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_xy'], label='zz_xy', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_xy'], label='xy_xy', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_xy'], label='xz_xy', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_xy'], label='yz_xy', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_xy'], label='xyz_xy', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_xy'], label='xxy_xy', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_xy'], label='xxz_xy', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_xy'], label='yyx_xy', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_xy'], label='yyz_xy', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_xy'], label='zzx_xy', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_xy'], label='zzy_xy', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_xy'], label='xxx_xy', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_xy'], label='yyy_xy', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_xy'], label='zzz_xy', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(5, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_xz'], label='s_xz', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_xz'], label='x_xz', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_xz'], label='y_xz', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_xz'], label='z_xz', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_xz'], label='xx_xz', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_xz'], label='yy_xz', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_xz'], label='zz_xz', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_xz'], label='xy_xz', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_xz'], label='xz_xz', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_xz'], label='yz_xz', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_xz'], label='xyz_xz', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_xz'], label='xxy_xz', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_xz'], label='xxz_xz', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_xz'], label='yyx_xz', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_xz'], label='yyz_xz', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_xz'], label='zzx_xz', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_xz'], label='zzy_xz', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_xz'], label='xxx_xz', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_xz'], label='yyy_xz', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_xz'], label='zzz_xz', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_yz'], label='s_yz', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_yz'], label='x_yz', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_yz'], label='y_yz', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_yz'], label='z_yz', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_yz'], label='xx_yz', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_yz'], label='yy_yz', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_yz'], label='zz_yz', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_yz'], label='xy_yz', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_yz'], label='xz_yz', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_yz'], label='yz_yz', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_yz'], label='xyz_yz', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_yz'], label='xxy_yz', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_yz'], label='xxz_yz', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_yz'], label='yyx_yz', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_yz'], label='yyz_yz', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_yz'], label='zzx_yz', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_yz'], label='zzy_yz', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_yz'], label='xxx_yz', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_yz'], label='yyy_yz', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_yz'], label='zzz_yz', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(6, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_xyz'], label='s_xyz', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_xyz'], label='x_xyz', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_xyz'], label='y_xyz', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_xyz'], label='z_xyz', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_xyz'], label='xx_xyz', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_xyz'], label='yy_xyz', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_xyz'], label='zz_xyz', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_xyz'], label='xy_xyz', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_xyz'], label='xz_xyz', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_xyz'], label='yz_xyz', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_xyz'], label='xyz_xyz', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_xyz'], label='xxy_xyz', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_xyz'], label='xxz_xyz', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_xyz'], label='yyx_xyz', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_xyz'], label='yyz_xyz', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_xyz'], label='zzx_xyz', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_xyz'], label='zzy_xyz', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_xyz'], label='xxx_xyz', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_xyz'], label='yyy_xyz', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_xyz'], label='zzz_xyz', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_xxy'], label='s_xxy', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_xxy'], label='x_xxy', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_xxy'], label='y_xxy', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_xxy'], label='z_xxy', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_xxy'], label='xx_xxy', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_xxy'], label='yy_xxy', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_xxy'], label='zz_xxy', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_xxy'], label='xy_xxy', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_xxy'], label='xz_xxy', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_xxy'], label='yz_xxy', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_xxy'], label='xyz_xxy', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_xxy'], label='xxy_xxy', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_xxy'], label='xxz_xxy', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_xxy'], label='yyx_xxy', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_xxy'], label='yyz_xxy', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_xxy'], label='zzx_xxy', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_xxy'], label='zzy_xxy', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_xxy'], label='xxx_xxy', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_xxy'], label='yyy_xxy', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_xxy'], label='zzz_xxy', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(7, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_xxz'], label='s_xxz', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_xxz'], label='x_xxz', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_xxz'], label='y_xxz', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_xxz'], label='z_xxz', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_xxz'], label='xx_xxz', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_xxz'], label='yy_xxz', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_xxz'], label='zz_xxz', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_xxz'], label='xy_xxz', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_xxz'], label='xz_xxz', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_xxz'], label='yz_xxz', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_xxz'], label='xyz_xxz', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_xxz'], label='xxy_xxz', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_xxz'], label='xxz_xxz', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_xxz'], label='yyx_xxz', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_xxz'], label='yyz_xxz', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_xxz'], label='zzx_xxz', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_xxz'], label='zzy_xxz', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_xxz'], label='xxx_xxz', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_xxz'], label='yyy_xxz', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_xxz'], label='zzz_xxz', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_yyx'], label='s_yyx', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_yyx'], label='x_yyx', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_yyx'], label='y_yyx', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_yyx'], label='z_yyx', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_yyx'], label='xx_yyx', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_yyx'], label='yy_yyx', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_yyx'], label='zz_yyx', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_yyx'], label='xy_yyx', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_yyx'], label='xz_yyx', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_yyx'], label='yz_yyx', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_yyx'], label='xyz_yyx', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_yyx'], label='xxy_yyx', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_yyx'], label='xxz_yyx', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_yyx'], label='yyx_yyx', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_yyx'], label='yyz_yyx', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_yyx'], label='zzx_yyx', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_yyx'], label='zzy_yyx', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_yyx'], label='xxx_yyx', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_yyx'], label='yyy_yyx', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_yyx'], label='zzz_yyx', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(8, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_yyz'], label='s_yyz', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_yyz'], label='x_yyz', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_yyz'], label='y_yyz', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_yyz'], label='z_yyz', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_yyz'], label='xx_yyz', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_yyz'], label='yy_yyz', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_yyz'], label='zz_yyz', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_yyz'], label='xy_yyz', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_yyz'], label='xz_yyz', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_yyz'], label='yz_yyz', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_yyz'], label='xyz_yyz', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_yyz'], label='xxy_yyz', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_yyz'], label='xxz_yyz', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_yyz'], label='yyx_yyz', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_yyz'], label='yyz_yyz', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_yyz'], label='zzx_yyz', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_yyz'], label='zzy_yyz', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_yyz'], label='xxx_yyz', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_yyz'], label='yyy_yyz', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_yyz'], label='zzz_yyz', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_zzx'], label='s_zzx', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_zzx'], label='x_zzx', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_zzx'], label='y_zzx', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_zzx'], label='z_zzx', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_zzx'], label='xx_zzx', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_zzx'], label='yy_zzx', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_zzx'], label='zz_zzx', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_zzx'], label='xy_zzx', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_zzx'], label='xz_zzx', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_zzx'], label='yz_zzx', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_zzx'], label='xyz_zzx', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_zzx'], label='xxy_zzx', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_zzx'], label='xxz_zzx', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_zzx'], label='yyx_zzx', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_zzx'], label='yyz_zzx', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_zzx'], label='zzx_zzx', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_zzx'], label='zzy_zzx', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_zzx'], label='xxx_zzx', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_zzx'], label='yyy_zzx', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_zzx'], label='zzz_zzx', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(9, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_zzy'], label='s_zzy', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_zzy'], label='x_zzy', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_zzy'], label='y_zzy', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_zzy'], label='z_zzy', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_zzy'], label='xx_zzy', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_zzy'], label='yy_zzy', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_zzy'], label='zz_zzy', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_zzy'], label='xy_zzy', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_zzy'], label='xz_zzy', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_zzy'], label='yz_zzy', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_zzy'], label='xyz_zzy', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_zzy'], label='xxy_zzy', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_zzy'], label='xxz_zzy', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_zzy'], label='yyx_zzy', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_zzy'], label='yyz_zzy', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_zzy'], label='zzx_zzy', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_zzy'], label='zzy_zzy', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_zzy'], label='xxx_zzy', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_zzy'], label='yyy_zzy', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_zzy'], label='zzz_zzy', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_xxx'], label='s_xxx', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_xxx'], label='x_xxx', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_xxx'], label='y_xxx', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_xxx'], label='z_xxx', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_xxx'], label='xx_xxx', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_xxx'], label='yy_xxx', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_xxx'], label='zz_xxx', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_xxx'], label='xy_xxx', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_xxx'], label='xz_xxx', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_xxx'], label='yz_xxx', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_xxx'], label='xyz_xxx', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_xxx'], label='xxy_xxx', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_xxx'], label='xxz_xxx', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_xxx'], label='yyx_xxx', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_xxx'], label='yyz_xxx', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_xxx'], label='zzx_xxx', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_xxx'], label='zzy_xxx', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_xxx'], label='xxx_xxx', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_xxx'], label='yyy_xxx', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_xxx'], label='zzz_xxx', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
figs.append(plt.figure(10, figsize=[6.8, 4.4]))
subplots.append(plt.subplot(1, 1, 1))
data_temp = data.iloc[:, [0]]
x_min = data_temp.min()
if (type(x_min) is not float and type(x_min) is not int):
    x_min = x_min.min()
if (x_min != 0):
    x_min = (int(x_min/5) - 1)  * 5
x_max = data_temp.max()
if (type(x_max) is not float and type(x_max) is not int):
    x_max = x_max.max()
if (x_max != 0):
    x_max = (int(x_max/5) + 1) * 5
plt.xlim((x_min, x_max))
data_temp = data.iloc[:, [378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417]]
y_min = data_temp.min()
if (type(y_min) is not float and type(y_min) is not int):
    y_min = y_min.min()
if(y_min != 0):
    y_min = (int(y_min/1e-14) - 1)  * 1e-14
y_max = data_temp.max()
if (type(y_max) is not float and type(y_max) is not int):
    y_max = y_max.max()
if(y_max != 0):
    y_max = (int(y_max/1e-14) + 1) * 1e-14
plt.ylim((y_min, y_max))
curves.append(plt.plot(data['IDX'], data['s_yyy'], label='s_yyy', color='xkcd:strong blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_yyy'], label='x_yyy', color='xkcd:toxic green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_yyy'], label='y_yyy', color='xkcd:windows blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_yyy'], label='z_yyy', color='xkcd:blue blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_yyy'], label='xx_yyy', color='xkcd:blue with a hint of purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_yyy'], label='yy_yyy', color='xkcd:booger', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_yyy'], label='zz_yyy', color='xkcd:bright sea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_yyy'], label='xy_yyy', color='xkcd:dark green blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_yyy'], label='xz_yyy', color='xkcd:deep turquoise', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_yyy'], label='yz_yyy', color='xkcd:green teal', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_yyy'], label='xyz_yyy', color='xkcd:strong pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_yyy'], label='xxy_yyy', color='xkcd:bland', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_yyy'], label='xxz_yyy', color='xkcd:deep aqua', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_yyy'], label='yyx_yyy', color='xkcd:lavender pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_yyy'], label='yyz_yyy', color='xkcd:light moss green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_yyy'], label='zzx_yyy', color='xkcd:light seafoam green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_yyy'], label='zzy_yyy', color='xkcd:olive yellow', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_yyy'], label='xxx_yyy', color='xkcd:pig pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_yyy'], label='yyy_yyy', color='xkcd:deep lilac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_yyy'], label='zzz_yyy', color='xkcd:desert', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['s_zzz'], label='s_zzz', color='xkcd:dusty lavender', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['x_zzz'], label='x_zzz', color='xkcd:purpley grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['y_zzz'], label='y_zzz', color='xkcd:purply', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['z_zzz'], label='z_zzz', color='xkcd:candy pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xx_zzz'], label='xx_zzz', color='xkcd:light pastel green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yy_zzz'], label='yy_zzz', color='xkcd:boring green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zz_zzz'], label='zz_zzz', color='xkcd:kiwi green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xy_zzz'], label='xy_zzz', color='xkcd:light grey green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xz_zzz'], label='xz_zzz', color='xkcd:orange pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yz_zzz'], label='yz_zzz', color='xkcd:tea green', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xyz_zzz'], label='xyz_zzz', color='xkcd:very light brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxy_zzz'], label='xxy_zzz', color='xkcd:egg shell', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxz_zzz'], label='xxz_zzz', color='xkcd:eggplant purple', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyx_zzz'], label='yyx_zzz', color='xkcd:powder pink', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyz_zzz'], label='yyz_zzz', color='xkcd:reddish grey', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzx_zzz'], label='zzx_zzz', color='xkcd:baby shit brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzy_zzz'], label='zzy_zzz', color='xkcd:liliac', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['xxx_zzz'], label='xxx_zzz', color='xkcd:stormy blue', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['yyy_zzz'], label='yyy_zzz', color='xkcd:ugly brown', marker='None', linestyle='solid', linewidth=1))
curves.append(plt.plot(data['IDX'], data['zzz_zzz'], label='zzz_zzz', color='xkcd:custard', marker='None', linestyle='solid', linewidth=1))
plt.legend()
plt.show()
