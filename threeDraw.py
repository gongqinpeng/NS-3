import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# 示例数据结构，你需要替换为真实的数据
data = {
    "AODV": {
        "速度10m/s": [9.84,	15.689	,20.789	,25.928,	30.485

],
        "速度20m/s": [11.952,	19.892,	25.19,	29.792,	36.419

]
    },
    "OLSR": {
        "速度10m/s": [24.65,	31.275	,36.849	,38.4052,	40.156

],
        "速度20m/s": [29.0862,	35.862	,40.92,	43.0781,	49.0592

]
    },

    "BATMAN - ADV": {
        "速度10m/s": [14.892,	19.485	,24.5872,	27.058,	29.429

],
        "速度20m/s": [19.0459,	23.07981,	26.813,	30.0792	,34.0494

]
    },
    "BATMAN - LINKSTATUS": {
        "速度10m/s": [10.5952,	13.482,	19.562	,21.985,	24.592

],
        "速度20m/s": [14.92	,18.0921,	22.892,	25.0495,	28.095

]
    }
}

# x轴坐标（节点数量）
x = [10, 20, 50, 80, 100]

# 定义不同节点速度对应的线条样式和标记，可根据喜好调整
styles = {
    "速度10m/s": {'linestyle': '-', 'marker': 'o'},
    "速度20m/s": {'linestyle': '--', 'marker': 's'}
}

# 设置中文字体为SimHei（黑体），解决中文显示问题
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

for method, speed_data in data.items():
    for speed, values in speed_data.items():
        plt.plot(x, values, label=f"{method} - {speed}", **styles[speed])

# 设置标题和坐标轴标签
plt.title('路由开销率对比图')
plt.xlabel('节点数量')
plt.ylabel('路由开销率(%)')

# 添加图例
handles, labels = plt.gca().get_legend_handles_labels()
# 重新构建图例，确保顺序正确且显示完整
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='best')

# 显示网格
plt.grid(True)
plt.tight_layout()
# 显示图形
plt.show()
