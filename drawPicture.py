import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# 数据
data = {
"MST":[89.464,	80.46	,55.161,	38.791,	16.4654
],
"PSO":[129.464	,114.56	,80.464,	58.464,	26.464
],
"DOA":[130.526	,115.46,	82.461,	57.49646,	27.1654
],
   "DOA - SIM":[130.859	,119.26,	85.79,	59.492,	29.495],
#     "FastDDS":[16.6	,16.4	,32.8	,31.9,	32.7	,40.2,	38.4	,34.5,	37.8
#
#
#
# ],
#     "CycloneDDS":[66.4,	65.9,	66.7,	66.2	,67.1,	70.2,	85.4	,98.5,	104.3
#
#
#
# ],
#     "ChangeDDS":[17.2,	16.9,30.5,	31.4,	32.8,	38.4,	38.5,	37.5	,38.4
#
#
#
# ]
}

# x轴坐标
x = [10,	20,	50,	80	,100

]

# 设置中文字体为SimHei（黑体）
font = FontProperties(fname=r"c:\windows\fonts\simhei.ttf", size=14)

# 绘制折线图
for method, values in data.items():
    plt.plot(x, values, marker='o', label=method)

# 设置标题和坐标轴标签
plt.title('拓扑持续时间', fontproperties=font)
plt.xlabel('节点数量', fontproperties=font)
plt.ylabel('持续时间(s)', fontproperties=font)

# 添加图例
plt.legend(prop=font)

# 显示网格
plt.grid(True)
plt.tight_layout()
# 显示图形
plt.show()
##plt.savefig("C:\\Users\\Administrator\\Desktop\\毕设\\研究生毕设\\研究生毕设 (2)\\研究生毕设\\python图片\\网络吞吐量对比图.jpg")