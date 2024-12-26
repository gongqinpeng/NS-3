import networkx as nx
import matplotlib.pyplot as plt
import json

# 读取JSON文件
with open('DOA-50.json', 'r') as file:
    topology = json.load(file)

# 创建图
G = nx.Graph()
for i in range(len(topology)):
    for j in range(len(topology[i])):
        if topology[i][j]:
            G.add_edge(i, j)

# 绘制图
nx.draw(G, with_labels=True)
plt.show()
