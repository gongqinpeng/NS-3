from ns import *
import networkx as nx

# 生成一个简单的图和MST
G = nx.Graph()
G.add_edge('A', 'B', weight=1)
G.add_edge('B', 'C', weight=2)
G.add_edge('A', 'C', weight=3)

# 使用Prim算法生成MST
MST = nx.minimum_spanning_tree(G)

# 在NS-3中创建节点
nodes = ns.network.NodeContainer()
nodes.Create(len(MST.nodes))

# 创建一个字典来映射网络节点到NS-3节点
node_mapping = {}
for i, node in enumerate(MST.nodes):
    node_mapping[node] = nodes.Get(i)

# 配置点对点链路
pointToPoint = ns.point_to_point.PointToPointHelper()

# 遍历MST的边来创建点对点链接
for (u, v) in MST.edges():
    u_node = node_mapping[u]
    v_node = node_mapping[v]
    link = ns.network.NodeContainer(u_node, v_node)
    devices = pointToPoint.Install(link)

# 安装互联网堆栈
internet = ns.internet.InternetStackHelper()
internet.Install(nodes)

# 分配IP地址（这里只是一个示例，实际需要根据网络结构进行配置）
address = ns.internet.Ipv4AddressHelper()
for (u, v) in MST.edges():
    u_node = node_mapping[u]
    v_node = node_mapping[v]
    interface = ns.network.NetDeviceContainer()
    interface.Add(u_node.GetDevice(0))
    interface.Add(v_node.GetDevice(0))
    address.Assign(interface)

# 启用IPV4路由
ns.internet.Ipv4GlobalRoutingHelper.PopulateRoutingTables()

# 设置仿真时间并运行
ns.core.Simulator.Stop(ns.core.Seconds(10.0))
ns.core.Simulator.Run()
ns.core.Simulator.Destroy()
